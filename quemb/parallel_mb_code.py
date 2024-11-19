import numpy as np
import re, pyscf, os, sys
from pyscf.lib import logger
from molbe.mbe import BE
from molbe.fragment import fragpart
from molbe.helper import get_scfObj, get_eri, get_veff
from molbe.solver import be_func
from pyscf import qmmm
from molbe import be_var
from multiprocessing import Pool

"""
Prepare Calculation Setup
"""

geom_file = sys.argv[1]  # Geometry File
big_basis_region = sys.argv[2]  # How large of an area has a big basis?
small_basis_inp = "sto-3g"  # Small basis set of the full system
big_basis_be_inp = sys.argv[3]  # BE level in big basis region
big_basis_inp = sys.argv[4]  # Basis set of big basis region
nprocs = int(sys.argv[5])  # Number of available processors
charge_inp = int(sys.argv[6])  # Charge of system
chempotopt = True if sys.argv[7] == "True" else False # Chemical Potential Optimization
if len(sys.argv) > 8:
    mm_charges = np.load(sys.argv[8])  # npy file with a list of MM charges (optional)
    mm_coords = np.load(sys.argv[9])  # npy file with a list of MM coordinates (optional)
else:
    mm_charges = None
    mm_coords = None

#Verbosity to control print levels
#Verbose levels: Quiet (0), Error (1), Warn (2), Note (3), Info (4), Debug (5)
verbose = 0

# QuEmb molbe be_var settings to control scratch
be_var.SCRATCH = "/scratch/lweisbur/"
be_var.CREATE_SCRATCH_DIR = True

# add numbers to all geometries
num_geom_file = geom_file[:-4] + "-numbered.xyz"

"""
Establish all functions
"""
log = logger.Logger(sys.stdout, verbose)

def mixed_basis_solver(
    idx,
    f,
    small_basis,
    big_basis,
    auto_frag,
    auto_cen,
    auto_hlist,
    auto_add_centers,
    big_basis_be,
    heavy_atoms,
    heavy_atoms_ind,
    add_center,
    geom,
    numbered_geom_file,
    charge,
    mm_charges,
    mm_coords,
    chempot,
    skip_unnecessary_frag_build=True,
):
    """Run Mixed-Basis Solver with QuEmb

    Parameters
    _________
    idx: int
        Fragment number corresponding to the big basis set region
    f: list
        List of heavy atom indices in fragment idx
    small_basis: str
        Small basis set outside the fragment region, typically 'STO-3G'
    big_basis: str
        Big basis set for the big basis region near the framgnet
    auto_frag: list
        List of lists, from QuEmb's autofrag, giving all fragments. Each 
        element contains the indices of the heavy atoms in each gragment
    auto_cen: list
        List of lists, from QuEmb's autofrag, giving the indices of the centers
        of each fragment
    auto_hlist: list
        List of lists, from QuEmb's autofrag, giving all hydrogen indices in 
        each fragment
    auto_add_centers: list
        List of lists, from QuEmb's autofrag, giving heavy atoms in each fragment
        which are not centers in any other fragment
    big_basis_be: str
        Size of big basis set region, options "BE1", "BE2", "BE3"
    heavy_atoms: list
        List of heavy atoms with numbers
    heavy_atoms_ind: list
        List of indices for all heavy atoms
    add_center: list
        List of indices for heavy atoms which are not "centers" in auto_cen
    geom: list
        List giving the geometry of the system
    numbered_geom_file: str
        String with the path to the numbered geometry file
    charge: int
        Charge of the system
    mm_charges: list
        List of MM point charges, ordinarily set to None
    mm_coords: list
        List of MM coordinates, ordinarily set to None
    chempot: list
        Chemical potential for the system, set for chemical potential matching
    skip_unnecessary_frag_build: bool, optional
        Do not build extra fragments, default True
    """

    log.note("Fragment number: %s", idx)
    basis_library = {}

    # Build basis for the full system in the small basis
    for a in all_atoms_list:
        basis_library[a] = small_basis

    # Build mixed basis library for all heavy atoms and hydrogens in fragment
    frag_atoms = []
    for ind in f:
        atom = str(geom[ind].split()[0])
        frag_atoms.append(atom)
        basis_library[atom] = big_basis
        for x in auto_hlist[ind]:
            atom = str(geom[x].split()[0])
            frag_atoms.append(atom)
            basis_library[atom] = big_basis

    log.info("Basis library for fragment %s: %s", idx, basis_library)
    # Identify centers of fragment to perform be1 on
    center = str(geom[auto_cen[idx]].split()[0])

    # Add all heavy atoms which are only in this big basis fragment
    centers = [center]
    for i in auto_add_centers[idx]:
        centers.append(str(geom[i].split()[0]))

    log.info("All centers for fragment %s: %s", idx, centers)
    log.debug("Heavy atoms in fragment %s: %s", idx, heavy_atoms)
    log.debug("Heavy atom indices: %s", heavy_atoms_ind)

    log.info("Atoms in fragment %s: %s", idx, frag_atoms)

    # Build HF for the full mixed-basis sysetm
    mol_mixed = pyscf.gto.Mole()
    mol_mixed.atom = numbered_geom_file
    mol_mixed.basis = basis_library
    mol_mixed.charge = charge
    mol_mixed.build()

    # Run HF for fragment
    log.note("Running HF...")
    if not (mm_charges is None) and not (mm_coords is None):
        log.note("with MM charges...")
        mf_mixed = qmmm.mm_charge(pyscf.scf.RHF(mol_mixed), mm_coords, mm_charges)
    else:
        mf_mixed = pyscf.scf.RHF(mol_mixed)
    mf_mixed.kernel()

    # Run BE(1) fragmentation for mixed basis system
    if verbose < 4:
        frag_mixed = fragpart(be_type=big_basis_be, mol=mol_mixed, print_frags=False)
    else:
        frag_mixed = fragpart(be_type=big_basis_be, mol=mol_mixed, print_frags=True)

    frag_auto_f = frag_mixed.Frag_atom 
    frag_auto_c = frag_mixed.center_atom
    frag_auto_h = frag_mixed.hlist_atom
    frag_add_centers = frag_mixed.add_center_atom

    center_atoms = [geom[i].split()[0] for i in frag_auto_c]
    frag_nums = []; center_atoms_retry = [] # for linear chains, this is fine (making sure center is included)
    # center_atoms_retry is used to ensure that the center is included in the
    # fragment list if it is not an "additional" center
    # for non-linear chains, there may be edge cases. This is a temporary fix.

    for i in centers:
        try:
            frag_nums.append(center_atoms.index(i))
        except:
            center_atoms_retry.append(i)
    
    for i in center_atoms_retry:
        # center is not an "additional" center: want to ensure its
        # fragment is included in frag_nums
        # This shouldn't be necessary, but will print if issue
        log.warn("Problem with center %s", i)
        add_center_ind = heavy_atoms_ind[heavy_atoms.index(i)]
        in_frag_list = []
        for frag in frag_auto_f:
            if add_center_ind in frag:
                in_frag_list.append(frag_auto_f.index(frag))
        for f in in_frag_list:
            if f not in frag_nums:
                frag_nums.append(f)

    # frag_nums contains the fragments that needs to be constructed / evaluated.
    # reduce frag_mixed so that it only contains the fragments we need
    if skip_unnecessary_frag_build:
        frag_mixed.Nfrag = len(frag_nums)
        frag_mixed.fsites = [frag_mixed.fsites[i] for i in frag_nums]
        frag_mixed.centerf_idx = [frag_mixed.centerf_idx[i] for i in frag_nums]
        frag_mixed.ebe_weight = [frag_mixed.ebe_weight[i] for i in frag_nums]

    log.info("BE(1) list of center_atoms and corresponding fragments:\n %s,\n %s", center_atoms, frag_nums)

    mybe = BE(
        mf_mixed,
        frag_mixed,
        lo_method="lowdin",
        eri_file="eri_file_" + str(idx) + ".h5",
    )

    energy_list = []
    ehf_list = []
    indexing = []
    e_count = 0.0
    mybe.pot = [chempot] # Set chemical potential
    if skip_unnecessary_frag_build:
        for idxx in range(len(mybe.Fobjs)):
            tot_e, e_comp = be_func(
                mybe.pot,
                [mybe.Fobjs[idxx]],
                mybe.Nocc,
                solver="CCSD",
                enuc=mybe.enuc,
                frag_energy=True,
                eeval=True,
                hf_veff=mybe.hf_veff,
                only_chem=True,
            )
            energy_list.append(tot_e)
            ehf_list.append(mybe.Fobjs[idxx].ebe_hf)
            indexing.append(centers[idxx])
            for i in mybe.Fobjs[idxx].efac[1]:
                e_count += mybe.Fobjs[idxx]._rdm1[i, i]
    else:
        for idxx, x in enumerate(frag_nums):
            tot_e, e_comp = be_func(
                mybe.pot,
                [mybe.Fobjs[x]],
                mybe.Nocc,
                solver="CCSD",
                enuc=mybe.enuc,
                frag_energy=True,
                eeval=True,
                hf_veff=mybe.hf_veff,
                only_chem=True,
            )
            energy_list.append(tot_e)
            ehf_list.append(mybe.Fobjs[x].ebe_hf)
            indexing.append(centers[idxx])
            for i in mybe.Fobjs[x].efac[1]:
                e_count += mybe.Fobjs[x]._rdm1[i, i]
    return energy_list, ehf_list, mybe.enuc + mybe.E_core, e_count

def number_geom(geom, new):
    """Make numbered geometry file
    Paramaters
    __________
    geom: str
        String with path to geometry file
    new: str
        String with path to write new, numbered geometry file
    """
    lines = []
    w = open(geom, "rt")
    i = -1
    atom_list = []
    for line in w:
        v = line.split()
        if len(v) == 4:
            if v[0] not in atom_list:
                atom_list.append(v[0])
            newline = (
                v[0]
                + str(i)
                + " "
                + str(v[1])
                + " "
                + str(v[2])
                + " "
                + str(v[3])
                + "\n"
            )
            lines.append(newline)
        i += 1
    w.close()

    n = open(new, "wt")
    n.write(str(len(lines)) + "\n")
    n.write("\n")
    for line in lines:
        n.write(line)
    n.close()
    return lines, atom_list


# store all geometry information -- return geometry and list of all atom types
geom_lines, all_atoms_list = number_geom(geom_file, num_geom_file)

# set up the structure
mol_small_basis = pyscf.gto.Mole()
mol_small_basis.atom = num_geom_file
mol_small_basis.basis = small_basis_inp
mol_small_basis.charge = charge_inp
mol_small_basis.build()

# fragment the structure, via the big basis set region BE level
log.note("Initially fragmenting the structure")

if verbose < 1:
    frag_au= fragpart(be_type=big_basis_region, mol=mol_small_basis, print_frags=False)
else:
    frag_au= fragpart(be_type=big_basis_region, mol=mol_small_basis, print_frags=True)

au_frag = frag_au.Frag_atom
au_cen = frag_au.center_atom
au_hlist = frag_au.hlist_atom
au_add_centers = frag_au.add_center_atom
log.note("Big Basis Regions: All Fragments: %s", au_frag)
log.info("Big Basis Regions: 'True' Centers: %s", au_cen)
log.info("Big Basis Regions: Hydrogen List: %s", au_hlist)
log.info("Big Basis Additional Centers: %s", au_add_centers)

# Generate a list of all heavy atoms
heavy_atom = []
heavy_atom_ind = []
for i in geom_lines:
    if i.split()[0][0] != "H":
        heavy_atom.append(i.split()[0])
        heavy_atom_ind.append(int((re.split("(\d+)", i.split()[0]))[1]) - 1)
log.debug("All Heavy Atoms and Indices: %s %s", heavy_atom, heavy_atom_ind)

# Identify atoms which are not considered a "center" in a the big basis region scheme
additional_center = []
for i in heavy_atom_ind:
    if i not in au_cen:
        additional_center.append(i)
log.debug("Atoms not a 'True' center: %s", additional_center)

"""
Set up and run fragments in series or parallel
"""

def costfn(chempot = chempotopt, debug001=False):
    print("POT", chempot, flush=True)
    energies = []
    ehfs = []
    enucs = []
    elec_count = 0.0
    if nprocs == 1:
        for index, fragment in enumerate(au_frag):
            energy, ehf, enuc_core_k, ect = mixed_basis_solver(
                index,
                fragment,
                small_basis_inp,
                big_basis_inp,
                au_frag,
                au_cen,
                au_hlist,
                au_add_centers,
                big_basis_be_inp,
                heavy_atom,
                heavy_atom_ind,
                additional_center,
                geom_lines,
                num_geom_file,
                charge_inp,
                mm_charges,
                mm_coords,
                chempot,
            )
            energies.append(energy)
            ehfs.append(ehf)
            enucs.append(enuc_core_k)
            elec_count += ect
    else:
        pool_mix = Pool(nprocs)

        results = []

        for index, fragment in enumerate(au_frag):
            result = pool_mix.apply_async(
                mixed_basis_solver,
                [
                    index,
                    fragment,
                    small_basis_inp,
                    big_basis_inp,
                    au_frag,
                    au_cen,
                    au_hlist,
                    au_add_centers,
                    big_basis_be_inp,
                    heavy_atom,
                    heavy_atom_ind,
                    additional_center,
                    geom_lines,
                    num_geom_file,
                    charge_inp,
                    mm_charges,
                    mm_coords,
                    chempot,
                ],
            )
            results.append(result)

        [energies.append(result.get()[0]) for result in results]
        [ehfs.append(result.get()[1]) for result in results]
        [enucs.append(result.get()[2]) for result in results]
        elec_count = sum([result.get()[3] for result in results])
        pool_mix.close()

    log.note("Correlation Energies are... %s", energies)
    log.note("HF (elec.)  Energies are... %s", ehfs)
    print("Final e_corr: ", sum(sum(i) for i in energies))
    assert np.all(
        np.isclose(enucs, enucs[0])
    ), "Nuclear+Core energies from different split basis BE fragments differ." + str(enucs)
    print("HF Energy:    ", sum(sum(i) for i in ehfs) + enucs[0])
    log.note("ecount %s; expected %s", elec_count, mol_small_basis.nelec[0])

    return abs(elec_count - mol_small_basis.nelec[0])


if chempotopt: # optimize chem pot
    import scipy
    pot = 0.
    scipy.optimize.minimize(costfn, pot) # dumb optimizer
else: costfn([0.])
