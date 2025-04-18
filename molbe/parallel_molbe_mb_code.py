import numpy as np
import re, pyscf, os, sys
from pbe.pbe import pbe
from pbe.fragment import fragpart
from pbe.autofrag import autogen
from pbe.helper import get_scfObj, get_eri, get_veff
from pbe.solver import be_func
from pyscf import qmmm
import pbe_var
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
    mm_coords = np.load(
        sys.argv[9]
    )  # npy file with a list of MM coordinates (optional)
else:
    mm_charges = None
    mm_coords = None

# Mol-BE pbe_var settings
pbe_var.SCRATCH = "/scratch/q4bio/"
pbe_var.CREATE_SCRATCH_DIR = True

# add numbers to all geometries -- Match mol-be
num_geom_file = geom_file[:-4] + "-numbered.xyz"

"""
Establish all functions
"""


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
    print("Fragment number:", idx)

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

    print("Basis library for frag ", idx, ": ", basis_library)

    # Identify centers of fragment to perform be1 on
    center = str(geom[auto_cen[idx]].split()[0])

    # Add all heavy atoms which are only in this big basis fragment
    centers = [center]
    for i in auto_add_centers[idx]:
        centers.append(str(geom[i].split()[0]))

    print("All centers for fragment ", idx, ": ", centers)
    print("Heavy atoms in fragment ", heavy_atoms)
    print("Heavy atom index", heavy_atoms_ind)

    print("Atoms in fragment", idx, ": ", frag_atoms)

    # Build HF for the full mixed-basis sysetm
    mol_mixed = pyscf.gto.Mole()
    mol_mixed.atom = numbered_geom_file
    mol_mixed.basis = basis_library
    mol_mixed.charge = charge
    mol_mixed.build()

    # Run HF for fragment
    print("Running HF...")
    if not (mm_charges is None) and not (mm_coords is None):
        print("with MM charges...")
        mf_mixed = qmmm.mm_charge(pyscf.scf.RHF(mol_mixed), mm_coords, mm_charges)
    else:
        mf_mixed = pyscf.scf.RHF(mol_mixed)
    mf_mixed.kernel()

    # Run BE1 fragmentation for mixed basis system
    frag_mixed = fragpart(
        mol_mixed.natm, be_type=big_basis_be, mol=mol_mixed, molecule=True
    )
    frag_auto_f, frag_auto_c, frag_auto_h, frag_add_centers = autogen(
        mol=mol_mixed,
        kpt=None,
        be_type=big_basis_be,
        molecule=True,
        split_basis_return=True,
        print_frags=False
    )

    center_atoms = [geom[i].split()[0] for i in frag_auto_c]
    frag_nums = []; center_atoms_retry = [] # for linear chains, this is fine (making sure center is included)
    # center_atoms_retry is used to ensure that the center is included in the fragment list if it is not an "additional" center
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
        print("Problem with center", i)
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

    print(
        "BE1 list of center_atoms and corresponding fragments", center_atoms, frag_nums
    )

    mybe = pbe(
        mf_mixed,
        frag_mixed,
        lo_method="lowdin",
        molecule=True,
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
            # print("tot_e, e_comp", tot_e, e_comp)
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
            # print("tot_e, e_comp", tot_e, e_comp)
            energy_list.append(tot_e)
            ehf_list.append(mybe.Fobjs[x].ebe_hf)
            indexing.append(centers[idxx])
            for i in mybe.Fobjs[x].efac[1]:
                e_count += mybe.Fobjs[x]._rdm1[i, i]
    return energy_list, ehf_list, mybe.enuc + mybe.E_core - mybe.ek, e_count

def number_geom(geom, new):
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
print("geom_lines", geom_lines)
print("all_atoms_list", all_atoms_list)
# set up the structure
mol_small_basis = pyscf.gto.Mole()
mol_small_basis.atom = num_geom_file
mol_small_basis.basis = small_basis_inp
mol_small_basis.charge = charge_inp
mol_small_basis.build()

# fragment the structure
print("Initially fragmenting the structure")
au_frag, au_cen, au_hlist, au_add_centers = autogen(
    mol=mol_small_basis,
    kpt=None,
    be_type=big_basis_region,
    molecule=True,
    split_basis_return=True,
    print_frags=True
)

print("Big Basis Regions: All Fragments", au_frag)
print("Big Basis Regions: 'True' Centers", au_cen)
print("Big Basis Regions: Hydrogen List", au_hlist)
print("Big Basis EBE Weight:", au_add_centers)

# Generate a list of all heavy atoms
heavy_atom = []
heavy_atom_ind = []
for i in geom_lines:
    if i.split()[0][0] != "H":
        heavy_atom.append(i.split()[0])
        heavy_atom_ind.append(int((re.split("(\d+)", i.split()[0]))[1]) - 1)
print("All Heavy Atoms and Indices", heavy_atom, heavy_atom_ind)

# Identify atoms which are not considered a "center" in a the big basis region scheme
additional_center = []
for i in heavy_atom_ind:
    if i not in au_cen:
        additional_center.append(i)
print("Atoms not a 'True' center", additional_center)

#sys.exit()
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

    print("Correlation Energies are...", energies, flush=True)
    print("HF (elec.)  Energies are...", ehfs, flush=True)
    print("e_corr ", sum(sum(i) for i in energies), flush=True)
    assert np.all(
        np.isclose(enucs, enucs[0])
    ), "Nuclear+Core energies from different split basis BE fragments differ." + str(enucs)
    print("ehf    ", sum(sum(i) for i in ehfs) + enucs[0], flush=True)
    print("ecount ", elec_count, "; expected ", mol_small_basis.nelec, flush=True)

    return abs(elec_count - mol_small_basis.nelec[0])


if chempotopt: # optimize chem pot
    import scipy
    pot = 0.
    scipy.optimize.minimize(costfn, pot) # dumb optimizer
else: costfn([0.])
