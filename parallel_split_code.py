import numpy as np
import re, pyscf, os, sys
from pbe.pbe import pbe
from pbe.fragment import fragpart
from pbe.autofrag import autogen
from pbe.helper import get_scfObj, get_eri, get_veff
from pbe.solver import be_func
import pbe_var
from multiprocessing import Pool

"""
Prepare Calculation Setup
"""

geom_file = sys.argv[1] # Geometry File
big_basis_region = sys.argv[2] # How large of an area has a big basis?
small_basis_inp = 'sto-3g' # Small basis set of the full system
big_basis_be_inp = 'be1' # BE level in big basis region
big_basis_inp = sys.argv[3] # Basis set of big basis region
nprocs = int(sys.argv[4]) # Number of available processors
charge_inp = int(sys.argv[5]) # Charge of system

# Mol-BE pbe_var settings
pbe_var.SCRATCH = "/scratch/q4bio/"
pbe_var.CREATE_SCRATCH_DIR = True

# add numbers to all geometries -- Match mol-be
num_geom_file = geom_file[:-4]+"-numbered.xyz"

"""
Establish all functions
"""
def mixed_basis_solver(idx, f, small_basis, big_basis, auto_frag, auto_cen, auto_hlist, big_basis_be,
			heavy_atoms, heavy_atoms_ind, add_center, geom, numbered_geom_file, charge):

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

	print("Basis library for frag ",idx,": ",basis_library)

	# Identify centers of fragment to perform be1 on
	center = str(geom[auto_cen[idx]].split()[0]) 

	# Add all atoms which are only in this fragment
	centers = [center]
	for i in add_center:
		if i in auto_frag[idx]:
			centers.append(heavy_atoms[heavy_atoms_ind.index(i)])
	print("All centers for fragment ",idx,": ", centers)

	print("Atoms in fragment",idx,": ", frag_atoms)

	# Build HF for the full mixed-basis sysetm
	mol_mixed = pyscf.gto.Mole()
	mol_mixed.atom = numbered_geom_file
	mol_mixed.basis = basis_library
	mol_mixed.charge = charge
	mol_mixed.build()

	# Run HF for fragment
	print("Running HF...")
	mf_mixed = pyscf.scf.RHF(mol_mixed)
	mf_mixed.kernel()

	# Run BE1 fragmentation for mixed basis system
	frag_mixed = fragpart(mol_mixed.natm, be_type=big_basis_be, mol=mol_mixed, molecule=True)
	frag_auto_f, frag_auto_c, frag_auto_h = autogen(mol=mol_mixed,kpt=None,be_type=big_basis_be,
							molecule=True,split_basis_return=True,print_frags=False)

	center_atoms = [geom[i].split()[0] for i in frag_auto_c]
	frag_nums = []
	for i in centers:
		frag_nums.append(center_atoms.index(i))

	print("BE1 list of center_atoms and corresponding fragments", center_atoms, frag_nums)

	mybe = pbe(mf_mixed, frag_mixed, lo_method="lowdin", molecule=True, eri_file="eri_file_"+str(idx)+".h5")

	energy_list = []
	indexing = []
	for idxx, x in enumerate(frag_nums):
		tot_e, e_comp = be_func(mybe.pot, [mybe.Fobjs[x]], mybe.Nocc, solver='CCSD', enuc=mybe.enuc, 
					frag_energy=True, eeval=True, hf_veff=mybe.hf_veff)
		#print("tot_e, e_comp", tot_e, e_comp)
		energy_list.append(tot_e)
		indexing.append(centers[idxx])
	return energy_list #, indexing)

def number_geom(geom, new):
	lines=[]
	w = open(geom,'rt')
	i = -1
	atom_list=[]
	for line in w:
		v = line.split()
		if len(v) == 4:
			if v[0] not in atom_list: atom_list.append(v[0])
			newline = v[0]+str(i)+" "+str(v[1])+" "+str(v[2])+" "+str(v[3])+"\n"
			lines.append(newline)
		i+=1
	w.close()

	n = open(new, 'wt')
	n.write(str(len(lines))+'\n')
	n.write('\n')
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

# fragment the structure
print("Initially fragmenting the structure")
au_frag, au_cen, au_hlist = autogen(mol=mol_small_basis,kpt=None,be_type=big_basis_region,
					molecule=True,split_basis_return=True,print_frags=True)

print("Big Basis Regions: All Fragments",au_frag)
print("Big Basis Regions: 'True' Centers", au_cen)
print("Big Basis Regions: Hydrogen List", au_hlist)

# Generate a list of all heavy atoms
heavy_atom=[]
heavy_atom_ind=[]
for i in geom_lines:
	if i.split()[0][0] != 'H':
		heavy_atom.append(i.split()[0])
		heavy_atom_ind.append(int((re.split('(\d+)',i.split()[0]))[1])-1)
print("All Heavy Atoms and Indices", heavy_atom, heavy_atom_ind)

# Identify atoms which are not considered a "center" in a the big basis region scheme
additional_center=[]
for i in heavy_atom_ind:
	if i not in au_cen:
		additional_center.append(i)
print("Atoms not a 'True' center", additional_center)


"""
Set up and run fragments in series or parallel
"""

energies = []
if nprocs == 1:
	for index, fragment in enumerate(au_frag):
		energy = mixed_basis_solver(index, fragment, small_basis_inp, big_basis_inp, au_frag,
						au_cen, au_hlist, big_basis_be_inp, heavy_atom, heavy_atom_ind,
						additional_center, geom_lines, num_geom_file, charge_inp)
		energies.append(energy)
else:
	pool_mix = Pool(nprocs)

	results = []

	for index, fragment in enumerate(au_frag):
		result = pool_mix.apply_async(mixed_basis_solver, [index, fragment, small_basis_inp,
			big_basis_inp, au_frag, au_cen, au_hlist, big_basis_be_inp, heavy_atom,
			heavy_atom_ind, additional_center, geom_lines, num_geom_file, charge_inp])
		results.append(result)

	[energies.append(result.get()) for result in results]
	pool_mix.close()

print("Energies are...", energies, flush=True)
#print("Indexes are...", indexing)
print("e_corr ", sum(sum(i) for i in energies), flush=True)
#np.savetxt("Energies.csv",energies)
#np.savetxt("Indexes.csv",indexing)


