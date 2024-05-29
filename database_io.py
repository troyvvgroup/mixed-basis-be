"""
Functions written for I/O with exported database
"""
import h5py, sys, importlib, numpy
import scine_database as db
import scine_utilities as utils

# Temporary dirty trick to import pipeline functions
path_to_pipeline = "/home/minsik/pipeline/scripts/get-capped-qm-region"
sys.path.append(path_to_pipeline)
pipeline = importlib.import_module("get-capped-qm-region")


# Adapted from Moritz Bensberg's code (pipeline.scripts.get-capped-qm-region.get-capped-qm-region.py)
def access_q4bio_db(
    db_name="q4bio-model2-ligand-solvent", structure_id="6623e69125142a3e6e2156ab"
):
    manager = db.Manager()
    credentials = db.Credentials("tc-dl380-login.ethz.ch", 27017, db_name)
    manager.set_credentials(credentials)
    manager.connect()

    # Retrieve structure
    structure_collection = manager.get_collection("structures")
    structure = structure_collection.get_structure(db.ID(structure_id))

    # Write full structure to XYZ
    utils.io.write("structure.xyz", structure.get_atoms())
    structure_file = "structure.xyz"

    # Write connectivity file
    property_collection = manager.get_collection("properties")
    pipeline.write_connectivity_file("connectivity.dat", property_collection, structure)
    connectivity_file_name = "connectivity.dat"

    # Retrieve all atomic charges and write to file
    property_id = structure.get_property("atomic_charges")
    atomic_charges = db.VectorProperty(property_id, property_collection)
    with open("atomic_charges.csv", "w") as f:
        for charge in atomic_charges.get_data():
            f.write(str(charge) + "\n")
    charge_file = "atomic_charges.csv"

    # Get indices of QM atoms
    property_id = structure.get_property("qm_atoms")
    qm_atoms = db.VectorProperty(property_id, property_collection)
    qm_atom_indices = [int(i) for i in qm_atoms.get_data()]

    # Write atom types files
    annotated_atom_collection = pipeline.annotated_structure_with_residue_labels(
        structure, property_collection
    )
    pipeline.write_atom_types_file(annotated_atom_collection)
    atom_type_file = "atom_types.dat"

    # Write force field files
    property_ids = structure.get_properties("openmm_xml_files")
    force_fields = pipeline.write_xml_parameter_files(property_ids, property_collection)

    manager.disconnect()

    settings = utils.ValueCollection(
        {
            "electrostatic_embedding": True,
            "gaff_atomic_charges_file": charge_file,
            "mm_connectivity_file": connectivity_file_name,
            "gaff_atom_types_file": atom_type_file,
            "openmm_xml_files": force_fields,
            "qm_atoms": qm_atom_indices,
            "molecular_charge": 0,
            "spin_multiplicity": 1,
            "program": "xtb/swoose",
            "method": "gfn2",
        }
    )
    qm_region, mm_charges = pipeline.get_capped_qm_region(settings, structure_file)
    utils.io.write(structure_id + ".xyz", qm_region)
    print("XYZ file for the capped QM region is written as " + structure_id + ".xyz")

    charges = []
    positions = []
    i = 0
    while i < len(mm_charges):
        charges.append(mm_charges[i])
        positions.append((mm_charges[i + 2], mm_charges[i + 3], mm_charges[i + 4]))
        i += 5

    numpy.save(structure_id + "_mmcharge", charges)
    numpy.save(structure_id + "_mmcoords", positions)


def process_q4bio_db_export(db_path, output="output.h5"):
    import re, os, numpy

    # Check that db_path is a file
    assert os.path.isfile(db_path), "Input should be a recognizable path to a file."
    h = h5py.File(output, "w")
    idx_id_pair = {}
    with open(db_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.lower() == "begin":
                pass
            elif line.lower() == "end":
                h5grp.create_dataset("mmcoords", data=mmcoords)
                h5grp.create_dataset("mmcharge", data=mmcharge)
                h5grp.create_dataset("atcoords", data=atcoords)
                h5grp.create_dataset("elements", data=elements)
                h5grp.create_dataset("energy", data=energy)
                h5grp.create_dataset("charge", data=charge)
                h5grp.create_dataset("mult", data=mult)
            elif not re.match("^comment", line) is None:  # comment block
                comment_type = re.findall("^comment ([a-z_]*)", line)[0]
                if comment_type == "structure_label" or comment_type == "global":
                    pass
                elif comment_type == "global":  # charge and multiplicity
                    charge = int(re.findall("^global ([0-9]) ([0-9])")[0])
                    mult = int(re.findall("^global ([0-9]) ([0-9])")[1])
                elif comment_type == "number":
                    index = re.findall("^comment number ([0-9]*)", line)[0]
                elif comment_type == "structure_id":
                    structure_id = re.findall("^comment structure_id (.*)$", line)[0]
                    h5grp = h.create_group(structure_id)
                    idx_id_pair.update({index: structure_id})
                    atcoords = []  #                               [(x, y, z), ...]
                    elements = []  #                               ['C', 'H', ...]
                    mmcoords = []  # pyscf.qmmm.mm_charge format   [(x, y, z), ...]
                    mmcharge = []  # pyscf.qmmm.mm_charge format   [chg1, chg2, ...]
                else:
                    raise ValueError("Unrecognized comment block: " + line)
            elif not re.match("^atom", line) is None:  # data block
                data = line.split()
                if data[5] == "2.0":  # mm charge
                    mmcoords.append((float(data[1]), float(data[2]), float(data[3])))
                    mmcharge.append(float(data[6]))
                elif data[5] == "1.0":  # atom
                    atcoords.append((float(data[1]), float(data[2]), float(data[3])))
                    elements.append(data[4])
                else:
                    raise ValueError("Unrecognized data block: " + line)
            elif not re.match("^energy", line) is None:  # energy line
                energy = float(re.findall("^energy ([0-9.-]*)", line)[0])
            elif not re.match("^charge", line) is None:  # charge line
                # charge = float(re.findall("^charge ([0-9.-]*)", line)[0])
                pass
            else:
                raise ValueError("Unrecognized line: " + line)
    h.attrs.update(idx_id_pair)
    h.close()


if __name__ == "__main__":
    if sys.argv[1] == "process":
        process_q4bio_db_export(sys.argv[2], output=sys.argv[3])
        print("Converted ", sys.argv[2], " and wrote ", sys.argv[3])
    elif sys.argv[1] == "demo":
        import pyscf
        from pyscf import qmmm

        print("Using ", sys.argv[2])
        h = h5py.File(sys.argv[2], "r")
        print("Index ", sys.argv[3])
        print("structure_id ", h.attrs[sys.argv[3]], flush=True)
        mol = pyscf.gto.M()
        mol.charge = h[h.attrs[sys.argv[3]]]["charge"][()]
        mol.atom = [
            [x.decode("utf-8"), y]
            for (x, y) in zip(
                h[h.attrs[sys.argv[3]]]["elements"][()],
                h[h.attrs[sys.argv[3]]]["atcoords"][()],
            )
        ]
        mol.basis = "sto-3g"
        mol.unit = "bohr"
        mol.build()
        mf = qmmm.mm_charge(
            pyscf.scf.RHF(mol),
            h[h.attrs[sys.argv[3]]]["mmcoords"][()],
            h[h.attrs[sys.argv[3]]]["mmcharge"][()],
        )
        mf.kernel()
        print("qm/mm energy ", h[h.attrs[sys.argv[3]]]["energy"][()])
    elif sys.argv[1] == "query":
        # Sample workflow:
        # > python /home/minsik/split-basis-be/database_io.py query q4bio-model2-ligand-solvent 1 > structure_id
        # > cat structure_id
        # > python /home/minsik/split-basis-be/database_io.py q4bio-model2-ligand-solvent 6623e67f25142a3e6e215679
        # See /home/q4bio/splitting_basis/MMtest_q4bio-model2-ligand-solvent for sample SLURM script
        manager = db.Manager()
        credentials = db.Credentials("tc-dl380-login.ethz.ch", 27017, sys.argv[2])
        manager.set_credentials(credentials)
        manager.connect()
        structure_collection = manager.get_collection("structures")
        from json import dumps
        import itertools

        selection = {
            "$and": [
                {"label": "user_guess"},
                {"properties.qm_atoms": {"$exists": True}},
            ]
        }
        n_ids = int(sys.argv[3])
        print("Printing first", n_ids, "structure IDs from", str(sys.argv[2]))
        print(
            [
                str(x.id())
                for x in itertools.islice(
                    structure_collection.iterate_structures(dumps(selection)), n_ids
                )
            ]
        )
        manager.disconnect()
    else:
        # Creates three outputs:
        # structure_id.xyz              XYZ of QM region
        # structure_id_mmcharge.npy     List of MM charges
        # structure_id_mmcoords.npy     List of MM coordinates
        # In pyscf:
        # > mol = pyscf.gto.M(); mol.atm = "structure_id.xyz"
        # > ... (Do split basis stuff)
        # > mmcoords = numpy.load("structure_id_mmcoords.npy"); mmcharge = numpy.load("structure_id_mmcharge.npy")
        # > mf = qmmm.mm_charge(pyscf.scf.RHF(mol), mmcoords, mmcharge)
        db_name = sys.argv[1]
        structure_id = sys.argv[2]
        access_q4bio_db(db_name=db_name, structure_id=structure_id)
        print("Processed " + str(db_name) + "/" + str(structure_id))
