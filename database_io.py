"""
Functions written for I/O with exported database
"""
import h5py, sys


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
                elif comment_type == "global": # charge and multiplicity
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
                #charge = float(re.findall("^charge ([0-9.-]*)", line)[0])
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
