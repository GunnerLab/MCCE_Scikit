"""
Module: mcce_io.py
Module contains a collection of functions to find/access mcce output files.
Module to be used in microstate analysis.
"""

from pathlib import Path
import base


def check_h3_file(h3_filepath):
    """TODO:
    check filename = head3.lst
    check path exists
    """
    pass


def read_conformers(head_3_path: str):
    """Create a list of conformers found in head3.lst file."""

    check_h3_file(head_3_path)
    conformers = []

    lines = open(head_3_path).readlines()
    lines.pop(0)
    for line in lines:
        conf = base.Conformer()
        conf.load_from_head3lst(line)
        conformers.append(conf)

    return conformers


def get_msout_filename(mcce_output_path: str, pH: float, Eh: float):
    """Return the ms_out filename from path, pH and eH values."""

    msout_dir = Path(mcce_output_path).joinpath("ms_out")
    if not msout_dir.exists():
        raise FileNotFound(f"Folder 'ms_out' not found in {msout_dir}")
    if not msout_dir.is_dir():
        raise TypeError(f"'ms_out' must be a directory in {mcce_output_path})")

    prec_ph = 0 if pH.is_integer() else 1
    prec_eh = 0 if Eh.is_integer() else 1

    ms_file = f"pH{pH:.{prec_ph}f}eH{Eh:.{prec_eh}f}ms.txt"
    fname = msout_dir.joinpath(ms_file)
    if not fname.exists():
        raise FileNotFound(f"File {fname} not found in {msout_dir}")

    return fname


def MS_to_PDB(selected_confs: list, ms_index: int, step2_path: str, output_folder: str):
    """
    Create a new pdb file in `output_folder` from the `selected_confs`
    ms_index: index of selected_ms, part of output pdb filename.
    step2_path: path to step2_out.pdb.
    output_folder: path to folder for pdb created from selected_ms.
    """

    file_name = Path(output_folder).joinpath(f"ms_pdb_{ms_index}.pdb")

    with open(step2_path) as pdb:
        # ATOM      1  CA  NTR A0001_001   2.696   5.785  12.711   2.000       0.001      01O000M000 "
        # 0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        #         10        20        30        40        50        60        70        80        90
        # len = 91
        with open(file_name, "w") as output_pdb:
            for line in pdb:
                iCode = line[26]
                if iCode == " ":
                    iCode = "_"

                confID = line[17:20] + line[80:82] + line[21:26] + iCode + line[27:30]
                if confID[3:5] == "BK":
                    output_pdb.write(line)
                    continue

                if confID in selected_confs:
                    output_pdb.write(line)
    return
