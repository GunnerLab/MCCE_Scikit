"""
Module: mcce_io.py
Module contains a collection of functions to find/access mcce output files.
Module to be used in microstate analysis.
"""

from pathlib import Path
from functools import partial
import pickle
import numpy as np
import base
import constants as cst


np_save = partial(np.save, allow_pickle=False, fix_imports=False)


def check_h3_file(h3_filepath):
    """Check file exists"""
    fp = Path(h3_filepath)
    if not fp.exists():
        raise FileNotFoundError(f"File {fp.name} not found in {fp.parent}")
    return


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


def MS_to_PDB(selected_confs: list, ms_index: int, step2_path: str, output_folder: str):
    """
    Create a new pdb file in `output_folder` from the `selected_confs`
    ms_index: index of selected_ms, part of output pdb filename.
    step2_path: path to step2_out.pdb.
    output_folder: path to folder for pdb created from selected_ms.
    """

    file_name = Path(output_folder).joinpath(f"ms_pdb_{ms_index}.pdb")
    sep = "_"
    with open(step2_path) as pdb:
        # ATOM      1  CA  NTR A0001_001   2.696   5.785  12.711   2.000       0.001      01O000M000 "
        # 0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        #         10        20        30        40        50        60        70        80        90
        # len = 91
        with open(file_name, "w") as output_pdb:
            for line in pdb:
                confID = line[17:20] + line[80:82] + line[21:26] + sep + line[27:30]
                if confID[3:5] == "BK":
                    output_pdb.write(line)
                    continue

                if confID in selected_confs:
                    output_pdb.write(line)
    return


def delete_folder(dir_path: str):
    import shutil

    shutil.rmtree(dir_path)
    return


def get_msout_filename(mcce_output_path: str, pH: float, Eh: float):
    """Return the ms_out filename from path, pH and Eh values."""
    if not Path(mcce_output_path).exists():
        raise FileNotFound(f"Folder not found: {mcce_output_path}")

    msout_dir = Path(mcce_output_path).joinpath("ms_out")
    if not msout_dir.exists():
        raise FileNotFound(f"Folder 'ms_out' not found in {mcce_output_path}")
    if not msout_dir.is_dir():
        raise TypeError(f"'ms_out' must be a directory in {mcce_output_path})")

    prec_ph = 0 if pH.is_integer() else 1
    prec_eh = 0 if Eh.is_integer() else 1

    ms_file = f"pH{pH:.{prec_ph}f}eH{Eh:.{prec_eh}f}ms.txt"
    fname = msout_dir.joinpath(ms_file)
    if not fname.exists():
        raise FileNotFound(f"File {fname} not found in {msout_dir}")

    return fname


def mkdir_from_msout_file(msout_file: str):
    """Create a directory with the same name as msout_file w/o the extension.
    Return its path."""
    msout_file_dir = msout_file.parent.joinpath(msout_file.stem)
    if not msout_file_dir.exists():
        Path.mkdir(msout_file_dir)
        print("Created msout_file_dir.")

    return msout_file_dir


def divide_msout_file(
    mcce_output_path: str, pH: float, Eh: float, overwite: bool = False
):
    """Split the msout file into a "header" portion (preceeding MC:0 line) and temporary MCi file
    for the MC records in a folder created with the name of the msout_file as per arguments.
    """

    fname = get_msout_filename(mcce_output_path, pH, Eh)
    msout_file_dir = mkdir_from_msout_file(fname)

    header_file = msout_file_dir.joinpath("header")
    if header_file.exists() and not overwrite:
        print("The header file already exists. Set `overwrite` to True to replace it.")
        return

    steps_done = {"exper": False, "method": False, "fixed": False, "free": False}

    with open(header_file, "w") as header:
        with open(fname) as fh:
            for nl, line in enumerate(fh):
                line = line.strip()
                if not line:
                    continue
                if line[0] == "#":
                    continue

                if not steps_done["exper"]:
                    if line[0] != "T":
                        raise ValueError(
                            f"The first data line (experiemental variables) must start with T.\n{line}"
                        )

                    # ms_exper
                    header.write(line)
                    steps_done["exper"] = True
                    continue

                if not steps_done["method"]:
                    key, value = line.split(":")
                    if (
                        key.strip() != "METHOD"
                        or value.strip() not in cst.VALID_MC_METHODS
                    ):
                        raise ValueError(
                            f"""This file: {fname} is not a valid Monte Carlo microstate file or the method is unknown.
                            Supported methods are: {cst_VALID_MC_METHODS}."""
                        )
                    # ms_method
                    header.write(line)
                    steps_done["method"] = True
                    continue

                if not steps_done["fixed"]:
                    # fixed conformer indicies
                    header.write(line)
                    steps_done["fixed"] = True
                    continue

                if not steps_done["free"]:
                    # free residues
                    header.write(line)
                    steps_done["free"] = True

                if line.startswith("MC"):
                    break

        return


def split_msout_MC(mcce_output_path: str, pH: float, Eh: float, overwite: bool = False):
    """Split the MC records of msout file into into their own files in
    a folder created with the name of the msout_file as per arguments.
    Note: uses constants.MONTE_RUNS (6) :: not read from run.prm :: TODO.
    """

    fname = get_msout_filename(mcce_output_path, pH, Eh)
    msout_file_dir = mkdir_from_msout_file(fname)

    MC_done = dict((i, False) for i in range(cst.MONTE_RUNS))

    # open files for each MC records:: MCi
    for k in MC_done.keys():
        locals()[f"MC{k}"] = open(f"MC{k}", "w")

    with open(fname) as fh:
        for nl, line in enumerate(fh):
            line = line.strip()
            if not line:
                continue

            if not line.startswith("MC:0"):
                continue

            if not MC_done[0]:
                MC0.write(line)
            if line.startswith("MC:1"):
                MC_done[0] = True
                MC0.close()

            if not MC_done[1]:
                MC1.write(line)
            if line.startswith("MC:2"):
                MC_done[1] = True
                MC1.close()

            if not MC_done[2]:
                MC2.write(line)
            if line.startswith("MC:3"):
                MC_done[2] = True
                MC2.close()

            if not MC_done[3]:
                MC3.write(line)
            if line.startswith("MC:4"):
                MC_done[3] = True
                MC3.close()

            if not MC_done[4]:
                MC4.write(line)
            if line.startswith("MC:5"):
                MC_done[4] = True
                MC4.close()

            if not MC_done[5]:
                MC5.write(line)

    MC5.close()
    MC_done[5] = True
    # temp return:
    return MC_done
