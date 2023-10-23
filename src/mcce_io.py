# MCCE_Scikit/src/mcce_io.py

"""Module contains a collection of functions to find/access mcce output files.

The module contains the following functions:

- `check_path(filepath:str)` - Returns None or FileNotFoundError.
- `read_conformers(head_3_path)` - Returns a tuple: conformers (list), iconf_by_confname (dict).
- `MS_to_PDB(selected_confs,ms_index, step2_path, output_folder) - Creates a pdb file (MCCE format) into output_folder.
- ``
- ``

TODO: Create mcce_pdb2pdb conversion function: mcce -> standard pdb w/REMARK
"""

from pathlib import Path
from typing import Union
from functools import partial
import numpy as np
import base
import constants as cst


def mcce_pdb2pdb(fname):
    """TODO"""
    raise NotImplementedError()


def check_path(filepath: str) -> None:
    """
    Returns:
        None if 'all clear', else error if either the parent folder
        or the file itself does not exists.
    """
    fp = Path(filepath)
    if fp.is_dir():
        if not fp.exists():
            raise FileNotFoundError(f"Folder not found: {fp}")
    else:
        if not fp.parent.exists():
            raise FileNotFoundError(f"Parent folder not found: {fp.parent}")
        if not fp.exists():
            raise FileNotFoundError(f"File {fp.name} not found in {fp.parent}")

    return


def read_conformers(head_3_path: str) -> tuple:
    """Return a list of conformers found in head3.lst file and
    a dictionnary, `iconf_by_confname` :: confid -> index."""

    check_path(head_3_path)
    conformers = []
    iconf_by_confname = {}

    with open(head_3_path) as h3:
        for nl, line in enumerate(h3):
            if nl == 0:
                continue
            if len(line) <= 80:
                continue

            conf = base.Conformer()
            conf.load_from_head3lst(line)
            conformers.append(conf)
            iconf_by_confname[conf.confid] = nl

    return conformers, iconf_by_confname


def MS_to_PDB(
    selected_confs: list,
    ms_index: int,
    mc_run: int,
    remark_data: str,
    step2_path: str,
    output_folder: str,
) -> None:
    """Create a new pdb file in `output_folder` from the `selected_confs`
    Args:
        ms_index (int): Index of selected ms, part of output pdb filename.
        mc_run (int): Index of MC record used, part of output pdb filename.
        step2_path: path to step2_out.pdb.
        output_folder: path to folder for pdb created from selected_ms.
        remark_data_exper (dict): Used to create pdb REMARK section: data from MS instance:
                                  experimental (T, PH, EH, METHOD).
    Returns:
        None: The created file names format is f"mc{mc_run}_ms{ms_index}.pdb".

    Pre-requisites:
        The user must acertain that the `selected_confs` come from the same MCCE output files directory
        as the one provided in `ste2_path`.

    TODO: Add REMARK section in pdb
    """

    sep = "_"
    file_name = Path(output_folder).joinpath(f"mc{mc_run}_ms{ms_index}.pdb")

    # step2_out.pdb format:
    # ATOM      1  CA  NTR A0001_001   2.696   5.785  12.711   2.000       0.001      01O000M000 "
    # 0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    #         10        20        30        40        50        60        70        80        90
    # len = 91
    with open(step2_path) as pdb:
        with open(file_name, "w") as output_pdb:
            output_pdb.write(remark_data)
            for line in pdb:
                confID = line[17:20] + line[80:82] + line[21:26] + sep + line[27:30]
                if confID[3:5] == "BK":
                    output_pdb.write(line)
                    continue

                if confID in selected_confs:
                    output_pdb.write(line)
    return


def clear_folder(dir_path: str):
    """Delete all files in folder."""
    if not dir_path:
        return
    for f in dir_path.iterdir():
        if not f.is_dir():
            f.unlink()
    return


def list_folder(dir_path: str):
    """List the files in folder."""
    for f in dir_path.iterdir():
        if not f.is_dir():
            print("\t", f)
    return


def get_msout_filename(
    mcce_output_path: str, pH: Union[float, int], Eh: Union[float, int]
) -> Path:
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


def mkdir_from_msout_file(msout_file: Path) -> tuple:
    """Create a directory with the same name as msout_file w/o the extension.
    Returns:
        tuple: path, created (bool).
    """

    msout_file_dir = msout_file.parent.joinpath(msout_file.stem)
    exists = msout_file_dir.exists()
    if not exists:
        Path.mkdir(msout_file_dir)
        print("Created msout_file_dir.")

    return msout_file_dir, not exists


def check_msout_split(msout_file_dir: Path) -> bool:
    """Return True if the header file exist.
    Assumed: if it does, all the MCi files exist as well."""

    return msout_file_dir.joinpath("header").exists()


def split_msout_file(
    mcce_output_path: str, pH: float, Eh: float, overwrite: bool = False
):
    """Split the msout file into a "header" portion (preceeding MC:0 line) and MCi files
    for the MC records in a folder created with the name of the msout_file as per arguments.
    Note: Each file created is free of comments or blank lines.
    """

    fname = get_msout_filename(mcce_output_path, pH, Eh)
    msout_file_dir, created = mkdir_from_msout_file(fname)

    if check_msout_split(msout_file_dir) and not overwrite:
        print(
            "The ms_out file is already split into header and MCi files. Set `overwrite` to True to replace them."
        )
        return

    steps_done = {"exper": False, "method": False, "fixed": False, "free": False}
    MC_done = dict((i, False) for i in range(cst.MONTE_RUNS))

    header_file = msout_file_dir.joinpath("header")
    header = open(header_file, "w")
    # open files for each MC records:: MCi
    for k in MC_done.keys():
        p = msout_file_dir.joinpath(f"MC{k}")
        globals()[f"MC{k}"] = open(p, "w")

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
                header.write(line + "\n")
                steps_done["exper"] = True
                continue

            if not steps_done["method"]:
                key, value = line.split(":")
                if key.strip() != "METHOD" or value.strip() not in cst.VALID_MC_METHODS:
                    raise ValueError(
                        f"""This file: {fname} is not a valid Monte Carlo microstate file or the method is unknown.
                        Supported methods are: {cst_VALID_MC_METHODS}."""
                    )
                header.write(line + "\n")
                steps_done["method"] = True
                continue

            if not steps_done["fixed"]:
                header.write(line + "\n")
                steps_done["fixed"] = True
                continue

            if not steps_done["free"]:
                header.write(line + "\n")
                steps_done["free"] = True
                continue

            if line.startswith("MC:0"):
                header.close()
                continue

            if not MC_done[0]:
                MC0.write(line + "\n")
                if line.startswith("MC:1"):
                    MC_done[0] = True
                    MC0.close()
                continue

            if not MC_done[1]:
                MC1.write(line + "\n")
                if line.startswith("MC:2"):
                    MC_done[1] = True
                    MC1.close()
                continue

            if not MC_done[2]:
                MC2.write(line + "\n")
                if line.startswith("MC:3"):
                    MC_done[2] = True
                    MC2.close()
                continue

            if not MC_done[3]:
                MC3.write(line + "\n")
                if line.startswith("MC:4"):
                    MC_done[3] = True
                    MC3.close()
                continue

            if not MC_done[4]:
                MC4.write(line + "\n")
                if line.startswith("MC:5"):
                    MC_done[4] = True
                    MC4.close()
                continue

            if not MC_done[5]:
                MC5.write(line + "\n")

    MC5.close()
    MC_done[5] = True

    return
