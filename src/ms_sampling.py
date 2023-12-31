"""
Module `ms_sampling`

The module contains the following functions:
 - get_pdb_remark
 - get_selected_confs
 - pdbs_from_ms_samples
 - sample_microstates
 - sort_microstate_list
"""

from pathlib import Path
from datetime import datetime
import numpy as np
import base
import mcce_io as io
import ms_sampling as sampling


def sort_microstate_list(ms_list: list, by: str = None, reverse=False):
    """Sort a list of Microstate objects by 'energy' or 'count'.
    Args:
        ms_list (list): list of Microstate objects ([base.MS.Microstate,..]).
        by (str): Sort key name, one of "energy" or "count", case insensitive.
        reverse (bool, False): Argument for `sorted` function.
    """

    if by is None:
        raise ValueError("Argument `by` is required; one of ['energy', 'count'].")

    by = by.lower()
    if by not in ["energy", "count"]:
        raise ValueError(f"Values for `by` are 'energy' or 'count'; Given: {by}")
    idx = 0
    if by == "count":
        idx = 1
    ms_values = [[m.E, m.count, m.state] for m in ms_list]

    return sorted(ms_values, key=lambda x: x[idx], reverse=reverse)


def sample_microstates(size: int, sorted_ms_list: list) -> tuple:
    """
    Implement a sampling of all microstates.
    Args:
        size (int): sample size
        sorted_ms_list: sorted [base.MS.Microstate,..]
    Returns:
        tuple: cumsum of ms.count in sorted_ms_list, array of indices for selection
    """

    n_counts = 0.0
    ms_count_values = []

    for ms in sorted_ms_list:
        n_counts += ms[1]
        ms_count_values.append(ms[1])

    X = n_counts - size
    Y = n_counts / size
    count_selection = np.arange(size, X, Y)
    ms_cumsum = np.cumsum(ms_count_values)

    return ms_cumsum, count_selection


def get_selected_confs(ms: base.MS, selected_ms):
    """Return the list of conformer ids for selected_ms.
    Args:
        ms (base.MS): class instance
        selected_ms (int?): A single ms from base.MS.microstates list.
    """

    return [
        conf.confid
        for conf in ms.conformers
        if conf.iconf in selected_ms[2]() or conf.iconf in ms.fixed_iconfs
    ]


def pdbs_from_ms_samples(
    ms: base.MS,
    mcce_dir: str,
    n_sample_size: int,
    ms_sort_by: str,
    output_dir: str = None,
    clear_pdbs_folder: bool = True,
    list_files: bool = False,
) -> None:
    """Create `n_sample_size` MCCE_PDB files in `output_dir/pdbs_from_ms`.

    Args:
        ms (base.MS): A microstate class instance.
        mcce_dir (str): MCCE simulation output folder.
        n_sample_size (int): How many microstates/pdbs.
        ms_sort_by (str): Either 'energy' or 'count'.
        output_dir (str): Output folder path;
                          Folder "output_dir/pdbs_from_ms" will be created if necessary.
        list_files (bool): Whether to list output folder contents.
    """
    ms_sort_by = ms_sort_by.lower()
    if ms_sort_by not in ["energy", "count"]:
        raise ValueError(
            f"Values for `ms_sort_by` are 'energy' or 'count'; Given: {ms_sort_by}"
        )

    mcce_dir = Path(mcce_dir)
    step2_path = mcce_dir.joinpath("step2_out.pdb")
    io.check_path(step2_path)

    if not output_dir or output_dir is None:
        output_dir = ms.msout_file_dir

    pdb_out_folder = Path(output_dir).joinpath("pdbs_from_ms")
    if not pdb_out_folder.exists():
        Path.mkdir(pdb_out_folder)
    elif clear_pdbs_folder:
        io.clear_folder(pdb_out_folder)

    mc_run = ms.selected_MC  # part of pdb name
    sorted_ms_list = sort_microstate_list(ms.microstates, by=ms_sort_by)
    ms_cumsum, count_selection = sample_microstates(n_sample_size, sorted_ms_list)

    # Summarize what's being done:
    print(
        f"Creating n={n_sample_size:,} MCCE_PDB files in {output_dir} from (n) microstates sorted by '{ms_sort_by}'.\n",
        "NOTE: the output pdb will be free of any water molecules in step2_out.pdb.",
    )

    for c in count_selection:
        ms_index = np.where((ms_cumsum - c) > 0)[0][0]
        ms_selection = sorted_ms_list[ms_index]

        confs_for_pdb = get_selected_confs(ms, ms_selection)

        # gather initial data for REMARK section of pdb:
        remark_data = get_pdb_remark(ms, ms_index)
        # write the pdb in the folder
        io.ms_to_pdb(
            confs_for_pdb, ms_index, mc_run, remark_data, step2_path, pdb_out_folder
        )
        # pdb names: = Path(pdb_out_folder).joinpath(f"mc{mc_run}_ms{ms_index}.pdb")

    print("PDB files creation over.")
    if list_files:
        print(f"Files in {pdb_out_folder}:\n")
        io.list_folder(pdb_out_folder)

    return


def get_pdb_remark(ms: base.MS, ms_index: int):
    """Return a REMARK 250 string to prepend in pdb.

    > REMARK 250 is mandatory if other than X-ray, NMR, neutron, or electron study.
    [Ref]: https://www.wwpdb.org/documentation/file-format-content/format33/remarks1.html
    """

    R250 = """
REMARK 250
REMARK 250 EXPERIMENTAL DETAILS
REMARK 250   EXPERIMENT TYPE               : MCCE simulation
REMARK 250   DATE OF DATA COLLECTION       : {DATE}
REMARK 250   REMARK: DATE OF DATA COLLECTION is the date this pdb was created.
REMARK 250 EXPERIMENTAL CONDITIONS
REMARK 250   TEMPERATURE                   : {T:.2f} (K)
REMARK 250   PH                            : {PH:.2f}
REMARK 250   EH                            : {EH:.2f}
REMARK 250   METHOD                        : {METHOD}
REMARK 250   SELECTED MONTERUN             : {MC}
REMARK 250   SELECTED MICROSTATE INDEX     : {MS:,}
REMARK 250   SELECTED MICROSTATE ENERGY    : {E:.2f} (kcal/mol)
REMARK 250
"""
    dte = datetime.today()
    sel_ms = ms.microstates[ms_index]
    remark = R250.format(
        DATE=dte.strftime("%d-%b-%y"),
        T=ms.T,
        PH=ms.pH,
        EH=ms.Eh,
        METHOD=ms.method,
        MC=ms.selected_MC,
        MS=ms_index,
        E=sel_ms.E,
    )

    return remark
