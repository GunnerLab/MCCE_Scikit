# MCCE_Scikit/src/base.py

"""Contains all classes for creating MCCE objects:
 - Conformer
 - Microstate
 - Micro_tup: Microste as a namedtuple, used for reading only
 - MS

"""

from collections import namedtuple
from pathlib import Path
from typing import Union
from datetime import datetime
import numpy as np
import zlib
import mcce_io as io
import constants as cst


class Conformer:
    def __init__(self):
        self.iconf = 0
        self.confid = ""
        self.resid = ""
        self.crg = 0.0

    def load_from_head3lst(self, line):
        fields = line.split()
        self.iconf = int(fields[0]) - 1
        self.confid = fields[1]
        self.resid = self.confid[:3] + self.confid[5:11]
        self.crg = float(fields[4])


class Microstate:
    def __init__(self, state, E, count):
        self.stateid = zlib.compress(" ".join([str(x) for x in state]).encode())
        self.E = E
        self.count = count

    def state(self):
        return [int(i) for i in zlib.decompress(self.stateid).decode().split()]


Micro_tup = namedtuple("Microstate", ["E", "count", "state"])


class MS:
    """Uses split ms_out files."""

    def __init__(
        self,
        mcce_output_path: str,
        pH: Union[int, float],
        Eh: Union[int, float],
        selected_MC: int = 0,
        overwrite_split_files: bool = False,
    ):
        """MS.init

        Parameters:
            mcce_output_path (str): A MCCE simulation output folder.
            pH (int or float): A pH point.
            Eh (int or float): A Eh point.
            selected_MC (int): The index of an MC run; one of `range(constants.MONTERUNS)`.
            overwrite_split_files (bool): whether to redo the splitting of msout_file.
        """

        self.mcce_out = Path(mcce_output_path)
        self.selected_MC = selected_MC
        self.overwrite_split_files = overwrite_split_files
        self.T = cst.ROOMT
        self.pH = pH
        self.Eh = Eh
        self.method = ""
        self.conformers = []
        self.iconf_by_confname = {}
        self.fixed_iconfs = []
        self.fixed_residue_names = []
        self.fixed_crg = 0.0
        self.fixed_ne = 0.0
        self.fixed_nh = 0.0
        self.free_residues = []  # free residues, referred by conformer indices
        self.free_residue_names = []
        self.ires_by_iconf = {}  # index of free residue by index of conf

        self.microstates = {}  # a dict of microstates
        # self.microstates_by_id = {}  # dict
        self.counts = 0
        self.N_ms = 0
        self.N_uniq = 0

        self.fname = io.get_msout_filename(self.mcce_out, self.pH, self.Eh)
        self.msout_file_dir, created = io.mkdir_from_msout_file(self.fname)
        if created:
            io.split_msout_file(self.mcce_out, self.pH, self.Eh)
        elif self.overwrite_split_files:
            io.clear_folder(self.msout_file_dir)
            io.split_msout_file(self.mcce_out, self.pH, self.Eh)

        self._get_data()

    def __repr__(self):
        return f"""{type(self).__name__}("{self.mcce_out}", {self.pH}, {self.Eh}, selected_MC={self.selected_MC}, overwrite_split_files={self.overwrite_split_files})"""

    def _get_conformer_data(self):
        """Populate class vars: conformers, iconf_by_confname."""
        head3_path = self.mcce_out.joinpath("head3.lst")
        io.check_path(head3_path)
        self.conformers, self.iconf_by_confname = io.read_conformers(head3_path)

        return

    def _get_header_data(self):
        """Populate class vars: T, pH, Eh, method, fixed_confs, free_residues,
        free_residue_names, and ires_by_iconf.
        """

        steps_done = {"exper": False, "method": False, "fixed": False, "free": False}

        header_file = self.msout_file_dir.joinpath("header")
        with open(header_file) as fh:
            for nl, line in enumerate(fh):
                if not steps_done["exper"]:
                    fields = line.split(",")
                    for field in fields:
                        parts = field.split(":")
                        key = parts[0].upper().strip()
                        try:
                            value = float(parts[1])
                        except ValueError:
                            print(
                                f"Unrecognized experimental value \
                                  (number expected), found: '{parts[1]}')."
                            )
                        if key == "T":
                            self.T = value
                        elif key == "PH":
                            self.pH = value
                        elif key == "EH":
                            self.Eh = value
                        else:
                            raise ValueError(
                                f"Unrecognized experimental condition part: {key}"
                            )
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
                    self.method = value.strip()
                    steps_done["method"] = True
                    continue

                if not steps_done["fixed"]:
                    _, fields = line.split(":")
                    self.fixed_confs = [int(x) for x in fields.strip("\n").split()]
                    # TODO: check this:
                    self.fixed_residue_names = [
                        self.conformers[fc].resid for fc in self.fixed_confs
                    ]
                    steps_done["fixed"] = True
                    continue

                if not steps_done["free"]:
                    n_res, fres = line.split(":")
                    self.free_residues = [
                        [int(n) for n in grp.strip().split()]
                        for grp in fres.strip(" ;\n").split(";")
                    ]
                    if len(self.free_residues) != int(n_res):
                        msg = "Mismatch between the number of free residues indicator"
                        msg = (
                            msg
                            + " and the number of residues listed on the same line.\n"
                        )
                        raise ValueError(msg + f"\t{line}")

                    self.free_residue_names = [
                        self.conformers[g[0]].resid for g in self.free_residues
                    ]
                    for ires, res in enumerate(self.free_residues):
                        for iconf in res:
                            self.ires_by_iconf[iconf] = ires
                    steps_done["fee"] = True

        return

    def _get_mc_data(self):
        """Populate class vars microstates with the data in a MC file identified
        in `self.selected_MC`.
        """

        microstates_by_id = {}

        MC_file = self.msout_file_dir.joinpath(f"MC{self.selected_MC}")
        # lines = open(MC_file).readlines()

        with open(MC_file) as fh:
            for nl, line in enumerate(fh):
                line = line.strip()

                if nl == 0:
                    _, confs = line.split(":")
                    current_state = [int(c) for c in confs.split()]
                    if not current_state:
                        msg = "The current ms state line cannot be empty.\n"
                        msg = msg + f"\tProblem line in {MC_file}: {nl}"
                        raise ValueError(msg)
                    continue

                fields = line.split(",")
                if len(fields) >= 3:
                    state_e = float(fields[0])
                    count = int(fields[1])
                    # flipped confs:
                    for ic in [int(c) for c in fields[2].split()]:
                        current_state[self.ires_by_iconf[ic]] = ic

                    ms = Microstate(current_state, state_e, count)

                    if ms.stateid in microstates_by_id.keys():
                        microstates_by_id[ms.stateid].count += ms.count
                    else:
                        microstates_by_id[ms.stateid] = ms
                    self.counts += ms.count

        self.microstates = list(microstates_by_id.values())

        return

    def _get_data(self):
        """Populate class variables from head3.lst, header and MC records files."""
        self._get_conformer_data()
        self._get_header_data()
        self._get_mc_data()

        return

    def get_occ(self, microstates: list) -> list:
        """Return the average occupancy of conformers in each microstates."""

        conf_occ = np.zeros(len(self.conformers))
        total_counts = 0
        for ms in microstates:
            total_counts += ms.count
            for iconf in ms.state():
                conf_occ[iconf] += ms.count

        return (conf_occ / total_counts).tolist()

    def confnames_by_iconfs(self, iconfs):
        """Return the conformers id given their indices."""

        return [self.conformers[ic].confid for ic in list(iconfs)]

    def select_by_conformer(
        self, microstates: list, conformer_selection: list = None
    ) -> tuple:
        """Select microstates with confomers in conformer_selection list.
        Args:
            microstates (list): List of microstates.
            conformer_selection (list): List of conformers.
        Returns:
            tuple: selected, unselected; if `conformer_selection` is empty,
            returns [], microstates.
        """

        selected = []
        unselected = []
        if conformer_selection is None:
            return selected, microstates

        iconf_in = {self.iconf_by_confname[confid] for confid in conformer_selection}
        for ms in microstates:
            if set(ms.state()) & iconf_in:
                selected.append(ms)
            else:
                unselected.append(ms)

        return selected, unselected

    def select_by_energy(self, microstates: list, energy_range: list = None) -> tuple:
        """Select microstates if their energies is in the range given in `energy_range`.

        Args:
            microstates (list): List of microstates.
            energy_range (list):  [low, high].
        Returns:
            tuple: selected, unselected; if `energy_range` is empty,
            returns [], microstates.
        """

        if len(energy_range) != 2:
            print(
                f"Provide two values for `energy_range`= [low, high]; invalid: {energy_range}"
            )
            return None, None

        selected = []
        unselected = []
        if energy_range is None:
            return selected, microstates
        energy_range.sort()

        for ms in microstates:
            if energy_range[0] <= ms.E < energy_range[1]:
                selected.append(ms)
            else:
                unselected.append(ms)

        return selected, unselected


# here bc circular import if in mcce_io
def get_pdb_remark(ms: MS, ms_index: int):
    """Return a REMARK 250 string to prepend in pdb.

    > REMARK 250 is mandatory if other than X-ray, NMR, neutron, or electron study.
    [Ref]: https://www.wwpdb.org/documentation/file-format-content/format33/remarks1.html
    """

    R250 = """
REMARK 250
REMARK 250 EXPERIMENTAL DETAILS
REMARK 250  EXPERIMENT TYPE                : MCCE simulation
REMARK 250  DATE OF DATA COLLECTION        : {DATE}
REMARK 250 REMARK: DATE OF DATA COLLECTION is the date this pdb was created.
REMARK 250 EXPERIMENTAL CONDITIONS
REMARK 250  TEMPERATURE                    : {T}
REMARK 250  PH                             : {PH}
REMARK 250  EH                             : {EH}
REMARK 250  METHOD                         : {METHOD}
REMARK 250  SELECTED MONTERUN              : {MC}
REMARK 250  SELECTED MICROSTATE INDEX      : {MS}
REMARK 250
REMARK 250 REMARK: DATE OF DATA COLLECTION : Date this pdb was created.
"""
    dte = datetime.today()
    remark = R250.format(
        DATE=dte.strftime("%d-%b-%y"),
        T=f"{ms.T:.2f}",
        PH=f"{ms.pH:.2f}",
        EH=f"{ms.Eh:.2f}",
        METHOD=f"{ms.method.upper()}",
        MC=f"{ms.selected_MC}",
        MS=ms_index,
    )
    return remark
