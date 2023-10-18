"""
Module: base.py
Contains all classes for creating MCCE objects:
 - Conformer
 - Microstate
 - MS
"""

from collections import namedtuple
from pathlib import Path
import numpy as np
import mcce_io as io
import constants as cst


class Conformer:
    def __init__(self):
        self.iconf = 0
        self.confid = ""
        self.resid = ""
        self.crg = 0.0
        # self.occ = 0.0  # not used

    def load_from_head3lst(self, line):
        fields = line.split()
        self.iconf = int(fields[0]) - 1
        self.confid = fields[1]
        self.resid = self.confid[:3] + self.confid[5:11]
        # self.occ = 0.0
        self.crg = float(fields[4])


class Microstate:
    def __init__(self, state, E, count):
        self.state = state
        self.E = E
        self.count = count


Micro_tup = namedtuple("Microstate", ["E", "count", "state"])


class MS:
    """
    FIXED confs = "always chosen", so each MS is comprised of
    all FIXED + 1 conf from each free residue.
    """

    def __init__(self, mcce_output_path, pH, Eh):
        self.T = cst.ROOMT  # 273.15
        self.pH = pH
        self.Eh = Eh
        self.N_ms = 0
        self.N_uniq = 0
        self.lowest_E = 0.0
        self.highest_E = 0.0
        self.average_E = 0.0

        self.fixed_iconfs = []
        self.fixed_crg = 0.0
        self.fixed_ne = 0.0
        self.fixed_nh = 0.0
        self.free_residues = []  # free residues, referred by conformer indices
        self.iconf2ires = {}  # from conformer index to free residue index
        self.microstates = {}
        self.conformers = []
        self.load_msout(mcce_output_path)

    @staticmethod
    def check_next_lines(lines: list):
        """Continue reading lines until non-comment line found.
        Return a data line and the remaining lines.
        """
        while True:
            line = lines.pop(0).strip()
            if line and line[0] != "#":
                break
        return line, lines

    def load_msout(self, mcce_output_path):
        """Populate class variables from mcce_output_path/ms_out/ms_file"""
        fname = io.get_msout_filename(mcce_output_path, self.pH, self.Eh)

        lines = open(fname).readlines()
        line, lines = self.check_next_lines(lines)

        fields = line.split(",")
        for field in fields:
            key, value = field.split(":")
            key = key.strip().upper()
            value = float(value)
            if key == "T":
                self.T = value
            elif key == "PH":
                self.pH = value
            elif key == "EH":
                self.Eh = value

        # second line, confirm this is from Monte Carlo sampling
        line, lines = self.check_next_lines(lines)

        key, value = line.split(":")
        if key.strip() != "METHOD" or value.strip() not in cst.VALID_MC_METHODS:
            raise ValueError(
                f"This file: {fname} is not a valid Monte Carlo microstate file."
            )

        # Third line, fixed conformer indicies
        line, lines = self.check_next_lines(lines)

        self.fixed_iconfs = [int(i) for i in line.split(":")[1].split()]

        # 4th line, free residues
        line, lines = self.check_next_lines(lines)

        residues = line.split(":")[1].split(" ;")
        self.free_residues = [[int(i) for i in group.split()] for group in residues]
        self.iconf2ires = {
            iconf: idx for idx, res in enumerate(self.free_residues) for iconf in res
        }

        # find the next MC record
        found_mc = False
        newmc = False
        self.N_ms = 0

        for line in lines:
            if line.find("MC:") == 0:  # ms starts
                found_mc = True
                newmc = True
                continue

            if newmc:
                current_state = [int(c) for c in line.split(":")[1].split()]
                newmc = False
                continue

            if found_mc:
                fields = line.split(",")
                if len(fields) < 3:
                    continue

                state_e = float(fields[0])
                count = int(fields[1])
                flipped = [int(c) for c in fields[2].split()]
                for ic in flipped:
                    current_state[self.iconf2ires[ic]] = ic

                ms = Microstate(current_state, state_e, count)

                key = ",".join([str(i) for i in ms.state])
                if self.microstates.get(key) is None:
                    self.microstates[key] = ms
                else:
                    self.microstates[key].count += ms.count

        # find N_ms, lowest, highest, averge E
        E_sum = 0.0
        msvalues = self.microstates.values()
        self.N_uniq = len(msvalues)

        self.lowest_E = next(iter(msvalues)).E
        self.highest_E = next(iter(msvalues)).E

        for ms in msvalues:
            self.N_ms += ms.count
            E_sum += ms.E * ms.count
            if self.lowest_E > ms.E:
                self.lowest_E = ms.E
            if self.highest_E < ms.E:
                self.highest_E = ms.E
        self.average_E = E_sum / self.N_ms
