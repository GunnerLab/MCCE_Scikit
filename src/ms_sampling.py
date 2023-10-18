"""
Module: ms_sampling

"""

from collections import namedtuple
from pathlib import Path
import numpy as np
import base
import mcce_io as io


Micro_tup = namedtuple("Microstate", ["E", "count", "state"])


def sort_microstate_list(ms_list: list, by: str = None):
    """Sort a list of ms objects ([base.MC.microstate,..]) by='energy' or by='count'."""
    if by is None:
        raise ValueError("Argument `by` is required.")

    by = by.lower()
    if by not in ["energy", "count"]:
        raise ValueError(f"Values for `by` are 'energy' or 'count'; Given: {by}")
    idx = 0
    if by == "count":
        idx = 1

    return sorted(ms_list, key=lambda x: x[idx])


def get_regular_samples(size: int, sorted_ms_list) -> tuple:
    """Implement a 'regular sampling' of all microstates.
    size: sample size
    sorted_ms_list: sorted list of base.MC.microstate
    Return: a tuple
        cum sum of ms.count in sorted_ms_list, array of indices for selection
    """
    n_counts = 0.0
    ms_count_values = []

    for m in sorted_ms_list:
        n_counts += m.count
        ms_count_values.append(m.count)

    full_count_list = list(np.arange(n_counts))
    ms_cum_sum = np.cumsum(ms_count_values)
    X = n_counts - size
    Y = n_counts / size

    return ms_cum_sum, np.arange(size, X, Y)


def get_selected_confs(ms, conformers, selected_ms):
    """Return the list of conformers for selected_ms.
    ms: a base.MC object instance
    conformers: list of Conformer objects as output by `io.read_conformers`.
    selected_ms: a single ms from base.MC.microstate list.
    """
    return [
        conf.confid
        for conf in conformers
        if conf.iconf in selected_ms.state or conf.iconf in ms.fixed_iconfs
    ]
