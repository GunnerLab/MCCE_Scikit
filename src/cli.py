"""
Module: MCCE_Scikit/src/cli.py

Command line interface to obtain a collection of pdbs from a sample of microstates.

"""

from pathlib import Path
from typing import Union
import base
import ms_sampling as sampling
import click


@click.command()
@click.argument(
    "mcce_output_path",
    type=click.Path(
        exists=True,
        file_okay=False,
        readable=True,
        path_type=Path,
    ),
    help="The folder with files from a MCCE simulation (MCCE executable not required).",
)
@click.argument(
    "pH",
    type=Union[int, float],
    help="The pH point; part of experiemntal variables defining a microstate.",
)
@click.argument(
    "Eh",
    type=Union[int, float],
    default=0,
    show_default=True,
    help="The Eh point; part of experiemntal variables defining a microstate.",
)
@click.argument(
    "MC",
    type=int,
    default=0,
    show_default=True,
    help="The (zero-based) index of the MONTERUNS to use.",
)
@click.option(
    "--overwrite",
    default=False,
    type=bool,
    show_default=True,
    help="`pdbs_from_ms` uses a split msout file (header and MCi files); overwrite means the file will be split anew.",
)
@click.argument(
    "n_sample_size",
    type=int,
    help="The number of microstates to sample (and pdb to create).",
)
@click.argument(
    "ms_sort_by",
    type=str,
    default="energy",
    help="The name referencing the Conformer variable to use when sorting; 'energy'-> 0: Conf.E, 'count'-> 1: Conf.count.",
)
@click.argument(
    "output_dir",
    type=click.Path(
        exists=True,
        file_okay=False,
        readable=True,
        path_type=Path,
    ),
    help="""The parent folder of the folder recieving the created pdb files.
                The actual output folder will be: Path(output_dir).joinpath("pdbs_from_ms")""",
)
@click.option(
    "--list_files",
    default=False,
    type=bool,
    show_default=True,
    help="Whether to list the files created.",
)
def pdbs_from_ms(
    mcce_output_path,
    pH,
    Eh,
    MC,
    overwrite,
    n_sample_size,
    ms_sort_by,
    output_dir,
    list_files,
):
    """ """
    ms = base.MS(
        mcce_output_path, pH, Eh, selected_MC=MC, overwrite_split_files=overwrite
    )
    click.echo(ms)

    sampling.pdbs_from_ms_samples(
        ms,
        mcce_output_path,
        n_sample_size,
        ms_sort_by,
        output_dir,
        list_files=list_files,
    )


if __name__ == "__main__":
    pdbs_from_ms()
