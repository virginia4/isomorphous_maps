# %matplotlib inline
from __future__ import absolute_import, division, print_function
import argparse
import numpy as np
import pandas as pd

# import stat_ded.structure_factors_from_data, stat_ded.scaling, stat_ded.weights
from stat_ded.kk import structure_factors_from_data
from stat_ded.scaling import *
from  stat_ded.weights import *


def main():

    """
    Author: Virginia Apostolopopoulou, virginia.apostolopoulou@cfel.de

    Code to calculate  difference structure factors and weights for
    isomorphous difference maps (FLIGHT_obs - FDARK_obs).

    The algorithm is based on the Marius Schmidt' FORTRAN script and here it has been
    converted to python.

    Marius Schmidt email: smarius@uwm.edu
    ***********************************************************************
    This script also uses a significant number of functions and libraries from the 
    Computational Crystallographic toolbox (cctbx).

    For documentations and further information: https://cctbx.github.io
    Functions has been altered to fit the needs of this script.

    * fun structure_factors_from_data is a modified version of
    cctbx_project/mmtbx/command_line/maps.py/function:run() (which is the
    equivalent to phenix.maps)
    """

    # options for weight estimation:
    # * no: for no weighted maps
    # * Marius: Marius weights
    # * Dalton: Dalton weights (probably will remove cause is the same as Marius)
    # * Zhong: Zhong weights

    # # Read input arguments from terminal
    parser = argparse.ArgumentParser()
    parser.add_argument("--PHI_file", type=str, help="enter file for dark phases")
    parser.add_argument(
        "--F_file1", type=str, help="enter data file for dark observable data"
    )
    parser.add_argument(
        "--F_file2", type=str, help="enter data file for light observable data"
    )

    parser.add_argument(
        "--w_method",
        default="q_weights",
        type=str,
        nargs="?",
        action="store",
        choices=["no_weights", "q_weights", "Zhong_weights"],
        help="options: no_weights, q_weights (default), Zhong_weights",
    )

    parser.add_argument(
        "--backgr_correction",
        default=1,
        type=float,
        nargs="?",
        action="store",
        help="use a Background Density Correction factor (BDC) when using no weights -- Default BDC=1",
    )

    parser.add_argument(
    "--res_max",
    default=None,
    type=float,
    nargs="?",
    action="store",
    help="set max resolution cutoff - default: take all reso",
    )

    parser.add_argument(
        "--res_min",
        default=None,
        type=float,
        nargs="?",
        action="store",
        help="set min resolution cutoff - default: take all reso",
    )

    parser.add_argument(
        "--scaling",
        default="aniso",
        type=str,
        nargs="?",
        action="store",
        choices=["linear", "iso", "aniso"],
        help="options: linear (linear), isotropic (iso), anisotropic (aniso), (default option)",
    )

    args = parser.parse_args()

    print("Read phases from:", args.PHI_file)
    print("Read dark structure factors from:", args.F_file1)
    print("Read light structure factors from:", args.F_file2)
    print("Calculate isomorphous density maps using", args.w_method, "weights")
    print("Resolution max cut off:", args.res_max)
    print("Resolution min cut off:", args.res_min)
    print("Scaling F_obs to F_calc:", args.scaling)

    res_cutoff_max = args.res_max
    res_cutoff_min = args.res_min
    BDC = args.backgr_correction

    # assign list with the names of the data files for which you want to calculate 
    # the structure factors first arg of the list will contain the file name that 
    # has the phases and the second file that has the observed intensities 
    # (or amplitudes) and sigmas.
    arguments_dark = [args.PHI_file, args.F_file1]

    # call fun structure_factors_from_data to get structure factors from F_file1
    # structure_factors_from_data() will produce an *.mtz file with the amplitudes form F_file1.
    # filename_test_dark carries the name of the product file.
    filename_test_dark = structure_factors_from_data(
        args=arguments_dark[0:], output_file_name="crystal_data_1"
    )
    # read stucture_factors_from_data() output file in a pd.Dataframe
    df_dark_fmodel = rs.read_mtz(filename_test_dark)

    # call fun structure_factors_from_data again to get structure factors from F_file2
    arguments_light = [args.PHI_file, args.F_file2]

    filename_test_light = structure_factors_from_data(
        args=arguments_light[0:], output_file_name="crystal_data_2"
    )
    # read stucture_factors_from_data() output file in a pd.Dataframe
    df_light_fmodel = rs.read_mtz(filename_test_light)

    # choose columns from file and assign to new dataframes
    df_dark_fc_phic = f2mtz_fun(df_dark_fmodel)

    # next separate the light and dark observable columns
    # I am keeping the name of the columns as in the original script
    df_dark_ampli = df_dark_fmodel[["FOBS", "SIGFOBS"]].copy()
    columns_to_rename_dark = df_dark_ampli.columns
    df_dark_ampli = df_dark_ampli.rename(
        columns={
            columns_to_rename_dark[0]: "F_DARK",
            columns_to_rename_dark[1]: "SIGF_DARK",
        }
    )

    df_light_ampli = df_light_fmodel[["FOBS", "SIGFOBS"]].copy()
    columns_to_rename_light = df_light_ampli.columns
    df_light_ampli = df_light_ampli.rename(
        columns={
            columns_to_rename_light[0]: "F_LIGHT",
            columns_to_rename_light[1]: "SIGF_LIGHT",
        }
    )

    # join dataframes together, arrange columns in a certain order
    df_all = cad_fun(df_dark_fc_phic, df_dark_ampli, df_light_ampli)
    df_all.compute_dHKL(inplace=True)
    if res_cutoff_max != None:
        df_all = df_all[df_all.RESOLUTION > res_cutoff_max]
    if res_cutoff_min != None:
        df_all = df_all[df_all.RESOLUTION < res_cutoff_min]

    F_DARK_scaled, SIGF_DARK_scaled, F_LIGHT_scaled, SIGF_LIGHT_scaled = calculate_scaled_param(args.scaling, df_all)

    cell_units = df_all.cell
    spacegroup_input = df_all.spacegroup

    # WEIGHTED DIFFERENCES
    # change the phases
    df_all.loc[df_all["PHI_D"] > 180.0, "PHI_D"] = df_all["PHI_D"] - 360.0
    df_all.loc[df_all["PHI_D"] < 180.0, "PHI_D"] = df_all["PHI_D"] + 360.0
    df_all.loc[(df_all["PHI_D"] > 180.0), "PHI_D"] = df_all["PHI_D"] - 360.0

    calculate_weights(
        args.w_method,
        BDC,
        df_all,
        F_LIGHT_scaled,
        F_DARK_scaled,
        SIGF_LIGHT_scaled,
        SIGF_DARK_scaled,
        cell_units, 
        spacegroup_input
    )
    print("----- Finished! -----")


if __name__ == "__main__":
    main()
