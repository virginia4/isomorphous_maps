import reciprocalspaceship as rs
import pandas as pd
import numpy as np

def compute_weights_Marius(df, dfsqm, sigdf, sigdfsqm):
    print("Use Marius weights")
    wm = 1 + (df**2 / dfsqm) + (sigdf / sigdfsqm)
    return wm**-1

def compute_weights_Zhong(df, adfm, sigdf, sigdfsqm):
    print("Use Zhong weights")
    wz = 1 + (df**2 / adfm**2) + (sigdf / sigdfsqm**2)
    return wz**-1

def no_weights(df_all, BDC, F_LIGHT_scaled, F_DARK_scaled, 
               SIGMF_LIGHT_scaled, SIGMF_DARK_scaled, 
               cell_units, spacegroup_input
               ):
    DF = F_LIGHT_scaled - BDC*F_DARK_scaled
    S12 = SIGMF_LIGHT_scaled**2 + SIGMF_DARK_scaled**2

    col1 = rs.DataSeries(DF, dtype="F")
    col2 = rs.DataSeries(S12, dtype="Q")
    col3 = rs.DataSeries(df_all["PHI_D"], dtype="P")

    df_all = pd.concat([col1, col2, col3], axis=1)
    df_all = df_all.rename(columns={0: "Diff_OBS", 1: "SIGDF", "PHI_D": "PHI"})
    df_all = df_all.dropna()
    df_all.loc[df_all["Diff_OBS"] < 0, "PHI"] = df_all["PHI"] + 180.0
    df_all.loc[(df_all["Diff_OBS"] < 0) & (df_all["PHI"] > 180.0), "PHI"] = (
        df_all["PHI"] - 360.0
    )
    df_all.loc[df_all["Diff_OBS"] < 0, "Diff_OBS"] = np.abs(df_all.Diff_OBS)
    df_all_to_mtz = rs.DataSet(
        df_all, cell=cell_units, spacegroup=spacegroup_input
    ).infer_mtz_dtypes()

    return df_all_to_mtz


def q_weights(df_all, F_LIGHT_scaled, F_DARK_scaled, 
              SIGMF_LIGHT_scaled, SIGMF_DARK_scaled, 
              cell_units, spacegroup_input
              ):
    DFSQ = (F_LIGHT_scaled - F_DARK_scaled) ** 2
    S12SQ = SIGMF_LIGHT_scaled**2 + SIGMF_DARK_scaled**2
    DFSQM = DFSQ.mean()
    S12SQM = S12SQ.mean()

    ## DETERMINE MEAN WEIGHT TO PRESERVE ABSOLUTE SCALE
    DF = F_LIGHT_scaled - F_DARK_scaled
    S12 = SIGMF_LIGHT_scaled**2 + SIGMF_DARK_scaled**2

    W = compute_weights_Marius(DF, DFSQM, S12, S12SQM)
    WMEAN = W.mean()
    print("WMEAN", WMEAN)
    DFDWM = DF.mean() / WMEAN
    print("MARIUS WEIGHTS")
    print("MEAN AMPLITUDE DIFFERENCE/<WEIGHT>, (M)", DFDWM)
    print("MEAN SQUARE AMPLITUDE DIFFERENCE (M)", DFSQM)
    print("MEAN SQUARE SIGMA OF DIFFERENCE AMPLITUDES (M):", S12SQM)

    col1 = rs.DataSeries(DF, dtype="F")
    col2 = rs.DataSeries(S12, dtype="Q")
    col3 = rs.DataSeries(W, dtype="W")
    col4 = rs.DataSeries(df_all["PHI_D"], dtype="P")

    df_all = pd.concat([col1, col2, col3, col4], axis=1)
    df_all = df_all.rename(
        columns={0: "Diff_OBS", 1: "SIGDF", 2: "FOM", "PHI_D": "PHI"}
    )
    df_all = df_all.dropna()
    df_all.loc[df_all["Diff_OBS"] < 0, "PHI"] = df_all["PHI"] + 180.0
    df_all.loc[(df_all["Diff_OBS"] < 0) & (df_all["PHI"] > 180.0), "PHI"] = (
        df_all["PHI"] - 360.0
    )
    df_all.loc[df_all["Diff_OBS"] < 0, "Diff_OBS"] = np.abs(df_all.Diff_OBS)

    df_all_to_mtz = rs.DataSet(
        df_all, cell=cell_units, spacegroup=spacegroup_input
    ).infer_mtz_dtypes()

    return df_all_to_mtz


def Z_weights(df_all, F_LIGHT_scaled, F_DARK_scaled, 
              SIGMF_LIGHT_scaled, SIGMF_DARK_scaled, 
              cell_units, spacegroup_input
              ):
    DF = F_LIGHT_scaled - F_DARK_scaled
    ADF = np.abs(F_LIGHT_scaled - F_DARK_scaled)
    S12 = np.sqrt(SIGMF_LIGHT_scaled**2 + SIGMF_DARK_scaled**2)

    DFM = DF.mean()
    S12M = S12.mean()
    ADFM = ADF.mean()
    WZ = compute_weights_Zhong(DF, ADFM, S12, S12M)
    WZMEAN = WZ.mean()

    print("ZHONG WEIGHTS")
    print("WZMEAN", WZMEAN)
    print("MEAN AMPLITUDE DIFFERENCE", DFM)
    print("MEAN AMPLITUDE DIFFERENCE SQUARED", DFM**2)
    print("MEAN ABSOLUTE AMPLITUDE DIFFERENCE", ADFM)
    print("MEAN ABSOLUTE AMPLITUDE DIFFERENCE SQUARED (Z)", ADFM**2)
    print("MEAN SIGMA OF DIFFERENCE AMPLITUDES", S12M)
    print("MEAN SIGMA OF DIFFERENCE AMPLITUDES SQUARED (Z)", S12M**2)

    #        col1 = rs.DataSeries(DF / WZMEAN, dtype="F")
    col1 = rs.DataSeries(DF, dtype="F")
    col2 = rs.DataSeries(S12, dtype="Q")
    col3 = rs.DataSeries(WZ, dtype="W")
    col4 = rs.DataSeries(df_all["PHI_D"], dtype="P")

    df_all = pd.concat([col1, col2, col3, col4], axis=1)
    df_all = df_all.rename(
        columns={0: "Diff_OBS", 1: "SIGDF", 2: "FOM", "PHI_D": "PHI"}
    )
    df_all = df_all.dropna()
    df_all.loc[df_all["Diff_OBS"] < 0, "PHI"] = df_all["PHI"] + 180.0
    df_all.loc[(df_all["Diff_OBS"] < 0) & (df_all["PHI"] > 180.0), "PHI"] = (
        df_all["PHI"] - 360.0
    )
    df_all.loc[df_all["Diff_OBS"] < 0, "Diff_OBS"] = np.abs(df_all.Diff_OBS)

    df_all_to_mtz = rs.DataSet(
        df_all, cell=cell_units, spacegroup=spacegroup_input
    ).infer_mtz_dtypes()

    df_all_to_mtz.to_csv(
        "isomorphous_maps_ZHONG.phs", header=None, sep="\t", float_format="%11.4f"
    )

    return df_all_to_mtz


def calculate_weights(method, BDC, df_all, F_LIGHT_scaled, 
                      F_DARK_scaled, SIGMF_LIGHT_scaled, 
                      SIGMF_DARK_scaled, cell_units, 
                      spacegroup_input
    ):
    if method == "no_weights":
        print("No weights")
        df_all_no_weights = no_weights(df_all, BDC, F_LIGHT_scaled, F_DARK_scaled, SIGMF_LIGHT_scaled, SIGMF_DARK_scaled, cell_units, spacegroup_input)
        if BDC == 1.0: 
            print("Saving different maps...")
            df_all_no_weights.write_mtz("isomorphous_maps_NOWEIGHT.mtz")
        else: 
            print("Saving different maps...")
            print("Background subtraction for BDC factor = ", BDC)
            df_all_no_weights.write_mtz("isomorphous_maps_NOWEIGHT_wBDC.mtz")

    elif method == "q_weights":
        print("Marius weights")
        ## for Marius' weights
        df_all_q_weights = q_weights(df_all, F_LIGHT_scaled, F_DARK_scaled, SIGMF_LIGHT_scaled, SIGMF_DARK_scaled, cell_units, spacegroup_input)
        print("Saving different q-weighted maps...")
        df_all_q_weights.write_mtz("isomorphous_maps_qweights.mtz")

    elif method == "Zhong_weights":
        print("Zhong weights")
        ## for Zhong' weights
        df_all_Z_weights = Z_weights(df_all, F_LIGHT_scaled, F_DARK_scaled, SIGMF_LIGHT_scaled, SIGMF_DARK_scaled, cell_units, spacegroup_input)
        print("Saving different Zhong-weighted maps...")
        df_all_Z_weights.write_mtz("isomorphous_maps_ZHONG.mtz")



        
