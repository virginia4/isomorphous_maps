from scipy.optimize import curve_fit, least_squares
import numpy as np

def f2mtz_fun(df_in):
    """
    f2mtz_fun just need to
    pick the columns we need for the calculations and assign them in a new Dataframe.
    It takes the 3 columns as inputs: FC, SIGFC and PHIC.
    When there is not a column SIGC in the .phs file, it creates the file and fills
    it with 1s (this is just to be consistent with what the f2mtz does)
    """
    SIGC_in_dataframe = "SIGF_C" in df_in
    if SIGC_in_dataframe == True:
        df_out = df_in[["FCALC", "SIGFCALC", "PHIFCALC", "RESOLUTION"]].copy()
    else:
        df_out = df_in[["FCALC", "PHIFCALC", "RESOLUTION"]].copy()
        #     dataseries_SIGFCAL = rs.DataSeries(1, dtype="Q")
        df_out["SIGFCALC"] = 1
        df_out["SIGFCALC"] = df_out["SIGFCALC"].astype("Stddev")
        df_out = df_out.loc[:, ["FCALC", "SIGFCALC", "PHIFCALC", "RESOLUTION"]]
        return df_out


def cad_fun(df1, df2, df3):
    df1 = df1.sort_index()
    df2 = df2.sort_index()
    df3 = df3.sort_index()

    df_12 = df1.join(df2)
    df_all = df_12.join(df3)
    df_all = df_all.dropna()

    columns_to_rename = df_all.columns
    # rename the columns of df_all
    df_all = df_all.rename(
        columns={
            columns_to_rename[0]: "FC",
            columns_to_rename[1]: "SIGF_C",
            columns_to_rename[2]: "PHI_D",
            columns_to_rename[3]: "RESOLUTION",
            columns_to_rename[4]: "F_DARK",
            columns_to_rename[5]: "SIGF_DARK",
            columns_to_rename[6]: "F_LIGHT",
            columns_to_rename[7]: "SIGF_LIGHT",
        }
    )
    return df_all


# SCALING FUNCTIONS
def lin_scale_fun(x, a):
    return a * x

def iso_scale_func(p, x1, x2, qs):
    r = x1 - (p[0]*np.exp(-p[1]*(qs**2)))*x2
    return r

def aniso_scale_func(p, x1, x2, H_arr):
    
    h = H_arr[:,0]
    k = H_arr[:,1]
    l = H_arr[:,2]
    
    h_sq = np.square(h)
    k_sq = np.square(k)
    l_sq = np.square(l)
    
    hk_prod = h*k
    hl_prod = h*l
    kl_prod = k*l

    t = - (h_sq * p[1] + k_sq * p[2] + l_sq * p[3] + 
        2*hk_prod * p[4] + 2*hl_prod * p[5] + 2*kl_prod * p[6])    
    expnt = np.exp( t )
    r = x1 - p[0] * expnt * x2  
    return r

def scale_lin(x_dataset, y_dataset):
    # scales y_dataset to x_dataset using linear scaling

    const, _ = curve_fit(lin_scale_fun, x_dataset, y_dataset)

    data_lin_scaled = y_dataset/const.item()
    return const.item(), data_lin_scaled

def scale_iso(x_dataset, y_dataset, indx):

    """
    AUTHOR: Alisia Fadini
    Isotropic resolution-dependent scaling of data2 to data1.
    (minimize [dataset1 - c*exp(-B*sintheta**2/lambda**2)*dataset2]
    Input :
    1. dataset1 in form of 1D numpy array
    2. dataset2 in form of 1D numpy array
    3. dHKLs for the datasets in form of 1D numpy array
    Returns :
    1. entire results from least squares fitting
    2. c (as float)
    3. B (as float)
    2. scaled dataset2 in the form of a 1D numpy array
    """

    # # scale observed data (for F_file_1 and F_file_2) against model data (F_C)
    # const1, _ = curve_fit(model_fun, df_all.F_DARK, df_all.FC)
    # const2, _ = curve_fit(model_fun, df_all.F_LIGHT, df_all.FC)

    # F_DARK_scaled = model_fun(df_all.F_DARK, const1.item())
    # SIGF_DARK_scaled = model_fun(df_all.SIGF_DARK, const1.item())
    # F_LIGHT_scaled = model_fun(df_all.F_LIGHT, const2.item())
    # SIGF_LIGHT_scaled = model_fun(df_all.SIGF_LIGHT, const2.item())

        
    # def iso_scale_func(p, x1, x2, qs):
    #     r = x1 - (p[0]*np.exp(-p[1]*(qs**2)))*x2
    #     return r
    
    p0 = np.array([1.0, -20])
    qs = 1/(2*indx)
    matrix = least_squares(iso_scale_func, p0, args=(x_dataset, y_dataset, qs))
    
    return matrix, matrix.x[0], matrix.x[1], (matrix.x[0]*np.exp(-matrix.x[1]*(qs**2)))*y_dataset


def scale_aniso(x_dataset, y_dataset, Miller_indx):
    
    p0 = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float32)
    matrix_ani = least_squares(aniso_scale_func, p0, args=(x_dataset, y_dataset, Miller_indx))

    h = Miller_indx[:,0]
    k = Miller_indx[:,1]
    l = Miller_indx[:,2]
    h_sq = np.square(h)
    k_sq = np.square(k)
    l_sq = np.square(l)

    hk_prod = h*k
    hl_prod = h*l
    kl_prod = k*l

    t = - (h_sq * matrix_ani.x[1] + k_sq * matrix_ani.x[2] + l_sq * matrix_ani.x[3] 
       + 2*hk_prod * matrix_ani.x[4] + 2*hl_prod * matrix_ani.x[5] + 2*kl_prod * matrix_ani.x[6])

    data_ani_scaled = (matrix_ani.x[0]*np.exp(t))*y_dataset
    
    return matrix_ani, data_ani_scaled

def calculate_scaled_param(scaling_method, df_all):
# scale observed data (for F_file_1 and F_file_2) against model data (F_C)
# linear scaling

    if scaling_method == 'linear':
        print('Linear scaling')
        lin_scl_factor_dark, F_DARK_lin_scaled = scale_lin(df_all.FC, df_all.F_DARK)
        SIGF_DARK_lin_scaled = df_all.SIGF_DARK/lin_scl_factor_dark

        lin_scl_factor_light, F_LIGHT_lin_scaled = scale_lin(df_all.FC, df_all.F_LIGHT)
        SIGF_LIGHT_lin_scaled = df_all.SIGF_LIGHT/lin_scl_factor_light

        return F_DARK_lin_scaled, SIGF_DARK_lin_scaled, F_LIGHT_lin_scaled, SIGF_LIGHT_lin_scaled

    elif scaling_method == 'iso':
        print('--- Isotropic scaling ---')
        res_matrix = df_all.dHKL.to_numpy(dtype=np.float32)
        qs = 1/(2*res_matrix)
        iso_scl_factors_dark, iso_scl_factors_dark_1, iso_scl_factors_dark_2, F_DARK_iso_scaled = scale_iso(df_all.FC, df_all.F_DARK, res_matrix)
        SIGF_DARK_iso_scaled = (iso_scl_factors_dark_1 * np.exp(-iso_scl_factors_dark_2*(qs**2)))  * df_all.SIGF_DARK
        iso_scl_factors_light, iso_scl_factors_light_1, iso_scl_factors_light_2, F_LIGHT_iso_scaled = scale_iso(df_all.FC, df_all.F_LIGHT, res_matrix)
        SIGF_LIGHT_iso_scaled = (iso_scl_factors_light_1 * np.exp(-iso_scl_factors_light_2*(qs**2)))
        return F_DARK_iso_scaled, SIGF_DARK_iso_scaled, F_LIGHT_iso_scaled, SIGF_LIGHT_iso_scaled

    elif scaling_method == 'aniso':
        print('--- Anisotropic scaling ---')
        H = df_all.index.to_frame()
        H_arr = H.to_numpy(dtype=np.float32)

        h = H_arr[:,0]
        k = H_arr[:,1]
        l = H_arr[:,2]
        h_sq = np.square(h)
        k_sq = np.square(k)
        l_sq = np.square(l)

        hk_prod = h*k
        hl_prod = h*l
        kl_prod = k*l

        aniso_scl_factors_dark, F_DARK_aniso_scaled = scale_aniso(df_all.FC, df_all.F_DARK, H_arr)
        t_dark = - (h_sq * aniso_scl_factors_dark.x[1] + k_sq * aniso_scl_factors_dark.x[2] + l_sq * aniso_scl_factors_dark.x[3] 
           + 2*hk_prod * aniso_scl_factors_dark.x[4] + 2*hl_prod * aniso_scl_factors_dark.x[5] + 2*kl_prod * aniso_scl_factors_dark.x[6])
        SIGF_DARK_aniso_scaled = (aniso_scl_factors_dark.x[0]*np.exp(t_dark)) * df_all.SIGF_DARK

        aniso_scl_factors_light, F_LIGHT_aniso_scaled = scale_aniso(df_all.FC, df_all.F_LIGHT, H_arr)
        t_light = - (h_sq * aniso_scl_factors_light.x[1] + k_sq * aniso_scl_factors_light.x[2] + l_sq * aniso_scl_factors_light.x[3] 
           + 2*hk_prod * aniso_scl_factors_light.x[4] + 2*hl_prod * aniso_scl_factors_light.x[5] + 2*kl_prod * aniso_scl_factors_light.x[6])
        SIGF_LIGHT_aniso_scaled = (aniso_scl_factors_light.x[0]*np.exp(t_light)) * df_all.SIGF_LIGHT

        return F_DARK_aniso_scaled, SIGF_DARK_aniso_scaled, F_LIGHT_aniso_scaled, SIGF_LIGHT_aniso_scaled
