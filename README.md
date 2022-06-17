# isomorphous_maps


## Requirements: 
* install cctbx library. Go to http://cci.lbl.gov/docs/cctbx/installation/ for directions
* need to python3.8 to run 

## Input files: 
* PHI_file: a PBD file from which the script will extract the phases to use in the calculation of the maps 
* F_file1: reflection data 1 (in our case from the **dark** state of the protein)
* F_file2: reflection data 2 (here: data after the protein crystals were excited using some **light** source)

## Example **latest version** : 
python3.8 VA_weights.py --PHI_file=photolyase_refine_87.pdb --F_file1=pl-100ps_partial-off.mtz --F_file2=pl-100ps_partial-on.mtz --w_method=Zhong_weights --backgr_correction=0.9 --res_max=.. --res_min=.. --scaling=aniso

## Key input arguments for --w_mehtod and --scaling
* Options for --w_methods: no_weights, q_weights, Zhong_weights
* Options for --scaling: linear, iso, aniso


## The script output three files **latest version** : 
* crystal_data_1_fmodel.mtz: a mtz file containing observed and calculated amplitudes and phases. Calculated based on the data on the F_file1 (in the example the dark data)
* crystal_data_2_fmodel.mtz: a mtz file containing observed and calculated amplitudes and phases. Calculated based on the data on the F_file2 (in the example the light data)
* isomorphous_maps_[name of ].mtz: output file containing the isomorphous maps: differences, phases, weights [Diff_OBS, PHI, FOM]. Use coot or other similar software for visualization 

## Current issues and future work
* The reflection data that I have tried so far are for IMEAN i.e. intensities, which I convert to amplitudes using the function for ccbtx. Maybe if one uses reflections data with amplitudes will cause a problem. I need to double check

