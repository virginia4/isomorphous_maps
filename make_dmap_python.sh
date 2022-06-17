#!/bin/tcsh

# source /Applications/ccp4-7.1/bin/ccp4.setup-sh


# edit here ===>
#
# signal in the map:  -0.19602 0.02184

## V: sets variables
## Four parameters: * dark model which should be the pdb file (dark_model)
##                  * Observable amplitudes from the dark model (mtz file) (dark_obs)
##                  * Observable amplitudes from the light model (light_obs)
##                  * model_tt.mtz structure factors, i.e. FCALC and PHICALC (model_F)
##                  * bin_nam (bin_nam)
set dark_model = ./dark_HCC_v2.pdb  
set dark_obs = ./FOBS_dark_V.mtz
set light_obs = ./FOBS_light_V.mtz
set model_F = ./model_tt_V.mtz
set bin_nam =  3ps



# pos 0.16 rmsd 0.02
## V: how is this parameters set?

set resmax = 1.6
set scalmin = 20.0
set mapmin = 20.0


## =============================

# phase file must be calculated from refined model in DARK
# folder
# dump the file and add column with 1.0 as sigmas 


# V: mtz2various (availble from CCP4) reads a MTZ file  and produces a ASCII file (assigned to HKLOUT) 
# in the suitable form. Parameter USER means that is accepts all inout columns. 
mtz2various HKLIN $model_F HKLOUT model_phs.hkl << mtz_phs
     LABIN FP=FC PHIC=PHIC
     OUTPUT USER '(3I5,F12.3,'  1.00  ',F12.3)' 
mtz_phs


# # get the phase file with added sig column into the system

# # V: f3mtz (availble from CCP4) convert a formatted reflection file to MTZ format
f2mtz HKLIN model_phs.hkl HKLOUT FC_dark.mtz << f2m_phs 
     ## SKIP 5  Miss out 5 lines of header
     CELL 66.9   66.9   40.9 90.0 90.0 120.0  # angles default to 90
     SYMM 173
     LABOUT H   K  L   FC_D SIG_FC_D PHI_D
     CTYPE  H   H  H   F     Q        P
f2m_phs


# # cad things together (collect and sort crystallographic reflection data from several files, to generate a single set.)

cad             \
HKLIN1 FC_dark.mtz HKLIN2 $dark_obs HKLIN3 $light_obs HKLOUT all.mtz << END-cad
     LABIN FILE 1 E1=FC_D E2=SIG_FC_D E3=PHI_D
     CTYP  FILE 1 E1=F E2=Q E3=P 
     LABIN FILE 2 E1=F_DARK E2=SIGF_DARK
     CTYP  FILE 2 E1=F E2=Q
     LABIN FILE 3 E1=F_LIGHT E2=SIGF_LIGHT
     CTYP  FILE 3 E1=F E2=Q
END
END-cad


# scale the things
# labels so far:
# H K L FC_D SIG_FC PHI_D F_DARK SIGF_DARK F_LIGHT SIGF_LIGHT
# 1 scale dark to FC dark
# 2 scale light to dark


echo " SCALEIT NUMBER 1, MARIUS CHECK "

## SCALEIT calculates and applies a derivative to native scaling 
## function using either (a) an overall scale factor, (b) a scale and 
## isotropic temperature factor or (c) a scale and anisotropic temperature factors. 
scaleit HKLIN all.mtz HKLOUT all_sc1.mtz << END-scaleit1
     TITLE FPHs scaled to FP
     reso $scalmin $resmax      # Usually better to exclude lowest resolution data
     #WEIGHT            Sigmas seem to be reliable, so use for weighting
     #refine anisotropic
     #Exclude FP data if: FP < 5*SIGFP & if FMAX > 1000000
     EXCLUDE FP SIG 4 FMAX 10000000
     REFINE ANISOTROPIC 
     #AUTO
     LABIN FP=FC_D SIGFP=SIGF_DARK  -
       FPH1=F_DARK SIGFPH1=SIGF_DARK -
       FPH2=F_LIGHT SIGFPH2=SIGF_LIGHT
     #LABIN FP=F_DARK SIGFP=SIGF_DARK  -
     #  FPH1=FC_D SIGFPH1=SIG_FC_D -
     #  FPH2=F_LIGHT SIGFPH2=SIGF_LIGHT
     CONV ABS 0.0001 TOLR  0.000000001 NCYC 150
END
END-scaleit1

echo " SCALEIT OVER, MARIUS CHECK "

scaleit HKLIN all_sc1.mtz HKLOUT all_sc2.mtz << END-scaleit2
     TITLE FPHs scaled to FP
     reso $scalmin $resmax      # Usually better to exclude lowest resolution data
     #WEIGHT    Sigmas seem to be reliable, so use for weighting
     #refine anisotropic
     #Exclude FP data if: FP < 5*SIGFP & if FMAX > 1000000
     REFINE ANISOTROPIC 
     EXCLUDE FP SIG 4 FMAX 10000000
     LABIN FP=F_DARK SIGFP=SIGF_DARK -
       FPH1=F_LIGHT SIGFPH1=SIGF_LIGHT
     CONV ABS 0.0001 TOLR  0.000000001 NCYC 40
END
END-scaleit2

# freerflag HKLIN all_sc2.mtz HKLOUT all_sc2_free.mtz << +freerfrac 0.05+
# freerflag HKLIN all_sc2.mtz HKLOUT all_sc2_free.mtz << END-freerflag
# END
# END-freerflag 

# echo "MARIUS unweighted maps"



fft HKLIN all_sc2.mtz MAPOUT ${bin_nam}_nonw.map << endfft
  RESO $mapmin  $resmax
  GRID 160 160 120
  BINMAPOUT
  LABI F1=F_LIGHT SIG1=SIGF_LIGHT F2=F_DARK SIG2=SIGF_DARK PHI=PHI_D
endfft


# # dump the scaled files to calculate the weighted map
# # -------------------------------------------------------------------------------
mtz2various HKLIN all_sc2.mtz  HKLOUT light_scaled.hkl << end_mtzv1
     LABIN FP=F_LIGHT SIGFP=SIGF_LIGHT
     OUTPUT USER '(3I5,2F12.3)'
     RESOLUTION 60.0 $resmax 
end_mtzv1

f2mtz hklin light_scaled.hkl hklout light_scaled.mtz << f2m_light
     CELL 66.9   66.9   40.8 90.0 90.0 120.0  # angles default to 90
     SYMM 173
     LABOUT H   K  L   F_LIGHT SIGF_LIGHT
     CTYPE  H   H  H   F      Q  
     END
f2m_light
# # -------------------------------------------------------------------------------


# # -------------------------------------------------------------------------------
mtz2various HKLIN all_sc2.mtz  HKLOUT dark_scaled.hkl << end_mtzv2
     LABIN FP=F_DARK SIGFP=SIGF_DARK
     OUTPUT USER '(3I5,2F12.3)'
     RESOLUTION 60.0 $resmax 
end_mtzv2

f2mtz hklin dark_scaled.hkl hklout dark_scaled.mtz << f2m_dark
     CELL 66.9   66.9   40.8 90.0 90.0 120.0  # angles default to 90
     SYMM 173
     LABOUT H   K  L   F_DARK SIGF_DARK
     CTYPE  H   H  H   F      Q  
     END
f2m_dark
# # -------------------------------------------------------------------------------



# # -------------------------------------------------------------------------------
mtz2various HKLIN all_sc2.mtz  HKLOUT dark_phase.hkl << end_mtzv3
     LABIN FP=FC_D SIGFP=SIG_FC_D PHIC=PHI_D
     OUTPUT USER '(3I5,3F12.3)'
     RESOLUTION 60.0 $resmax 
end_mtzv3

f2mtz hklin dark_phase.hkl hklout dark_phase.mtz << f2m_phase
     CELL 66.9   66.9   40.8 90.0 90.0 120.0  # angles default to 90
     SYMM 173
     LABOUT H   K  L   FC_D SIG_FC_D PHI_D
     CTYPE  H   H  H   F      Q      P
     END
f2m_phase
# # -------------------------------------------------------------------------------


# this is the wmar.inp file

echo light_scaled.mtz > wmar_py.inp
echo dark_scaled.mtz >> wmar_py.inp
echo dark_phase.mtz >> wmar_py.inp
echo ${bin_nam}_dark.mtz >> wmar_py.inp

# echo light_scaled.hkl > wmar.inp
# echo dark_scaled.hkl >> wmar.inp
# echo dark_phase.hkl >> wmar.inp
# echo ${bin_nam}_dark.phs >> wmar.inp

# # this will produce a difference structure factor file
# # h k l DF weight Phase
# # run weighting program
# # =======================>

# echo "Marius weighting"

# /Users/vapostolop/Desktop/scripts/Marius_scaling/test1/weight_zv2 < wmar.inp
echo "Marius weighting - python version"

python3 VA_weights.py wmar_py.inp

f2mtz HKLIN ${bin_nam}_dark_python.phs HKLOUT ${bin_nam}_dwt.mtz << end_weight 
CELL 66.9   66.9   40.8 90.0 90.0 120.0  # angles default to 90
SYMM 173
##  LABOUT H   K  L   DOBS_3ps  FOM_3ps  PHI
##  Change the line above to include the Z weights
LABOUT H   K  L   DOBS_3ps SIGDF  FOM_3ps  FOM_3ps_Z  PHI
CTYPE  H   H  H   F        Q      W        WZ         P
END
end_weight


#calculate weighted difference map

fft HKLIN ${bin_nam}_dwt.mtz MAPOUT ${bin_nam}_wd.map << END-wfft
  RESO $mapmin  $resmax 
  GRID 160 160 120 
  BINMAPOUT
  LABI F1=DOBS_3ps W=FOM_3ps PHI=PHI
END
END-wfft

fft HKLIN ${bin_nam}_dwt.mtz MAPOUT ${bin_nam}_wdz_python.map << END-wfft
  RESO $mapmin  $resmax 
  GRID 160 160 120 
  BINMAPOUT
  LABI F1=DOBS_3ps W=FOM_3ps_Z PHI=PHI
END
END-wfft

# mapmask mapin ${bin_nam}_wd.map mapout ${bin_nam}_wdex.map xyzin $dark_model << ee
# extend xtal
# border 0.0
# ee

# rm model_phs.hkl FC_dark.mtz
# rm light_dark.phs
# rm light_scaled.hkl dark_scaled.hkl dark_phase.hkl 
# rm all.mtz all_sc1.mtz all_sc2.mtz

