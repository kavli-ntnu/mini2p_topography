### Figure 6B - Obenhaus et al. 
# Cross correlation of
# topographic tuning maps

import sys, os
import os.path
from pathlib import Path
import numpy as np 
import datajoint as dj
import cmasher as cmr

# Make plots pretty 
import seaborn as sns
sns.set(style='white')

# Prevent bug in figure export as pdf: 
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

##### IMPORTS #########################################################################################
from matplotlib import pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from helpers_topography.plotting_helpers import plot_summary_xcorr_scoremaps
from tabulate import tabulate

##### LOAD SCHEMA COMPONENTS ########################################################################## 
from dj_schemas.dj_conn import *

##### EXPORT LOCATION #################################################################################
figure_dir = 'YOUR_EXPORT_DIRECTORY/'
figure_dir = Path(figure_dir)





if __name__ == "__main__":

    #### PRINT AVAILABLE PROJECTION PARAMETERS 
    # These are now important,
    # since they come into play for the cross correlation calculations

    all_params = ScoremapFOVProjParams.fetch(as_dict=True)
    print(tabulate(all_params, headers='keys'))



    #### PLOT XCORR FOV EXAMPLES #######################################################################################################################

    animal_name = '88592'
    binning_parameters = ['grid','border_bvs','mvl'] # ovc, info_content
    proj_param_id = 'A' #ScoremapFOVProjParams() - printed above
    scalebar_length = 50. # microns
    param_hash_id_cell ='standard'
    param_hash_session ='cf83e1357eefb8bd'




    # ONLY FOR PRINTING OUT INFORMATION:
    statistic_, sigma_ = ((ScoremapFOVProjParams 
                          & f'proj_param_id = "{proj_param_id}"').fetch1('proj_statistic','proj_sigma'))
    print(f'\n\nPlotting FOV XCORR for {animal_name}')
    print(f'Choosing projection statistic: {statistic_}')
    print(f'and projection sigma: {sigma_}\n')

    # Check binning 
    bins_sizes_scoremap_params = ScoremapFOVParams.fetch('bin_size_microns')
    if len(set(bins_sizes_scoremap_params)) > 1:
        print(f'More than one bin size found ({set(bins_sizes_scoremap_params)}). Stopping.')
        sys.exit()
    else:
        bins_size = bins_sizes_scoremap_params[0]



    ######## FILTER SESSIONS / ENTRIES ################################################################################################################
    sessions_animal = Session.proj(animal_name= 'animal_name', metasession_ref='metasession_name') & f'animal_name = {animal_name}' 

    scoremap_corr_binning = (ScoremapCorr 
                            * ScoremapFOVParams.proj(binning_param_A='binning_param_short', binning_param_set_A='binning_param_set') 
                            * ScoremapFOVParams.proj(binning_param_B='binning_param_short', binning_param_set_B='binning_param_set')
                            )

    scoremap_corr_entries_animal = (scoremap_corr_binning & [f'binning_param_set_A = "{param}"' for param in binning_parameters] 
                                                          & [f'binning_param_set_B = "{param}"' for param in binning_parameters]
                                                          & f'proj_param_id = "{proj_param_id}"'
                                                          & f'param_hash_id_cell= "{param_hash_id_cell}"'
                                                          & f'param_hash_session= "{param_hash_session}"'
                                                          & sessions_animal)

    ######## PLOT ALL PAIRWISE COMBINATIONS ############################################################################################################

    for entry in scoremap_corr_entries_animal:
        plot_summary_xcorr_scoremaps(animal_name,
                                    entry,
                                    invert_order=True,
                                    save_path=figure_dir,
                                    bin_size_microns=bins_size,
                                    scalebar_length=scalebar_length,
                                    )



    print(figure_dir)
    print('Success.')