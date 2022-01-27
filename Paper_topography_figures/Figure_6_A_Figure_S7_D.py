### Figure 6 A and Figure S7 D - Obenhaus et al. 
# Aligned FOV projection = Topographic score maps 


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
from helpers_topography.plotting_helpers import plot_summary_scoremaps
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

    all_params_ = ScoremapFOVProjParams.fetch(as_dict=True)
    print(tabulate(all_params_, headers='keys'))





    #### PLOT SCOREMAP EXAMPLES #######################################################################

    animal_name = '88592'
    binning_parameters = ['grid','border_bvs','mvl','info_content'] # ovc # info_content
    proj_param_id = 'B' #ScoremapFOVProjParams() - printed above
    top_percentile_cmap = None
    scalebar_length = 50. # microns
    param_hash_id_cell ='standard'
    param_hash_session ='cf83e1357eefb8bd'
    # Draw user annotated mask only? 
    anat_mask_label = None # pas_label or mec_label or None (if None, no mask will be drawn)
    anat_mask_sigma = 0.




    # How many unique reference metasessions are found in AlignmentFOV - there may only be one right now
    meta_ref_animal = set((AlignmentFOV * Session.proj(..., metasession_ref='metasession_name') 
                           & f'animal_name = "{animal_name}"').fetch('metasession_ref'))
    print(f'Found {len(meta_ref_animal)} reference metasession(s) for animal {animal_name}')
    if len(meta_ref_animal) != 1:
        print('None / too many reference metasessions found. Stopping.')
        sys.exit()

    # Check binning 
    bins_sizes_scoremap_params = ScoremapFOVParams.fetch('bin_size_microns')
    if len(set(bins_sizes_scoremap_params)) > 1:
        print(f'More than one bin size found ({set(bins_sizes_scoremap_params)}). Stopping.')
        sys.exit()
    else:
        bins_size = bins_sizes_scoremap_params[0]

    # ONLY FOR PRINTING OUT INFORMATION:
    statistic_, sigma_ = ((ScoremapFOVProjParams 
                          & f'proj_param_id = "{proj_param_id}"').fetch1('proj_statistic','proj_sigma'))


    print('\n\nPlotting FOV scoremap collage')
    print(f'Choosing projection statistic: {statistic_}')
    print(f'and projection sigma: {sigma_}\n')

    plot_summary_scoremaps(animal_name, 
                           binning_parameters=binning_parameters, 
                           proj_param_id=proj_param_id,
                           metasession_ref=None, 
                           param_hash_id_cell=param_hash_id_cell,
                           param_hash_session=param_hash_session,
                           percentile=top_percentile_cmap,
                           save_path=figure_dir,
                           bin_size_microns=bins_size,
                           scalebar_length=scalebar_length,
                           anat_mask_label=anat_mask_label,
                           anat_mask_sigma=anat_mask_sigma,
                           )




    print(figure_dir)
    print('Success.')