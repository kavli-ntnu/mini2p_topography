# Cell maps in Fig. 5B and Fig. S6B - Obenhaus et al. 

import sys, os
import os.path 
import numpy as np 
import pandas as pd 
import datajoint as dj
import cmasher as cmr

# Make plots pretty 
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='white')

# Prevent bug in figure export as pdf: 
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

##### IMPORTS ###########################################################################
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from Figure_3_A_and_B import plot_cellmaps

##### LOAD SCHEMA COMPONENTS #############################################################
from dj_schemas.dj_conn import * 


##### EXPORT LOCATION ####################################################################
if __name__ == "__main__":
  
    #session_name = '7e888f1d8eaab46b' # (animal 88592) Grid
    session_name  = 'dd100b60c09bfdbe' # (animal 88106) OV or 0ed19ecd643fdafa
    
    # Cell maps
    param_hash_id_cell = 'standard' # FilteredCells ! 
    plot_cellmaps(session_name, 
                  param_hash_id_cell=param_hash_id_cell, 
                  plot_ov=True,     
                  draw_pixels=False, 
                  draw_centers=True, 
                  dot_size=160
                  )

    print('Success.')