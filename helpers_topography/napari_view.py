### Helpers for viewing sessions with Napari 

from datetime import datetime

import numpy as np 
import napari # www.napari.org

import sys
sys.path.append('..')
from dj_schemas.utils import session_title_string # Create title in Napari viewer



def label_anatomy_napari(projection_corr_data, anatomical_mask_dict, session_name):
    '''
    Create a Napari session and display Max image projections, etc... 
    Create a label layer that enables user to paint in anatomical boundaries. 
    What the "label" (ID) corresponds to anatomically is set by the entries in anatomical_mask_dict 
    [retrieved from AnatomicalMaskParams()].
    
    The function returns the viewer object itself (so masks, etc. can be extracted) and 
    the name of the label layer (so that the layers can be referenced).
    '''
    
    
    anatomical_mask_dict = {k: v for k, v in anatomical_mask_dict.items() if v is not None}

    # Make a label layer name, consisting of anatomical keywords 
    label_layer_name = 'MEC{},PAS{},RSA{},PRH{}'.format(anatomical_mask_dict.get('mec_label','-'), \
                                                  anatomical_mask_dict.get('pas_label','-'),\
                                                  anatomical_mask_dict.get('rsa_label','-'),\
                                                  anatomical_mask_dict.get('prh_label','-'))
    with napari.gui_qt():
        viewer = napari.Viewer()
        viewer.window.resize(800,850)
        viewer.add_image(projection_corr_data['max_img_corr'], rgb=False, name='[Ch 1] Max image', opacity=1, colormap='green', gamma=.5)
        viewer.add_image(projection_corr_data['mean_image_corr'], rgb=False, name='[Ch 1] Mean image', colormap='green', visible=False, blending='additive')
        viewer.add_image(projection_corr_data['vmap_image_corr'], rgb=False, name='[Ch 1] Vmap image', colormap='green', visible=False, blending='additive')
        if not np.isnan(projection_corr_data['mean_image_second_corr']).all():
            viewer.add_image(projection_corr_data['mean_image_second_corr'], rgb=False, name='[Ch 2] image', opacity=1, colormap='red', gamma=.85, blending='additive')

            # Increase contrast of mean_image_second_corr:
            viewer.layers['[Ch 2] image'].contrast_limits_range = [0,18000]
            viewer.layers['[Ch 2] image'].contrast_limits       = [0,12000]

        # Add empty label layer
        viewer.add_labels(np.zeros_like(projection_corr_data['mean_image_corr']), name=label_layer_name)

        # Scale and translate 
        for layer_no in range(len(viewer.layers)):
            viewer.layers[layer_no].scale = [1.8,1.8]
            viewer.layers[layer_no].translate = [-240,-10]
            
        # Set up brush layer to add anatomical labels 
        viewer.layers[label_layer_name].mode = 'paint'
        viewer.layers[label_layer_name].selected_label = 8
        viewer.layers[label_layer_name].opacity = .2
        viewer.layers[label_layer_name].brush_size = 35
        
        viewer.title = session_title_string(session_name)

    return viewer, label_layer_name