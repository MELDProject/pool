##Generate spins for spatial permutation testing, used to compare 2 maps


import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
import os

from scipy.stats import special_ortho_group
import time 
from scipy.spatial import cKDTree
import pool.paths as paths
import pool.mesh_tools as mesh_tools
base_dir=paths.BASE_PATH
#load lh sphere of fsaverage_sym
sphere = nb.freesurfer.io.read_geometry(base_dir+'fsaverage_sym/surf/lh.sphere')

#Load cortex label
cortex = np.sort(nb.freesurfer.read_label(os.path.join(base_dir,"fsaverage_sym/label/lh.cortex.label")))

coords=sphere[0]
coords=coords[cortex]
tree=cKDTree(coords)
n_perm=1000
indices=np.zeros((n_perm,len(coords))).astype('int16')
for x in np.arange(n_perm):
    rotation = special_ortho_group.rvs(3)
    new_coords = coords @ rotation
    distance, indices[x]=tree.query(new_coords,k=1)
    
np.save(os.path.join(paths.data_dir,'spins_1000.npy',indices))