#Import packages
import pool.hdf5_io as hio
import pool.paths as paths
import pool.data_params as data_params
import pandas as pd
import numpy as np
import os
from pool.preprocessing import load_lesions_and_hemis

# Load list of included patients
listids= np.loadtxt(os.path.join(paths.data_dir,'patients.txt'),dtype=str)

#Load lesions and lesional hemisphere from included patients
# You can smooth the lesion labels

lesions,lesional,hemis=load_lesions_and_hemis(listids,smoothing=3)

np.savez_compressed(os.path.join(paths.data_dir,'lesions.npz'), lesions, dtype=bool)

#Loads demographics csv and create column lesion_masked which identifies which patients have lesion masks
df=pd.read_csv(os.path.join(paths.data_dir,'demographics.csv'))
df['lesion_masked'] = lesional
df.to_csv(os.path.join(paths.data_dir,'demographics.csv'),index=False)


