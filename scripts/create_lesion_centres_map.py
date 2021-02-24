import os
import nibabel as nb
import pandas as pd
import numpy as np
import pool.paths as paths
from pool.preprocessing import com_lesion_map


#import demographics and lesions
lesions=np.load(os.path.join(paths.data_dir,'lesions.npy'))
demographics = pd.read_csv(os.path.join(paths.data_dir,'demographics_qc.csv'))
demographics = demographics[demographics['lesion_masked'].astype(bool)]
lesions_com=com_lesion_map(demographics['ID'],lesions, lesion_size = np.mean(np.sum(lesions,axis=1)))

np.save(os.path.join(paths.data_dir,'lesions_coms.npy'), lesions_com)
