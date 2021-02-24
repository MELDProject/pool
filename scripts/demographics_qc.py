#Import packages
import pandas as pd
import os
import pool.paths as paths
import matplotlib.pyplot as plt
import seaborn as sns
import ptitprince as pt
import numpy as np

#Load demographics csv file
df=pd.read_csv(os.path.join(paths.data_dir,'demographics.csv'))

#problem - age of onset > age at preoperative
problem_2 = df['Age at preoperative'] <            df['Age of onset']
#swap round data for two patients in H26, where entered age of onset and age at preop are in the wrong columns
swap_columns=['MELD_H26_15T_FCD_0010','MELD_H26_3T_FCD_0011']
for ids in df[problem_2]["ID"]:
    if ids in swap_columns:
        df['Age of onset'][df["ID"]==ids] = 9.0
        df['Age at preoperative'][df["ID"]==ids] = 12.0
    else:
        df['Age of onset'][df["ID"]==ids]=np.nan
        df['Age at preoperative'][df["ID"]==ids] = np.nan

#replace recorded duration with difference between onset and age at scan
df['Duration'] = df['Age at preoperative'] -            df['Age of onset']

#Save qc'd demographics csv as demographics_qc.csv
df.to_csv(os.path.join(paths.data_dir,'demographics_qc.csv'),index=False)

