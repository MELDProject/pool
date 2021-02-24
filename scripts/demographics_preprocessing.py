# import necessary packages
import pool.hdf5_io as hio
import pool.paths as paths
import pool.data_params as data_params
import pandas as pd
import numpy as np
import os
from pool.preprocessing import tidy_features


#This creates a list of all the patient IDs (listids) and a list of IDs separated by scanner type (3T vs 1.5T)
listids, listids_scanner=hio.list_ids(data_params.site_codes,'patient')
data_path=paths.BASE_PATH

#Any patients to exclude
#patients confirmed as to be excluded by collaborating sites
listids.remove('MELD_H23_15T_FCD_0001')
listids.remove('MELD_H23_15T_FCD_0002')
listids.remove('MELD_H10_3T_FCD_0008')
#save out list of included patients into data folder
np.savetxt(os.path.join(paths.data_dir,'patients.txt'),listids,fmt='%s')

#loads in list of included patients
listids= np.loadtxt(os.path.join(paths.data_dir,'patients.txt'),dtype=str)

#load in features. Engel outcome will be replaced by seizure freedom further down
features=['ID','Age of onset','Duration','Age at preoperative','Sex','Ever reported MRI negative',
         'Engel Outcome','Surgery','f/u']

data=[]
for ids in listids:
    #create dictionaries of the features, then make into pandas dataframe
    #this retrieves the demographic features
    features_dict=dict(zip(features,hio.get_demographic_features(ids,features)))
    data.append(features_dict)

    
#This calculates the surface area of the lesion, the affected hemisphere and the affected lobes in each patient  
areas,hemis, lobes=hio.lesion_areas(listids)
#This retrieves the histology for each participant
histology=hio.histology_per_subject(listids)

#This creates a dataframe that included the demographic features, lesion area, hemisphere, lobe and histology
df=pd.DataFrame(data)
df['Lesion area']=areas
df['Hemisphere']=hemis
df['Lobe']=lobes
df['Histology']=histology

demographic_variables=list(df.columns)
#demographic_variables.remove('ID')
demographic_variables.remove('Engel Outcome')
demographic_variables.append('Seizure free')

# This tidies the features - e.g. replacing NO (not operated) with NaN
features_to_tidy=['Age of onset','Duration','Age at preoperative','Sex','Ever reported MRI negative',
         'Engel Outcome','Surgery','f/u']
for feature in features_to_tidy:
    print(feature)
    try:
        df[feature]=tidy_features(np.array(df[feature]))
    except ValueError:
      #  print(np.array(df[feature]))
        df[feature]=tidy_features(np.array(df[feature]).astype(float))

#switch engel for seizure free
df['Seizure free'] = (df['Engel Outcome']==1).astype(int)
df['Seizure free'][pd.isna(df['Engel Outcome'])]=np.nan

#remove 2.0 from surgery (some sites misscoded no surgery as 2 in the .csv files)
df['Surgery'][df['Surgery']==2.0]=0
df['Sex'][df['Sex']==2.0]=0
df['Ever reported MRI negative'][df['Ever reported MRI negative']==2.0]=0

df.to_csv(os.path.join(paths.data_dir,'demographics.csv'),index=False)

