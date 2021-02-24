import h5py
import os
import numpy as np
import pool.paths as paths
import pandas as pd
import nibabel as nb


def list_ids(site_codes, group='patient',base_path=paths.BASE_PATH):
    """Outputs a list of the patients included in the hdf5 file
    site_code eg H1, hospital site code
    group , patient , control or both"""
    if group == "both":
        groups = ['patient',  'control']
    else:
        groups = [group]
    #check if only one site or multiple
    if isinstance(site_codes, str):
        site_codes=[site_codes]
    subject_ids=[]
    subject_id_scanner={}
    for site_code in site_codes:
        for group in groups:
            if os.path.isfile(os.path.join(base_path,'MELD_'+site_code,site_code+"_"+group+"_featurematrix.hdf5")):
                f=h5py.File(os.path.join(base_path,'MELD_'+site_code,site_code+"_"+group+"_featurematrix.hdf5"),'r')
                scanners = f[site_code].keys()
                for scanner in scanners:
                    if scanner in subject_id_scanner.keys():
                        subject_id_scanner[scanner]+=list(f[os.path.join(site_code,scanner,group)].keys())
                    else:
                        subject_id_scanner[scanner]=list(f[os.path.join(site_code,scanner,group)].keys())
                    subject_ids+=list(f[os.path.join(site_code,scanner,group)].keys())
    return subject_ids, subject_id_scanner

  


def get_les_hemi(patient_id, base_path=paths.BASE_PATH,verbose=True):
    """Outputs which hemisphere is lesional when given a patient ID"""
    _,site_code,scanner,group,ID=patient_id.split('_')
    if group == 'C':
        print('Error: this is a control!')
        return 
    elif group == 'FCD':
        f=h5py.File(os.path.join(base_path,'MELD_'+site_code,site_code+"_patient_featurematrix.hdf5"),'r')
        keys=f[os.path.join(site_code,scanner+'/patient/'+patient_id+'/lh/')].keys()
        surf_dir=f.require_group(os.path.join(site_code,scanner,'patient',patient_id,'lh'))
        if '.on_lh.lesion.mgh' in surf_dir.keys():
            return 'lh'
        elif '.on_lh.lesion.mgh'in f.require_group(os.path.join(site_code,scanner,'patient',patient_id,'rh')) :
            return 'rh'
        else:
            if verbose:
        #        print('no lesion label found in patient ' + patient_id+' folders')
                return
            else:
                return('none')
    else:
        print('Error: incorrect naming scheme used. Unable to determine if patient or control.')
        return 
                 
def get_feature_list(patient_id, hemi='lh', base_path=paths.BASE_PATH):
    """Outputs a list of the features a participant has for each hemisphere"""
    
    _,site_code,scanner,group,ID=patient_id.split('_')
    if group =='C':
        groups ='control'
    else:
        groups='patient'
    f=h5py.File(os.path.join(base_path,'MELD_'+site_code,site_code+"_"+groups+"_featurematrix.hdf5"),'r')
    keys=list(f[os.path.join(site_code,scanner,groups,patient_id,hemi)].keys())
    if '.on_lh.lesion.mgh' in keys:
        keys.remove('.on_lh.lesion.mgh')
    return keys
                 
                 
def get_feature_values(patient_id, hemi='lh', feature='.on_lh.thickness.mgh', base_path=paths.BASE_PATH):
    """Outputs the values of a particular feature from a participant for one hemisphere"""
    _,site_code,scanner,group,ID=patient_id.split('_')
    if group =='C':
        groups ='control'
    else:
        groups='patient'
    f=h5py.File(os.path.join(base_path,'MELD_'+site_code,site_code+"_"+groups+"_featurematrix.hdf5"),'r')
    surf_dir=f[os.path.join(site_code,scanner,groups,patient_id,hemi)]
    overlay = surf_dir[feature][:]
    return overlay


def get_demographic_feature(patient_id,title,base_path=paths.BASE_PATH):
    """Outputs the demographic info of choice for a particular patient"""
    _,site_code,scanner,group,ID=patient_id.split('_')
    csv=pd.read_csv(os.path.join(base_path,'MELD_'+site_code,"MELD_"+site_code+"_participants.csv"), header=0,encoding='unicode_escape')
    Columns=csv.keys()
    new_title=[]
    for key in csv.keys():
        
        if title in key:
            new_title.append(key)
        elif 'ID' in key:
            ID_column=key
    if len(new_title)>1:
        print('Multiple columns matching string found, please make search more specific')
        return
    elif len(new_title)==0:
        print('Unable to find column matching this description. Please double check for typos')
        return
    else: 
        feature=csv[new_title][csv[ID_column]==patient_id].values
       
        return feature
    
def get_demographic_features(patient_id,titles,base_path=paths.BASE_PATH):
    """Outputs multiple demographic features of choice for a particular patient"""
    _,site_code,scanner,group,ID=patient_id.split('_')
    csv=pd.read_csv(os.path.join(base_path,'MELD_'+site_code,"MELD_"+site_code+"_participants.csv"), header=0,
                   encoding='unicode_escape')
    Columns=csv.keys()
    features=[]
    for title in titles:
        new_title=[]
        for key in csv.keys():
            if title in key:
                new_title.append(key)
            if 'ID' in key:
                ID_column=key
        if len(new_title)>1:
            print('Multiple columns matching string found, please make search more specific')
            return
        elif len(new_title)==0:
            print('Unable to find column matching this description. Please double check for typos')
            return
        else: 
            try :
                features.append(np.squeeze(csv[new_title][csv[ID_column]==patient_id].values).item())    
            except ValueError:
                pass

    return features
    


def get_csv_data(site_code,title,base_path=paths.BASE_PATH):
    """Outputs the columns of the csv file"""
    csv=pd.read_csv(os.path.join(base_path,'MELD_'+site_code,"MELD_"+site_code+"_participants.csv"), header=0, encoding='unicode_escape')
    Columns=csv.keys()
    new_title=[]
    for key in csv.keys():
        
        if title in key:
            new_title.append(key)
        elif 'ID' in key:
            ID_column=key
    if len(new_title)>1:
        print('Multiple columns matching string found, please make search more specific')
        return
    elif len(new_title)==0:
        print('Unable to find column matching this description. Please double check for typos')
        return
    else:   
        column=csv[new_title]
        csv_ids=csv[ID_column]
        return csv_ids, column
    
    
def get_subject_features(patient_id, features, hemi='lh',
                          base_path=paths.BASE_PATH):
    """Function to load all patient's data into memory
    Inputs: patient_id 
    features: list of features to be loaded
    hemi: 'lh' or 'rh' 
    Returns: feature_data, label
    """
    n_vert=163842
    _,site_code,scanner,group,ID=patient_id.split('_')
    # establish file name
    if group =='C':
        groups ='control'
    else:
        groups ='patient'
    f=h5py.File(os.path.join(base_path,'MELD_'+site_code,site_code+"_"+groups+"_featurematrix.hdf5"),'r')
    #load cortical mask if requested
    
    #establish if one hemisphere or both are needed
    feature_data=np.zeros((n_vert,len(features)))
    lesion_data=np.zeros(n_vert)
    surf_dir=f[os.path.join(site_code,scanner,groups,patient_id,hemi)]
        #check if lesional hemisphere
    if '.on_lh.lesion.mgh' in surf_dir.keys():
        lesion_data = np.ceil(surf_dir['.on_lh.lesion.mgh'][:]).astype(int)
    for i,feature in enumerate(features):
        #missing features set to zero
        if feature in surf_dir.keys():
            feature_data[:,i]=surf_dir[feature][:]
        else:
            pass
      #      print('missing feature : '+feature+' set to zero')
    return feature_data, lesion_data

from time import time
def load_subject_combined_hemisphere_data(patient_id, features, 
                          cortex=os.path.join(paths.BASE_PATH,'fsaverage_sym','label','lh.cortex.label'),
                         lesion_only=False,base_path=paths.BASE_PATH, normalise=False):
    """combine features from both hemispheres into single matrix and mask cortex and non-lesional data if necessary"""
    if isinstance(cortex,str) :
        cortex_label = nb.freesurfer.io.read_label(cortex)
    else:
        cortex_label=cortex
    for hemi in ['lh','rh']:
        hemi_data,hemi_label = get_subject_features(patient_id,features,hemi=hemi)
        #mask out non-cortical data
       # if cortex:
        hemi_data = hemi_data[cortex_label]
        hemi_label = hemi_label[cortex_label]
        #mask out healthy tissue from lesional hemisphere
        if lesion_only and sum(hemi_label)>0:
            hemi_data = hemi_data[hemi_label>0]
            hemi_label = hemi_label[hemi_label>0]
        if hemi=='lh':
            hd1=hemi_data.copy()
            hl1=hemi_label.copy()
    combined_data = np.vstack((hd1,hemi_data))
    if normalise:
        combined_data=normalise_subject_data(combined_data)
    combined_label = np.hstack((hl1,hemi_label))
    return combined_data, combined_label.astype(int)

def normalise_subject_data(subject_data):
    """normalise subject feature data to fix between 0 and 1"""
#    subject_data = np.divide(subject_data-np.min(subject_data,axis=0),np.max(subject_data,axis=0)-np.min(subject_data,axis=0),out=np.zeros_like(subject_data), where=subject_data!=0)
    subject_data = np.divide(subject_data-np.mean(subject_data,axis=0),np.std(subject_data,axis=0),out=np.zeros_like(subject_data), where=subject_data!=0)
    return subject_data

def get_histology(listids):
    """return histological classifications for list of subjects"""
    histology={'FCD_2A':[],
              'FCD_2B':[],
              'FCD_3':[],
               'FCD_1':[],
              'all_patients':[]}
    histology['all_patients']=listids
    no_histology=[]
    for ids in listids:
        if ids=='MELD_H10_3T_FCD_0008':
            #missing subject to be followed up
            pass
        else:
            histo_id=get_demographic_feature(ids,'Histo')
            try:
                histo_id=histo_id[0][0]
            except IndexError:
                spl=ids.split('_')
                spl[4]='0'+spl[4]
                n_ids='_'.join(spl)
                histo_id=get_demographic_feature(n_ids,'Histo')
                histo_id=histo_id[0][0]
            if '2A' in histo_id or 'IIA' in histo_id:
                histology['FCD_2A'].append(ids)
            elif '2B' in histo_id or 'IIB' in histo_id:
                histology['FCD_2B'].append(ids)
            elif '3A' in histo_id or 'IIIA' in histo_id:
                histology['FCD_3'].append(ids)
            elif '3B' in histo_id or 'IIIB' in histo_id:
                histology['FCD_3'].append(ids)
            elif '3C' in histo_id or 'IIIC' in histo_id:
                histology['FCD_3'].append(ids)
            elif '3D' in histo_id or 'IIID' in histo_id:
                histology['FCD_3'].append(ids)
            elif '1A' in histo_id or 'IA' in histo_id:
                histology['FCD_1'].append(ids)
            elif '1B' in histo_id or 'IB' in histo_id:
                histology['FCD_1'].append(ids)
            else:
                no_histology.append(ids)
    return histology, no_histology

def histology_per_subject(listids):
    """return histological classification for subjects, in order of input subjects"""
    histology=[]
    for ids in listids:

        if ids=='MELD_H10_3T_FCD_0008':
            #missing subject to be followed up
            histology.append('NaN')
        else:
            histo_id=get_demographic_feature(ids,'Histo')
            try:
                histo_id=histo_id[0][0]
            except IndexError:
                spl=ids.split('_')
                spl[4]='0'+spl[4]
                n_ids='_'.join(spl)
                histo_id=get_demographic_feature(n_ids,'Histo')
                histo_id=histo_id[0][0]
            if '2A' in histo_id or 'IIA' in histo_id:
                histology.append('FCD_2A')
            elif '2B' in histo_id or 'IIB' in histo_id:
                histology.append('FCD_2B')
            elif '3A' in histo_id or 'IIIA' in histo_id:
                histology.append('FCD_3')
            elif '3B' in histo_id or 'IIIB' in histo_id:
                histology.append('FCD_3')
            elif '3C' in histo_id or 'IIIC' in histo_id:
                histology.append('FCD_3')
            elif '3D' in histo_id or 'IIID' in histo_id:
                histology.append('FCD_3')
            elif '1A' in histo_id or 'IA' in histo_id:
                histology.append('FCD_1')
            elif '1B' in histo_id or 'IB' in histo_id:
                histology.append('FCD_1')
            else:
                histology.append('NaN')
    return histology


def lesion_areas(list_ids):
    """calculate lesion areas for all patients, will also return hemisphere and lobe"""
    """ lesion areas are proportion of the hemisphere that is lesion """
    data_path=paths.BASE_PATH
    area = nb.freesurfer.read_morph_data(os.path.join(data_path,'fsaverage_sym/surf/lh.area'))
    cortex = nb.freesurfer.read_label(os.path.join(data_path,'fsaverage_sym/label/lh.cortex.label'))
    lobes_i,_,lobes_labels = nb.freesurfer.read_annot(os.path.join(data_path,'fsaverage_sym/label/lh.lobes.annot'))
    total_area=np.sum(area[cortex])
    areas=[]
    hemis=[]
    lobes=[]
    
    for subject in list_ids:
        hemi=get_les_hemi(subject)
        if hemi is not None:
            lesion = get_feature_values(subject,hemi,'.on_lh.lesion.mgh').astype(bool)
            lesion_area=np.sum(area[lesion])/total_area
            hemis.append(hemi)
            areas.append(lesion_area)
            locations=np.unique(lobes_i[lesion],return_counts=True)
            lobe=lobes_labels[locations[0][np.argmax(locations[1])]].decode()
            #set cingulate and insula to second most prevelant
            if lobe in ['cingulate','insula']:
                try:
                    lobe=lobes_labels[locations[0][np.argsort(locations[1])[-2]]].decode()
                except IndexError:
                    #one fail case was cingulate next to frontal, so frontal
                    lobe='frontal'
            #lobe=lobe.decode('UTF-8')
            lobes.append(lobe)
        #lobes.append(lobes_labels[locations[0][np.argmax(locations[1])]])
        else:
            areas.append(float('NaN'))
            hemis.append(str('NaN'))
            lobes.append(str('NaN'))
    return areas, hemis, lobes


def glasser_areas(list_ids):
    """calculate glasser roi for all patients
    if no lesion mask, roi is 0"""
    
    data_path=paths.BASE_PATH
    glasser_i,_,glasser_labels = nb.freesurfer.read_annot(os.path.join(data_path,'fsaverage_sym/label/lh.HCP-MMP1.annot'))
    glasser_list=[]
    for subject in list_ids:
        hemi=get_les_hemi(subject)
        glasser_vector = np.zeros(len(glasser_labels))
        if hemi is not None:
            lesion = get_feature_values(subject,hemi,'.on_lh.lesion.mgh').astype(bool)
            locations=np.unique(glasser_i[lesion])
            glasser_vector[locations]=1
        glasser_list.append(glasser_vector)
        
    return glasser_list, glasser_labels

