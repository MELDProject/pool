#Functions for the logistic regression analysis. Called by some notebooks and scripts.
import numpy as np
import pool.paths as paths
import os
import nibabel as nb
import pandas as pd
from sklearn.linear_model import LogisticRegression
import h5py
import pool.data_loader as dl

#Calculate the coefficients at a particular vertex
def calculate_log_res_stats(vertex_vector, demographics_matrix):
    """vertex vector - lesion or non-lesion for each patient
    demographics - demographic variables for patient eg age of onset for each patient
    """
    clf = LogisticRegression(solver='liblinear',tol=0.1)
    coefs = clf.fit(demographics_matrix,vertex_vector).coef_
    return coefs

#Calcualte coefficients at a set of given vertices
def calculate_vertex_coefs(demo_features,lesions,vertices,min_number_subjects=1, sf = False):
    """calculate coefficients for all vertices of interest
    input matrix of features
    lesion masks
    vertex indices
    minimum number of lesions needed to calculate 
    sf - seizure freedom. switch predictor and variables so that can predict sf"""
    
    
    coefs = np.zeros((len(vertices),demo_features.shape[1]))
    for k,vertex in enumerate(vertices):
    #set minimum number of subjects
        vertex_vector=lesions[:,vertex]
        if np.sum(vertex_vector)<min_number_subjects:
            pass
        else:
            if sf:
                coefs[k,]=calculate_log_res_stats(demo_features.ravel(), vertex_vector.reshape(-1,1))
            else:
                coefs[k,]=calculate_log_res_stats(vertex_vector,demo_features)
    return coefs


def one_hot_code_features(demo_features):
    for feature in demo_features.columns:
        if demo_features[feature].dtypes =='O':
            if feature=='Hemisphere':
                demo_features['Hemisphere']= (demo_features['Hemisphere']=='lh').astype(int)
            else:
                demo_features = pd.concat([demo_features,pd.get_dummies(demo_features[feature], prefix=feature)],axis=1)
                demo_features.drop([feature],axis=1, inplace=True)
    return demo_features

# Load and organise data for regression
def prepare_data(random_cohort_index = None, features=['Age of onset', 'Sex','Hemisphere','Ever reported MRI negative','Duration'], lesions_file=os.path.join(paths.data_dir,'lesions.npz'),
                 demographics_file = os.path.join(paths.data_dir,'demographics_qc.csv'),
                 one_hot=True):
    #Load lesions
    lesions=dl.load_lesions(lesions_file)
    #load demographics
    demographics = pd.read_csv(demographics_file)
    demographics = demographics[demographics['lesion_masked'].astype(bool)]
    #Normalise lesion area
    if 'Lesion area' in features:
        demographics['Lesion area'] = (demographics['Lesion area'] - np.min(demographics['Lesion area']))/(np.max(demographics['Lesion area']) - np.min(demographics['Lesion area']))
    demo_features = demographics[features]
    vector=demo_features.isnull().any(axis=1)
    #Drop missing features
    demo_features=demo_features.dropna()
    #one hot code features
    if one_hot:
        demo_features = one_hot_code_features(demo_features)
    feature_names = demo_features.columns
    demo_features = demo_features.values
    input_lesions=lesions[~vector]
    #Load cortex labels
    cortical_vertices = nb.freesurfer.io.read_label(os.path.join(paths.BASE_PATH,'fsaverage_sym','label','lh.cortex.label'))
    cortical_vertices = np.sort(cortical_vertices)
    #shuffle cohort 
    if random_cohort_index is not None:
        shuffled_indices = np.random.RandomState(seed=random_cohort_index).permutation(len(input_lesions))
        input_lesions=input_lesions[shuffled_indices]
    return input_lesions, demo_features, cortical_vertices, feature_names



def save_coefs(coefs, h5file,feature_names, random_cohort_index=None, n_perm=1000):
    #create if doesn't exist
    if not os.path.isfile(h5file):
        try:
            f = h5py.File(h5file,'w')
            f.create_dataset('coefs', data=np.zeros_like(coefs))
            f.create_dataset('permuted_coefs', data = np.zeros((n_perm,coefs.shape[0],coefs.shape[1])))
            f.create_dataset('feature_names', data = np.string_(feature_names))
            f.close()
        except OSError:
            pass
    result=None
    while result is None:
        try:
        # connect
            f1 = h5py.File(h5file, 'r+')
            if random_cohort_index is not None:
                f1['permuted_coefs'][random_cohort_index,:,:]=coefs
            else:
                f1['coefs'][:,:] = coefs[:,:]  
            f1.close()
            result = 1
        except OSError:
             pass
    return
