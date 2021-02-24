#Run regression analysis
#Calculating main coefficients is done leaving "random_cohort_index" empty. 
#Calculating random permutations is done by rerunning and specifying random_cohort_index with an integer
#This is to allow parallel execution on the supercomputer.

import numpy as np
import pool.paths as paths
import os
import nibabel as nb
import pandas as pd
from sklearn.linear_model import LogisticRegression
import h5py
import pool.regression_functions as rf
import argparse



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run logistic regression analyses on MELD with permutation testing')
    parser.add_argument('--random_cohort_index', type=int, default=None)
    parser.add_argument('--n_perm', type=int, default=1000)
    parser.add_argument('--feature_list', default='postsurgical')
    args = parser.parse_args()
    n_perm=args.n_perm
    feature_list = args.feature_list
    hdf5_name=f'full_{feature_list}.hdf5'
    random_cohort_index = args.random_cohort_index
    feature_lists = {'presurgical':['Age of onset', 'Sex','Hemisphere','Ever reported MRI negative','Duration', 'Lesion area'],
                      'presurgical_reduced':['Age of onset','Duration'],
                       'duration':['Duration'],
                     'postsurgical':['Seizure free']}
    feature_set=feature_lists[feature_list]
    if not os.path.isdir(os.path.join(paths.BASE_PATH,'regression_models')):
        os.makedirs(os.path.join(paths.BASE_PATH,'regression_models'))
    filename=os.path.join(paths.BASE_PATH,'regression_models',hdf5_name)
    #import demographics and lesions
    lesions,demo_features,cortical_vertices, feature_names=rf.prepare_data(random_cohort_index=random_cohort_index,
                            features=feature_set)
    #run regression
    coefs=rf.calculate_vertex_coefs(demo_features,lesions,cortical_vertices,min_number_subjects=1, sf = feature_list=='postsurgical')
    #save coefficients
    rf.save_coefs(coefs,filename,feature_names,random_cohort_index=random_cohort_index, n_perm=n_perm)
