import numpy as np


def load_lesions(lesions_file):
    """load in lesions file"""
    with np.load(lesions_file,'r') as d:
        lesions=d['arr_0'].astype(int)
    return lesions