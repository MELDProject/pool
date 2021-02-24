import numpy as np
import nibabel as nb
import pool.hdf5_io as hio
import os
import pool.paths as paths
import subprocess
import nibabel as nb
from scipy.spatial import distance



def tidy_features(feature_list):
    """function to tidy up values of feature"""
    filter_values=[555,'555',666,'666','NO','No','No f/u']
    for value in filter_values:
        feature_list[feature_list==value]=float('NaN')
    return np.array(feature_list,dtype=np.float)


def load_lesions_and_hemis(listids, smoothing = 10):
    """load lesions based on subject ids
    returns cortex only data binary for where lesions are
    smoothing 0 means no smoothing.
    """
    cortex=nb.freesurfer.read_label(os.path.join(paths.BASE_PATH,
                                                 'fsaverage_sym','label','lh.cortex.label'))
    cortex_bin=np.zeros(163842).astype(bool)
    cortex_bin[cortex]=1
    basic_overlays=[]
    lesional_ids=[]
    hemispheres=[]
    hemis=['lh','rh']
    lesional_ids = []
    for ids in listids:
        hemi=hio.get_les_hemi(ids)
        if hemi:
                try:
                    overlay= np.round(hio.get_feature_values(ids,hemi,'.on_lh.lesion.mgh')[:])
                    #smooth a bit to adjust for registration issues and, conservative masks etc.
                    if smoothing>0:
  		#try without smoothing for smoothed output maps for clustering
                        overlay = np.ceil(smoothing_fs(overlay,fwhm=smoothing)).astype(int)
                         #overlay = smoothing_fs(overlay,fwhm=smoothing)
                    basic_overlays.append(overlay)
                    lesional_ids.append(1)
                    hemispheres.append(hemis.index(hemi))
                except TypeError:
                    print('skipping ',ids)  
                    lesional_ids.append(0)
                    hemispheres.append(np.nan)
        else:
            lesional_ids.append(0)
    basic_overlays=np.array(basic_overlays).astype(int) #[:,cortex_bin].astype(int)             
    return basic_overlays, np.array(lesional_ids), np.array(hemispheres)


def smoothing_fs(overlay,fwhm, subject='fsaverage_sym',hemi='lh'):
    subjects_dir=paths.BASE_PATH
    os.environ["SUBJECTS_DIR"]= paths.BASE_PATH
    
    """smooth surface overlay on fsaverage_sym"""
    tmpdir='/tmp/' + str(np.random.randint(10000))
    os.mkdir(tmpdir)
    dum=nb.load(os.path.join(subjects_dir,subject,'surf',hemi+'.white.avg.area.mgh'))
    save_mgh(os.path.join(tmpdir,hemi+'.tmp.mgh'),overlay,dum)
    subprocess.call("mris_fwhm --s " + subject + " --hemi " + hemi + " --cortex --smooth-only --fwhm " + str(fwhm) + " --i "
                            + os.path.join(tmpdir,hemi+'.tmp.mgh') + " --o " + os.path.join(tmpdir,hemi+'.sm_tmp.mgh'), shell=True)
    overlay_smoothed=load_mgh(os.path.join(tmpdir,hemi+".sm_tmp.mgh"))
    subprocess.call("rm -r " + tmpdir, shell =True)
    return overlay_smoothed


def load_mgh(filename):
    """ import mgh file using nibabel. returns flattened data array"""
    mgh_file=nb.load(filename)
    mmap_data=mgh_file.get_data()
    array_data=np.ndarray.flatten(mmap_data)
    return array_data;


def save_mgh(filename,array, demo):
    """ save mgh file using nibabel and imported demo mgh file"""
    mmap=np.memmap('/tmp/tmp', dtype='float32', mode='w+', shape=demo.get_data().shape)
    mmap[:,0,0]=array[:]
    output=nb.MGHImage(mmap, demo.affine, demo.header)
    nb.save(output, filename)
    
    
def com_lesion_map(list_ids,lesions,lesion_size=5000):
    """calculate lesion map where lesions have fixed, circular size"""
    
    data_path=paths.BASE_PATH

    lesions = lesions.astype(bool)
    lesions_com=np.zeros_like(lesions)
    verts,faces=nb.freesurfer.io.read_geometry(os.path.join(paths.BASE_PATH,
                                                    'fsaverage_sym','surf','lh.sphere'))
    centers=np.zeros((len(list_ids),3))
    for k,subject in enumerate(list_ids):
        centers[k]= np.mean(verts[lesions[k]],axis=0)
    distances=distance.cdist(centers,verts)
    for k,subject in enumerate(list_ids):
        vertices=distances[k].argsort()[:int(lesion_size)]
        lesions_com[k,vertices]=1
        
    return lesions_com