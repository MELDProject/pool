import numpy as np
import nibabel as nb
import pool.hdf5_io as hio
import os
import pool.paths as paths
import subprocess
import nibabel as nb
from scipy.spatial import distance


def spherical_np(xyz):
    pts=np.zeros(xyz.shape)
    xy=xyz[:,0]**2 + xyz[:,1]**2
    pts[:,0]=np.sqrt(xy+xyz[:,2]**2)
    pts[:,1]=np.arctan2(np.sqrt(xy),xyz[:,2])
    pts[:,2]=np.arctan2(xyz[:,1], xyz[:,0])
    return pts

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
    


def f7(seq):
    #returns uniques but in order to retain neighbour triangle relationship
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))];

    
def get_neighbours_from_tris(tris, label=None):
    """Get surface neighbours from tris
        Input: tris
         Returns Nested list. Each list corresponds 
        to the ordered neighbours for the given vertex"""
    n_vert=np.max(tris+1)
    neighbours=[[] for i in range(n_vert)]
    for tri in tris:
        neighbours[tri[0]].extend([tri[1],tri[2]])
        neighbours[tri[2]].extend([tri[0],tri[1]])
        neighbours[tri[1]].extend([tri[2],tri[0]])
    #Get unique neighbours
    for k in range(len(neighbours)):      
        if label is not None:
            neighbours[k] = set(neighbours[k]).intersection(label)
        else :
            neighbours[k]=f7(neighbours[k])
    return np.array(neighbours);





def compute_islands(area_txtfile,neighbours):
    """calculates islands with same label value"""
    islands=np.zeros_like(area_txtfile).astype(int)
        #array for quick indexing
    indices=np.arange(len(area_txtfile))
    #k is island counter
    k=0
    #while some vertices haven't been assigned
    while np.sum(islands==0)>0:
        k+=1
        #start with lowest unassigned vertex
        cluster=[np.min(np.where(islands==0)[0])]
        #set vertex to island value
        islands[cluster]=k
        #get label value (i.e. inside or outside label, 1/0)
        v_seed_label=area_txtfile[cluster[0]]

        old_cluster=0
        #while size of island k increases
        while np.sum(islands==k)>np.sum(old_cluster):
            #calculate new vertices
            added_vertices=islands==k-old_cluster
            #store for next comparison
            old_cluster=islands==k
            # for the new vertices
            for v in indices[added_vertices]:
                #get the neighbours
                neighbourhood=np.array(neighbours[v])
                #of the neighbours, which are in the same label and not yet part of the island.
                #set these to the island value
                islands[neighbourhood[np.logical_and(area_txtfile[neighbourhood]==v_seed_label,
                                                     islands[neighbourhood]!=k)]]=k
    return islands

##TODO tidy any label file??

def tidy_holes_binary(area_txtfile, neighbours,threshold_area=50, iterations=2):
    """fills small holes in binary surface annotation file.
    threshold_area is maximum area to be switched.
    iterations is the number of times looped through in case of nested holes"""
    #create empty surf file to be filled with island indices
    for i in range(iterations):
        islands=compute_islands(area_txtfile,neighbours)
        #then fill the holes
        new_area=np.copy(area_txtfile).astype(int)
        island_index,counts=np.unique(islands,return_counts=True)
        ordered=np.argsort(counts)
        for ordered_index in ordered:
            if counts[ordered_index]<threshold_area:
                island_i=island_index[ordered_index]
                new_area[islands==island_i]=(area_txtfile[islands==island_i][0]-1)*-1
        area_txtfile=new_area
    return new_area