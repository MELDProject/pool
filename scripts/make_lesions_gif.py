## script to make a GIF of lesions
#requires freesurfer
import numpy as np
import pool.paths as paths
import nibabel as nb
import pandas as pd
import pool.preprocessing as pp
import matplotlib.cm as cm
import os
import imageio
import subprocess
def colorise_ribbon(ribbon_section,t1_section,r_threshold,t_threshold):
    t1_gray = cm.gray(np.array(t1_section/t_threshold))
    mask =ribbon_section>0
    coloured_cortex = cm.turbo(np.array(ribbon_section/r_threshold))
    t1_gray[mask]=coloured_cortex[mask]
    return t1_gray*255

lesions=np.load(os.path.join(paths.data_dir,'lesions.npy'))

combined_overlay=np.sum(lesions,axis=0)
demo = nb.load(os.path.join(paths.BASE_PATH,'fsaverage_sym','surf','lh.white.avg.area.mgh'))
pp.save_mgh(os.path.join(paths.data_dir,'lh.lesion_probability_map.mgh'),combined_overlay,demo)

os.environ["SUBJECTS_DIR"]= paths.BASE_PATH
subprocess.call("mri_surf2vol --identity fsaverage_sym --template {} --o {} --hemi lh --surfval {} --fill-projfrac -1 1 0.05".format(
    os.path.join(paths.BASE_PATH,'fsaverage_sym/mri/T1.mgz'),
    os.path.join(paths.data_dir,'lh.lesion_probability_map.mgz'),
    os.path.join(paths.data_dir,'lh.lesion_probability_map.mgh')),
                shell=True)
subprocess.call("mris_apply_reg --src {} --trg {} --streg {} {}".format(
                os.path.join(paths.data_dir,'lh.lesion_probability_map.mgh'),
                os.path.join(paths.data_dir,'rh.lesion_probability_map.mgh'),
                       os.path.join(paths.BASE_PATH,'fsaverage_sym/surf/rh.sphere.left_right'),
                     os.path.join(paths.BASE_PATH,'fsaverage_sym/surf/lh.sphere.reg'),
)
                       ,shell=True)
    

subprocess.call("mri_surf2vol --identity fsaverage_sym --template {} --o {} --hemi rh --surfval {} --fill-projfrac -1 1 0.05".format(
    os.path.join(paths.BASE_PATH,'fsaverage_sym/mri/T1.mgz'),
    os.path.join(paths.data_dir,'rh.lesion_probability_map.mgz'),
    os.path.join(paths.data_dir,'rh.lesion_probability_map.mgh')),
                shell=True)



subprocess.call("mri_convert {} {}".format(os.path.join(paths.data_dir,'rh.lesion_probability_map.mgz'),
                                          os.path.join(paths.data_dir,'rh.lesion_probability_map.nii')),shell=True)
subprocess.call("mri_convert {} {}".format(os.path.join(paths.data_dir,'lh.lesion_probability_map.mgz'),
                                          os.path.join(paths.data_dir,'lh.lesion_probability_map.nii')),shell=True)

subprocess.call("fslmaths {} -add {} {}".format(os.path.join(paths.data_dir,'rh.lesion_probability_map.nii'),
                                               os.path.join(paths.data_dir,'lh.lesion_probability_map.nii'),
                                               os.path.join(paths.data_dir,'lesion_probability_map.nii')),shell=True)

t1=np.array(nb.load(os.path.join(paths.BASE_PATH,'fsaverage_sym/mri/T1.mgz')).dataobj)
p_map=np.array(nb.load(os.path.join(paths.data_dir,'lesion_probability_map.nii.gz')).dataobj)

stack=[]
r_threshold=np.percentile(combined_overlay,95)
t_threshold=np.max(t1)
for section in np.arange(np.int(256)):
    #saggital
    sag=colorise_ribbon(p_map[section],t1[section],r_threshold,t_threshold
                   )
    #horizontal
    hoz=colorise_ribbon(p_map[:,section],t1[:,section],r_threshold,t_threshold
                   )
    #coronal
    coronal=colorise_ribbon(p_map[:,:,section],t1[:,:,section],r_threshold,t_threshold
                   )
    #flip coronal
    coronal = np.flipud(np.rot90(coronal))
    
    combi = np.hstack((sag,coronal,hoz))
    stack.append(combi.astype(np.uint8))
    
imageio.imwrite(os.path.join(paths.fig_dir,'lesion_locations_volume.png'),stack[100])
imageio.mimwrite(os.path.join(paths.fig_dir,'lesion_locations.gif'), stack, fps=20)
subprocess.call('ffmpeg -r 25 -i {} -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" {}'.format(
  os.path.join(paths.fig_dir,'lesion_locations.gif'),
  os.path.join(paths.fig_dir,'lesion_locations.mp4')),shell=True)



