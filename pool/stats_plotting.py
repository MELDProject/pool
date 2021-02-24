import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pool.paths as paths
import nibabel as nb
#from lesion_characteristics.demographic_location_regression import load_data
import seaborn as sns
from PIL import Image
import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import pool.mesh_tools as mt
import matplotlib
import pool.matplotlib_surface_plotting as msp


def plot_coefficients(coefficients_overlay,filename='tmp.png',plot='coefs'):
    """plot medial and lateral using nilearn."""
    if plot=='coefs':
        vmax=np.mean(coefficients_overlay)+1.2*np.std(coefficients_overlay)
#                     np.abs(np.mean(coefficients_overlay)-1*np.std(coefficients_overlay))])
        vmin=np.mean(coefficients_overlay)-1.2*np.std(coefficients_overlay)
        vmin = np.percentile(coefficients_overlay[coefficients_overlay>np.min(coefficients_overlay)],10)
        vmax = np.percentile(coefficients_overlay[coefficients_overlay>np.min(coefficients_overlay)],90) 
        cmap='viridis'
    elif plot=='perm_pval':
        vmax=0.2
        vmin=0
#        coefficients_overlay = coefficients_overlay<0.05
        cmap='magma'
    surf = nb.freesurfer.io.read_geometry(os.path.join(paths.BASE_PATH,'fsaverage_sym','surf','lh.partial_inflated'))
    msp.plot_surf(surf[0],surf[1], coefficients_overlay, rotate=[90,270],vmax=vmax,vmin=vmin,cmap=cmap,filename=filename)
    
    return


def interpolate_values(values):
 #   """interpolate out to full mesh resolution"""
    import stripy
    sphere = nb.freesurfer.io.read_geometry(os.path.join(paths.BASE_PATH,'fsaverage_sym','surf','lh.sphere'))
    spherical=mt.spherical_np(sphere[0])
    lats,lons = spherical[:,1]-np.pi/2,spherical[:,2]
    mesh=stripy.sTriangulation(lons[:len(values)],lats[:len(values)])
    interpolated_values=mesh.interpolate_nearest(lons,lats,values)[0]
    return interpolated_values

def plot_pvals(permed_pvals,orig_pval,filename, p_sig=0.05, n_factors=1):
    """plot distributions of pvalues across both hemispheres"""
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 35}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots(tight_layout=True,figsize=(15,10))
    #count how many above significance
    #two tailed so p_sig/2
    sig_perm = np.mean(permed_pvals<p_sig/2, axis=1)*100
   # ax.hist(sig_perm,100)
    sig_orig = np.mean(orig_pval<p_sig/2)*100
    thresh_sig = 100*(1 - 0.05/n_factors)
#    print('thresh',thresh_sig)
    threshold=np.percentile(sig_perm,thresh_sig)
    sns.distplot(sig_perm,ax=ax ,label="%  in null")
    ax.plot([threshold,threshold],[0,1],'gray',lw=5,alpha=0.3, label='{}th percentile'.format(int(thresh_sig)))
    ax.plot([sig_orig,sig_orig],[0,1],'r',label='%  in MELD')
    ax.set_ylabel('Count')
    ax.set_xlabel('Permuted fractions of significant vertices')
    ax.legend(loc=[1.01,0.6])
    fig.savefig(filename)
    plt.close()
    return

def smooth_cluster_pvals(pvals,surf, fwhm=5):
    """smooth and cluster pvals"""
    neighbours=mt.get_neighbours_from_tris(surf[1])
    tidied=mt.tidy_holes_binary(pvals<0.025,neighbours,threshold_area=100)
    smoothed=mt.smoothing_fs(tidied,fwhm=5)
    return smoothed

def plot_combined_coefficients_pvals(coefficients_overlay,pvals,filename='tmp.png',mask=None,fwhm=5,
                                    vmin=None,vmax=None, return_p=False):
    """plot medial and lateral using nilearn."""
    if not vmin:
        vmin = np.percentile(coefficients_overlay[coefficients_overlay>np.min(coefficients_overlay)],10)
    if not vmax:
        vmax = np.percentile(coefficients_overlay[coefficients_overlay>np.min(coefficients_overlay)],90) 
    surf = nb.freesurfer.io.read_geometry(os.path.join(paths.BASE_PATH,'fsaverage_sym','surf','lh.partial_inflated'))
    smoothed = smooth_cluster_pvals(pvals,surf,fwhm=fwhm)
    msp.plot_surf(surf[0],surf[1], 
       coefficients_overlay,pvals=smoothed, rotate=[90,270],vmax=vmax,vmin=vmin,filename=filename,cmap='viridis', mask=mask)
    if return_p:
        return smoothed
    return


def load_p_vals(dataset):
    """read in coefficients and pvalues"""
    with h5py.File(dataset,'r') as f:
        factors=f['feature_names'][:]
        coefs=f['coefs'][:,:]
        permuted_coefs=f['permuted_coefs'][:,:]
        feature_names=f['feature_names'][:]
    mask=np.all(permuted_coefs[:,:,:]==0,axis=(0,2))
    #calculate pvalues
    n_perm=permuted_coefs.shape[0]
    tiled_coefs=np.tile(coefs,(n_perm,1,1))
    pvals=np.mean(tiled_coefs>permuted_coefs,axis=0)
    pvals=0.5-np.abs(pvals-0.5)
    #print(permuted_coefs.shape)
    #loop over factors
    #get p value by ranking and dividing by n_perm to give p value relative to other permutations
    permuted_pvals=0.5-np.abs(np.argsort(permuted_coefs,axis=0)/n_perm-0.5)
    permuted_pvals[:,mask,:]=0.3 
    return coefs, pvals, permuted_pvals, mask, factors

def plot_dataset(dataset, p_sigs=[0.05,0.01], plots=['coefs','perm_pval'], plot_combined=True,
                outdir=os.path.join('..','figures','regression_plots')):
    """read in coefficients and pvalues.
    plot dataset and associated figures"""
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
        
    coefs, pvals, permuted_pvals, mask, factors =load_p_vals(dataset)
    
    cortex=np.sort(nb.freesurfer.read_label(os.path.join(paths.BASE_PATH,
                                             'fsaverage_sym','label','lh.cortex.label')))
    cortex_bin = np.zeros(163842).astype(bool)
    cortex_bin[cortex] = 1
    # mask out all permuted coefs which were set to 0
    
    for k,factor in enumerate(factors):
        factor=factor.decode("utf-8")
        print('plotting {}'.format(factor))
        factor_pvalues=pvals[:,k]
        factor_coefs=coefs[:,k]
        #for plotting, hide medial wall
        factor_coefs[mask]=np.min(coefs)
        factor_pvalues[mask]=0.3
        #expand them out
        factors_expanded=np.zeros(np.max(cortex)+1)
        factor_pvalues_expanded =np.zeros(np.max(cortex)+1)
        mask_expanded =np.zeros(np.max(cortex)+1)
        factors_expanded[cortex]=factor_coefs
        factor_pvalues_expanded[cortex]=factor_pvalues
        mask_expanded[cortex]=mask
        mask_expanded[~cortex_bin] = 1
        for plot in plots:
            if plot=='perm_pval':
                values=factor_pvalues_expanded
            else:
                values = factors_expanded

            plot_coefficients(values,
                              filename=os.path.join(outdir,'{}_{}.png'.format(factor,plot)),
                             plot =plot)
        if plot_combined:
            filename=os.path.join(outdir,'{}_{}.png'.format(factor,'combined_overlays'))
            plot_combined_coefficients_pvals(factors_expanded,factor_pvalues_expanded,filename=filename, mask=mask_expanded)

        if p_sigs is not None:
                for p_sig in p_sigs:
                    #correction for multiple comparisons
                    plot_pvals(permuted_pvals[:,~mask,k],
                               factor_pvalues[~mask],  filename=os.path.join(outdir,
                                                                      '{}_comparison_vs_permuted_{}.pdf'.format(factor,p_sig)),
                               p_sig=p_sig, n_factors=len(factors))
    return



def plot_perm(r_s,r_perms,filename):
    p_s = 2*(0.5-np.abs(np.mean(r_s>r_perms,axis=0)-0.5))
    font = {'family' : 'times new roman',
        'weight' : 'normal',
        'size'   : 40}
    matplotlib.rc('font', **font)
    if r_s>0:
        thr=np.percentile(r_perms,97.5)
    else:
        thr=np.percentile(r_perms,2.5)
    levels = np.array([0.05,0.01,0.001])
    index_p = np.sum(levels>p_s).astype(int)-1
    if index_p>=0:
        p_sig = levels[index_p]
    else:
        p_sig = 'NA'
    plt.figure(figsize=(15,10))
    sns.distplot(r_perms,bins=20, label='Permuted statistics')
    plt.plot([r_s,r_s],[0,1.0],c='r',label='Test statistic,\nr$_{rank}$'+'={:.2f}, p<{}'.format(r_s,p_sig),linewidth=4)
    plt.plot([thr,thr],[0,1.0],c='grey',label='Threshold r',lw=5,alpha=0.3,)
    plt.xlabel('Correlation coefficients')
    plt.ylabel('Relative frequency')
    plt.legend(loc=[1.01,0.6])
    plt.tight_layout()
    plt.savefig(os.path.join(paths.fig_dir,filename))