import numpy as np
import nibabel as nb
import scipy.stats as st
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
import pool.paths as paths

def calculate_spun_stats(map1,map2,test,spins):
    """calculate statistic between pairs of maps with a function specified in test
    spin map2 and recalculate
    test must return stat of interest only"""

    base_stat = test(map1,map2)
    n_spins=len(spins)
    spun_stats = np.zeros(n_spins)
    for spin in np.arange(n_spins):
        spun_m2 = map2[spins[spin]]
        spun_stats[spin] = test(map1,spun_m2)
        
    return base_stat, spun_stats

def spearman_r(map1,map2):
    return st.spearmanr(map1,map2)[0]

def pearson_r(map1,map2):
    return np.corrcoef(map1,map2)[0,1]

def chi2(map1,map2):
    c=np.bincount(2 * (map1 ) + (map2 ), minlength=4).reshape(2, 2)
    #odds ratio
    #stat = c[0, 0] * c[1, 1] / (c[1, 0] * c[0, 1])
    try:
        stat=st.chi2_contingency(c)[0]
    except ValueError:
        stat=1

    return stat

def ttest(map1,map2):
    in_ = map1[map2]
    out_ = map1[~map2]
    stat = st.ttest_ind(in_,out_)[0]
    return stat

def neurosynth_binary_annotation(map1,spins,test='chi2'):
    ns_dir = '/data1/bigbrain/phate_testing/neurosynth/fsaverage_sym/'
    cortex_l=nb.freesurfer.read_label(os.path.join(paths.BASE_PATH,
                                         'fsaverage_sym','label','lh.cortex.label'))
    cortex_bin=np.zeros(163842).astype(bool)
    cortex_bin[cortex_l]=1
    cortex=cortex_bin
    topics = os.listdir(ns_dir)
    perm_tstat= np.zeros((len(topics),len(spins)))
    t_stat = np.zeros(len(topics))
    topic_terms=np.zeros(len(topics),dtype=object)
    for k,topic in enumerate(topics):
        topic_terms[k]='_'.join(topic[:-32].split('_')[2:])
        z_map=np.loadtxt(os.path.join(ns_dir,topic),skiprows=1)[:len(cortex)]
        binary_neuro= (z_map>0)[cortex]
        
        t_stat[k],perm_tstat[k]=calculate_spun_stats(map1,binary_neuro,test=eval(test),
spins=spins)
    pvals= np.mean(np.greater_equal(t_stat,perm_tstat.T),axis=0)
    pvals = np.abs(0.5-np.abs(0.5-pvals))
    return t_stat, pvals, topic_terms, perm_tstat


def plot_neurosynth(t_stat,pvals,topic_terms,title='',test='chi2',ax=None):
    if not ax:
        fig = plt.figure(figsize=(4,6))
        ax = fig.add_subplot(1,1,1)
    colors=np.ones(len(topic_terms))*0.7
    colors[pvals<0.05/len(topic_terms)]=1
    cols=cm.bwr(colors)
    
    ordered_pvals=pvals[np.argsort(t_stat)]
    ax.barh(np.arange(len(t_stat)),np.sort(t_stat),color=cols[np.argsort(t_stat)])
    topic_terms[np.argsort(pvals)],t_stat[np.argsort(pvals)],np.sort(pvals)
    ax.set_yticks(np.arange(len(topic_terms)))
    ax.set_yticklabels(topic_terms[np.argsort(t_stat)])
    #plt.yticks(np.arange(len(topics))[ordered_pvals<0.05],topic_terms[np.argsort(t_stat)][ordered_pvals<0.05])
    ax.set_title(title)
    ax.set_xlabel(test)
    return ax
