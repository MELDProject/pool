import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
#subdivide into 2 random groups 1000 time. See if they match.
import ptitprince as pt
import pandas as pd
import numpy as np
from matplotlib import cm
import os
import pool.paths as paths
import nibabel as nb

def permute(n):
    return np.random.permutation(n)

def rank_vector(vector):          
    temp = vector.argsort()
    ranks = np.empty_like(temp)
    ranks[temp] = np.arange(len(vector))
    return ranks
verts,faces=nb.freesurfer.io.read_geometry(os.path.join(paths.BASE_PATH,
                                                        'fsaverage_sym','surf','lh.partial_inflated'))
cortex=nb.freesurfer.read_label(os.path.join(paths.BASE_PATH,
                                                 'fsaverage_sym','label','lh.cortex.label'))
cortex_mask=np.zeros(len(verts)).astype(bool)
cortex_mask[cortex]=1
lesions=np.load(os.path.join(paths.data_dir,'lesions_smoothed.npy'))
n_iter=1000
#sample sizes
cohort_sizes=[20,40,60,80,100,150,200,250]

r_n_values = [] 


for n in np.arange(n_iter):
    bool_i=np.zeros(len(lesions),dtype=bool)
    #randomly shuffle all subjects
    indices=permute(len(lesions))
    #larger cohort of fixed size drawn
    cohort2=rank_vector(np.sum(lesions[indices[250:]],axis=0)[cortex_mask])
    for k,n_subs in enumerate(cohort_sizes):
        bool_i[indices[:n_subs]]=1
        #first smaller cohort drawn
        v=np.sum(lesions[bool_i],axis=0)[cortex_mask]
        cohort1=rank_vector(v)   
        r_n_values.append([np.corrcoef(cohort1,cohort2)[0,1],n_subs]) 

def func(x, a,b, c):
    return (1-a)-b*x**c #+b-(b-0)*np.exp(-d*x)
                      
                      
from scipy.optimize import curve_fit
r_n_values=pd.DataFrame(r_n_values,columns=['R','n_subs'])
mu=np.array(r_n_values.groupby('n_subs').mean()['R'])
std=np.array(r_n_values.groupby('n_subs').std()['R'])
print("MEDIANS")
print(np.array(r_n_values.groupby('n_subs').median()['R']))

lb,up=stats.t.interval(alpha=0.95, df=n_iter-1, loc=mu, scale=
             std) 

colors=cm.Set2(np.arange(len(cohort_sizes)))
my_pal={}
k=-1
for c in np.arange(251):
    if c in cohort_sizes:
        k+=1
        my_pal[c]=colors[k]
    else:
        my_pal[c]=colors[k]
        
import matplotlib
font = {'family' : 'normal',
        'size'   : 20}

matplotlib.rc('font', **font)
fig, ax = plt.subplots(figsize=(8,5))

popt, pcov = curve_fit(func, r_n_values['n_subs'],r_n_values['R'],maxfev=10000)
x_vals=np.arange(400)+10
sigma_ab = np.sqrt(np.diagonal(pcov))
ax.plot(x_vals, func(x_vals, *popt), 'r-',)

popt, pcov = curve_fit(func, cohort_sizes,lb,maxfev=10000)
l_b=func(x_vals, *(popt))
popt, pcov = curve_fit(func, cohort_sizes,up,maxfev=10000)
u_b=func(x_vals, *(popt))
ax.fill_between(x_vals, l_b, u_b,
                 color = 'black', alpha = 0.05)
#plt.scatter(r_n_values['n_subs'],r_n_values['R'])
pt.RainCloud(y='R',x='n_subs',data=r_n_values,ax=ax,order=np.arange(251),width_viol = 20, width_box = 10,
            jitter=5,palette=my_pal)
ax.set_xlim([0,411])
ax.set_xticks(cohort_sizes+[400])
ax.set_xticklabels(cohort_sizes+[400],rotation=45)

ax.set_xlabel('Number of subjects')
ax.set_ylabel('r$_{rank}$')
fig.savefig(os.path.join(paths.fig_dir,'pattern_consitency_cohorts.pdf'))
