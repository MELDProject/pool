
import os
import pool.stats_plotting as sp
import pool.paths as paths

#sp.plot_dataset(os.path.join(paths.BASE_PATH,'regression_models','full_presurgical.hdf5'),p_sigs=[0.05],plots=['perm_pval'],
#             outdir=os.path.join(paths.fig_dir,'regression_plots','full_presurg'))

sp.plot_dataset(os.path.join(paths.BASE_PATH,'regression_models','full_postsurgical.hdf5'),p_sigs=[0.05],plots=['perm_pval'],
             outdir=os.path.join(paths.fig_dir,'regression_plots','full_postsurg'))

sp.plot_dataset(os.path.join(paths.BASE_PATH,'regression_models','full_presurgical_reduced.hdf5'),p_sigs=[0.05],plots=['perm_pval'],
             outdir=os.path.join(paths.fig_dir,'regression_plots','full_presurg_red'))
#sp.plot_dataset(os.path.join(paths.BASE_PATH,'regression_models','full_presurgical_area.hdf5'),p_sigs=[0.05],plots=['perm_pval'],
 #            outdir=os.path.join(paths.fig_dir,'regression_plots','full_presurg_area'))
