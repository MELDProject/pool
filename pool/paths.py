import os


for BASE_PATH in ['/Users/sophie/Dropbox/meld/Data/',
'/Users/konradwagstyl/Dropbox (Cambridge University)/MELD/Data/',
'/Users/sophieadler/Dropbox/meld/Data/',
'/home/kwagstyl/Dropbox/MELD/Data/',
'/rds/project/kw350/rds-kw350-meld/meld_data/Data/'
]:
    if os.path.exists(BASE_PATH):
        print('Setting BASE_PATH to ' + BASE_PATH)
        break
if not os.path.exists(BASE_PATH):
    print('WARNING: BASE_PATH not found, setting to None, need to add it to pool/paths.py')
    BASE_PATH = None
    

for data_dir in ['/data1/MELD/pool/data',
                 '/Users/sophieadler/Desktop/MELD_github_repos/pool/data',
'/home/kw350/software/pool/data','../data'
]:
    if os.path.exists(data_dir):
        print('Setting data_dir to ' + data_dir)
        break
if not os.path.exists(data_dir):
    print('WARNING: data_dir not found, setting to None, need to add it to pool/paths.py')
    data_dir = None

for fig_dir in ['/data1/MELD/pool/figures',
                 '/Users/sophieadler/Desktop/MELD_github_repos/pool/figures',
'/home/kw350/software/pool/figures']:
    if os.path.exists(fig_dir):
        print('Setting fig_dir to ' + fig_dir)
        break
if not os.path.exists(data_dir):
    print('WARNING: fig_dir not found, setting to None, need to add it to pool/paths.py')
    fig_dir = None
