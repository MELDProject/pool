from setuptools import setup, find_packages

setup(name='pool',
      version='0.1',
      packages=find_packages(),
      install_requires=['nibabel','h5py','pillow','pandas','scikit-learn','seaborn','ptitprince','statsmodels','scipy',
                        'matplotlib>=3.3.2'],
      package_dir={'pool':'pool'},
     )

