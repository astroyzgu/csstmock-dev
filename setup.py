from setuptools import setup, find_packages
import glob
setup(
  name="csstmock",
  version="1.0.0",
  description="Construct mock catalog and apply foreground masks for CSST.",
  author="SJTU group",
  author_email="guyizhou@sjtu.edu.cn",
  python_requires=">=3.7",
  install_requires=[
    "matplotlib",
    "numpy",
    "astropy",
    "healpy",
    "h5py",
    "corrfunc",
    "imageio",
    "pytest"
  ],
  scripts=glob.glob('bin/*'),
  packages=find_packages("src"),
  package_dir={"": "src"},
  package_data={"mypkg": ["*.npy"]},
  include_package_data=True
)