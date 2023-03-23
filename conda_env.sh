conda create -n md-automation python=3.10 ambertools=22.3 compilers -c conda-forge -y
conda activate md-automation
conda install -c conda-forge openmm=8.0 cudatoolkit=11 -y
# using the newest version of pdbfixer, since the bugfix & improvement are not yet released
pip install git+https://github.com/openmm/pdbfixer
# for debugging
conda install ipython jupyterlab -y


