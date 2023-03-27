conda create -n md-automation python=3.10 ambertools=22.3 compilers -c conda-forge -y
conda activate md-automation
conda install -c conda-forge openmm=8.0 cudatoolkit=11 -y
# using the newest version of pdbfixer, since the bugfix & improvement are not yet released
pip install git+https://github.com/openmm/pdbfixer
# small molecule parameterization, protein pH adjustment, and ligand preparation
conda install -c conda-forge acpype openbabel rdkit pdb2pqr -y
# small molecule download
pip install beautifulsoup4
# for debugging
conda install ipython jupyterlab -y
pip install 'ipywidgets>=7.6.0,<8'
# for beautify logger
conda install -c conda-forge ascii-art -y
