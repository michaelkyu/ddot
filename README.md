# The Data-Driven Ontology Toolkit (DDOT)
A Python library for constructing, analyzing, and visualizing hierarchical models.

# Documentation

http://the-data-driven-ontology-toolkit-ddot.readthedocs.io/en/latest/?

# Installation

Python >=2.7
Python >=3.6

## by anaconda (recommended)

This is the recommended choice for installing ddot because it has a
reliable means to install dependencies and is more lightweight than
docker. Install the [Anaconda
distribution]((https://conda.io/docs/user-guide/install/download.html))
of Python . It is optional but recommended that you create a conda
virtual environment

..

   conda create -n <environment_name>

   # Activate virtual environment
   source activate <environment>

..
   
   # Update conda
   conda update conda
   
   # Install dependencies
   conda install pandas numpy scipy networkx=1.11
   conda install -c conda-forge python-igraph
   pip install tulip-python ndex-dev   
   
   # Install ddot
   conda install -c conda-forge ddot

## by pip

pip install ddot

This will attempt to install the dependencies, but it is recommended that you install dependencies separately to more easily debug any installation errors. 

## from the repository

Download the repository and then run

`python setup.py install --user`

This will attempt to install the dependencies, but it is recommended that you install dependencies separately to more easily debug any installation errors. 

# Dependencies

* [networkx=1.11](https://networkx.github.io/) (You may be using networkx>=2.0, but this is incompatible with ddot)
* [python-igraph](http://igraph.org/python/) (Install through conda or pip)
* [pandas>=0.20](http://pandas.pydata.org/)
* [numpy](https://docs.scipy.org/doc/)
* [scipy](https://docs.scipy.org/doc/)
* [ndex-dev](https://github.com/ndexbio/ndex-python) (Install through pip)
* [tuplip-python](https://pypi.python.org/pypi/tulip-python) (Install through pip)

# Docker image

A docker image of DDOT exists

# Getting started

The <examples> subdirectory contains Jupyter notebooks for example analyses that can be done with DDOT.

# Citing DDOT

If you find DDOT helpful in your research, please cite

Michael Ku Yu, Jianzhu Ma, Keiichiro Ono, Fan Zheng, Samson Fong,
Aaron Gary, Jing Chen, Barry Demchak, Dexter Pratt, Trey Ideker. "A
swiss-army knife for hierarchical modeling of biological systems." (in
preparation)

# Help

TBD: Google Groups Forum
