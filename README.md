# The Data-Driven Ontology Toolkit (DDOT)
A Python library for constructing, analyzing, and visualizing hierarchical models.

# Documentation

http://the-data-driven-ontology-toolkit-ddot.readthedocs.io/en/latest/?

# Installation

Python >=2.7
Python >=3.6

- by anaconda (recommended)

  This is the recommended method for installing ddot because it has a
  reliable means to install dependencies and is more lightweight than
  docker. Install the [Anaconda
  distribution]((https://conda.io/docs/user-guide/install/download.html))
  of Python . It is optional but recommended that you create a conda
  virtual environment

  ..

     # Create environment (specify a name)
     conda create -n <environment_name>

     # Activate virtual environment
     source activate <environment_name>

  ..
   
     # Update conda
     conda update conda
   
     # Install dependencies
     conda install pandas numpy scipy networkx=1.11
     conda install -c conda-forge python-igraph
     conda install libiconv # Needed for igraph to run properly
     pip install tulip-python ndex-dev
   
     # Install ddot
     conda install -c conda-forge ddot

- by pip

  ..
     pip install ddot

  This will attempt to install the dependencies, but it is recommended that you install dependencies separately to more easily debug any installation errors. 

- from the repository

  Download the repository and then run

  ..
     python setup.py install

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

A docker image of DDOT can be pulled from Docker Hub.

..

    # image with DDOT installed in anaconda3 (Python 3.6)
    docker pull michaelkyu/ddot-anaconda3
    docker run -i -t michaelkyu/ddot-anaconda3
    
    # image with DDOT installed in anaconda2 (Python 2.7)
    docker pull michaelkyu/ddot-anaconda2
    docker run -i -t michaelkyu/ddot-anaconda2

# Getting started

The `examples` folder contains a tutorial (`Tutorial.ipynb`) for
learning basic usage of DDOT. See the other example Jupyter notebooks.

# Citing DDOT

If you find DDOT helpful in your research, please cite

Michael Ku Yu, Jianzhu Ma, Keiichiro Ono, Fan Zheng, Samson Fong,
Aaron Gary, Jing Chen, Barry Demchak, Dexter Pratt, Trey Ideker. "A
swiss-army knife for hierarchical modeling of biological systems." (in
preparation)

# Help

TBD: Google Groups Forum
