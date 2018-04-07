# The Data-Driven Ontology Toolkit (DDOT)
A Python library for constructing, analyzing, and visualizing hierarchical models.

Supports Python v2.7 and >=3.6

# Documentation

http://the-data-driven-ontology-toolkit-ddot.readthedocs.io/en/latest/?

# Installation

There are multiple ways to install DDOT.

* Install by anaconda (recommended)

  This is the easiest method to install dependencies. First, install the [Anaconda
  distribution](https://conda.io/docs/user-guide/install/download.html)
  of Python.

  ```bash
  # Create anad activate a virtual environment (optional but recommended)
  conda create -n <environment_name>
  source activate <environment_name>

  # Update conda
  conda update conda
   
  # Install dependencies
  conda install pandas numpy scipy networkx=1.11
  conda install -c conda-forge python-igraph
  conda install libiconv # Needed for igraph to run properly
  pip install tulip-python ndex-dev
   
  # Install ddot
  conda install -c conda-forge ddot
  ```

* Install by pip

  ```bash
  pip install ddot
  ```

  This will attempt to install the dependencies, but it is recommended that you install dependencies separately to more easily debug any installation errors. 

* Install by downloading the repository and then running

  ```bash
  python setup.py install
  ```

  This will attempt to install the dependencies, but it is recommended that you install dependencies separately to more easily debug any installation errors. 

# Dependencies

* Python v2.7 or >=3.6
* [numpy](https://docs.scipy.org/doc/)
* [scipy](https://docs.scipy.org/doc/)
* [pandas>=0.20](http://pandas.pydata.org/)
* [networkx=1.11](https://networkx.github.io/) You might have networkx>=2.0 already installed, but this is incompatible with ddot.
* [python-igraph](http://igraph.org/python/) Recommend installing through [conda](https://anaconda.org/conda-forge/python-igraph) or [pip](https://pypi.python.org/pypi/python-igraph/0.7).
* [ndex-dev](https://github.com/ndexbio/ndex-python) Recommend installing through [pip](https://pypi.python.org/pypi/ndex-dev).
* [tuplip-python](https://pypi.python.org/pypi/tulip-python) Recommend installing through [pip](https://pypi.python.org/pypi/tulip-python).

# Docker image

A docker image of DDOT can be pulled from Docker Hub.

```bash
# image with DDOT installed in anaconda3 (Python 3.6)
docker pull michaelkyu/ddot-anaconda3
docker run -i -t michaelkyu/ddot-anaconda3
   
# image with DDOT installed in anaconda2 (Python 2.7)
docker pull michaelkyu/ddot-anaconda2
docker run -i -t michaelkyu/ddot-anaconda2
```

# Getting started

The [examples](examples) folder contains a tutorial (`Tutorial.ipynb`) for
learning basic usage of DDOT. See the other example Jupyter notebooks.

# Citing DDOT

If you find DDOT helpful in your research, please cite

Michael Ku Yu, Jianzhu Ma, Keiichiro Ono, Fan Zheng, Samson Fong,
Aaron Gary, Jing Chen, Barry Demchak, Dexter Pratt, Trey Ideker. "A
swiss-army knife for hierarchical modeling of biological systems." (in
preparation)

# Help

TBD: Google Groups Forum
