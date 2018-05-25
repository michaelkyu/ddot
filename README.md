# The Data-Driven Ontology Toolkit (DDOT)
A Python library for constructing, analyzing, and visualizing hierarchical models.

# Documentation

To get a quick start on DDOT, please see the tutorial (`examples/DDOT_tutorial.ipynb`) in the [examples](examples) folder. Please see the other example Jupyter notebooks.

Further documentation is provided at http://ddot.readthedocs.io/. It includes a description of the [Ontology class](http://ddot.readthedocs.io/en/latest/ontology.html) and a list of [utility functions](http://ddot.readthedocs.io/en/latest/utils.html)

# Installation

DDOT requires the following

* Python v2.7 or >=3.6
* [numpy](https://docs.scipy.org/doc/)
* [scipy](https://docs.scipy.org/doc/)
* [pandas>=0.20](http://pandas.pydata.org/)
* [networkx=1.11](https://networkx.github.io/) networkx>=2.0 is incompatible with ddot.
* [python-igraph](http://igraph.org/python/) Recommend installing through [anaconda](https://anaconda.org/conda-forge/python-igraph) or [pip](https://pypi.python.org/pypi/python-igraph/0.7).
* [ndex-dev](https://github.com/ndexbio/ndex-python) Recommend installing through [pip](https://pypi.python.org/pypi/ndex-dev).
* [tuplip-python](https://pypi.python.org/pypi/tulip-python) Recommend installing through [pip](https://pypi.python.org/pypi/tulip-python).

The recommended method for installing these dependencies is to use the [Anaconda](https://conda.io/docs/user-guide/install/download.html) distrubution of Python, and then install Python packages via the conda and pip repositories.

  ```bash
  # Create anad activate a virtual environment (optional, but recommended).
  # Read more about virtual environments at https://conda.io/docs/user-guide/tasks/manage-environments.html
  conda create -n <environment_name>
  source activate <environment_name>

  # Update conda
  conda update conda
   
  # Install dependencies
  conda install pandas numpy scipy networkx=1.11
  conda install -c conda-forge python-igraph
  conda install libiconv # Needed for igraph to run properly
  pip install tulip-python
  pip install ndex-dev
   
Once the dependencies are satisfied, you can install ddot as follows

## Install `ddot` Python package with `pip`
  ```bash
  pip install git+git://github.com/michaelkyu/ontology.git@v0.2rc1
  conda install -c conda-forge ddot
  ```

## Docker image

A docker image of DDOT can be pulled from Docker Hub.

```bash
# image with DDOT installed in anaconda3 (Python 3.6)
docker pull michaelkyu/ddot-anaconda3
docker run -i -t michaelkyu/ddot-anaconda3
   
# image with DDOT installed in anaconda2 (Python 2.7)
docker pull michaelkyu/ddot-anaconda2
docker run -i -t michaelkyu/ddot-anaconda2
```


# Citing DDOT

If you find DDOT helpful in your research, please cite

Michael Ku Yu, Jianzhu Ma, Keiichiro Ono, Fan Zheng, Samson Fong,
Aaron Gary, Jing Chen, Barry Demchak, Dexter Pratt, Trey Ideker. "A
swiss-army knife for hierarchical modeling of biological systems." (in
review)

# Help

Please post any questions or issues to the DDOT forum at https://groups.google.com/forum/#!forum/ontology
