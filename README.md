# The Data-Driven Ontology Toolkit (DDOT)

The Data-Driven Ontology Toolkit (DDOT) facilitates the inference, analysis, and visualization of biological hierarchies using a data structure called an ontology.

* Open source under MIT license.
* Supports Python 2.7 or >=3.6.

# Documentation

For a quick start, please see the [tutorial](examples/Tutorial.ipynb) and other Jupyter notebooks in the [examples](examples) folder.

For further documentation, please see http://ddot.readthedocs.io/. This includes a description of the [Ontology class](http://ddot.readthedocs.io/en/latest/ontology.html) and a list of [utility functions](http://ddot.readthedocs.io/en/latest/utils.html).

Please post questions or issues to the [Google Groups forum](https://groups.google.com/forum/#!forum/ontology).

# Installation

DDOT requires the following software

* Python v2.7 or >=3.6
* [numpy](https://docs.scipy.org/doc/)
* [scipy](https://docs.scipy.org/doc/)
* [pandas>=0.20](http://pandas.pydata.org/)
* [networkx=1.11](https://networkx.github.io/). Note that networkx>=2.0 is incompatible with ddot.
* [python-igraph](http://igraph.org/python/). Recommend installing through [anaconda](https://anaconda.org/conda-forge/python-igraph) or [pip](https://pypi.python.org/pypi/python-igraph/0.7).
* [ndex-dev](https://github.com/ndexbio/ndex-python). Recommend installing through [pip](https://pypi.python.org/pypi/ndex-dev).
* [tulip-python](https://pypi.python.org/pypi/tulip-python). Recommend installing through [pip](https://pypi.python.org/pypi/tulip-python).

The recommended method for installing these dependencies is to use the [Anaconda](https://conda.io/docs/user-guide/install/download.html) distrubution of Python, and then install Python packages via the conda and pip repositories.

  ```bash
  # Create and activate a virtual environment (optional, but recommended).
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
  ```   

## Install the `ddot` Python package

After dependencies are satisfied, install ddot with pip.

  ```bash
  pip install git+git://github.com/michaelkyu/ddot.git
  ```

If a previous version of ddot is already installed, then force an upgrade with the `--upgrade` option.

  ```bash
  pip install --upgrade git+git://github.com/michaelkyu/ddot.git
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

Yu MK, Ma J, Ono K, Zheng F, Fong S, Gary A, Chen J, Demchak B, Pratt
D, Ideker T. "A swiss-army knife for hierarchical modeling of
biological systems." (in preparation)
