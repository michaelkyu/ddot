# The Data-Driven Ontology Toolkit (DDOT)

The Data-Driven Ontology Toolkit (DDOT) facilitates the inference, analysis, and visualization of biological hierarchies using a data structure called an ontology. 

* Open-source Python package under MIT license. Supports Python 2.7 or >=3.6.
* The [HiView](http://hiview.ucsd.edu) web application visualizes hierarchical structure and the biological evidence for that structure.

# Documentation

For a quick start on DDOT's functionality, please see the [tutorial](examples/Tutorial.ipynb) and other Jupyter notebooks in the [examples](examples) folder.

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
  # Learn more about virtual environments at https://conda.io/docs/user-guide/tasks/manage-environments.html
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

### Install the `ddot` Python package

After dependencies are satisfied, download or clone this repository

  ```bash
  git clone https://github.com/michaelkyu/ddot.git
  ```
  
Next, compile C++ files for running [CliXO v0.3](https://github.com/mhk7/clixo_0.3) and an [ontology alignment algorithm](https://mhk7.github.io/alignOntology/). 

  ```bash
  cd /path/to/ddot_repository
  make
  ```
  
Finally, install ddot using `pip`. If you are installing within a conda virtual environment, remember to enter the environment with `source activate <environment_name>` before running `pip`.

  ```bash
  pip install /path/to/ddot_repository
  ```

# Docker image

A Docker image of DDOT is located online at Docker Hub. To learn more about Docker, see https://docs.docker.com/get-started/

### Download and run image from Docker Hub

For Python 3.6,

```bash
# Download image installed with DDOT in anaconda3 (Python 3.6)
docker pull michaelkyu/ddot-anaconda3
# Run image in a container
docker run -it -p 8888:8888 michaelkyu/ddot-anaconda3
```

For Python 2.7,

```
# Download image installed with DDOT in anaconda2 (Python 2.7)
docker pull michaelkyu/ddot-anaconda2
# Run image in a container
docker run -it -p 8888:8888 michaelkyu/ddot-anaconda2
```

### Using DDOT in Docker

After running the image, you will be inside the container's command line. Here, you can run DDOT in a basic Python terminal

```
(base) root@<container>:/$ python

Python 2.7.14 |Anaconda, Inc.| (default, Dec  7 2017, 17:05:42) 
[GCC 7.2.0] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> import ddot
```

Alternatively, you can run DDOT in [example Jupyter notebooks](examples). To do so, start a Jupyter server in the container's command line

```
(base) root@<container>:/$ jupyter notebook --no-browser --all-root --ip 0.0.0.0 --NotebookApp.token=''
```

Next, open up your web browser and access the notebooks at http://0.0.0.0:8888/notebooks/ddot_examples/. We recommend starting with the tutorial `Tutorial.ipynb`.

# Citing DDOT

If you find DDOT helpful in your research, please cite

Yu MK, Ma J, Ono K, Zheng F, Fong S, Gary A, Chen J, Demchak B, Pratt
D, Ideker T. "A swiss-army knife for hierarchical modeling of
biological systems." (in preparation)
