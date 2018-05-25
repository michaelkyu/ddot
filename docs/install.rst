Installation
============

DDOT supports Python 2.7 or >=3.5. DDOT can be installed through standard Python environments or through a Docker image.

Installation Methods
--------------------

There are multiple ways to install DDOT.

-  Install by anaconda (recommended)

   This is the easiest method to install dependencies. First, install
   the `Anaconda distribution`_ of Python.

   .. code:: bash

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

-  Install by pip

   .. code:: bash

      pip install ddot

   This will attempt to install the dependencies, but it is recommended
   that you install dependencies separately to more easily debug any
   installation errors.

-  Install from source. Download the `source code`_ and then run

   .. code:: bash

      python setup.py install

   This will attempt to install the dependencies, but it is recommended
   that you install dependencies separately to more easily debug any
   installation errors.

Dependencies
------------

-  Python 2.7 or >=3.6
-  `numpy`_
-  `scipy`_
-  `pandas>=0.20`_
-  `networkx=1.11`_ You might have networkx>=2.0 already installed, but this is incompatible with ddot.
-  `python-igraph`_ Recommend installing through `conda`_ or `pip`_.
-  `ndex-dev`_ Recommend installing through `pip <https://pypi.python.org/pypi/ndex-dev>`__.
-  `tuplip-python`_ Recommend installing through `pip <https://pypi.python.org/pypi/tulip-python>`__.

Docker image
------------

A docker image of DDOT can be pulled from Docker Hub.

.. code:: bash

   # image with DDOT installed in anaconda3 (Python 3.6)
   docker pull michaelkyu/ddot-anaconda3
   docker run -i -t michaelkyu/ddot-anaconda3

   # image with DDOT installed in anaconda2 (Python 2.7)
   docker pull michaelkyu/ddot-anaconda2
   docker run -i -t michaelkyu/ddot-anaconda2


.. _source code: https://github.com/michaelkyu/ddot   
.. _Anaconda distribution: https://conda.io/docs/user-guide/install/download.html
.. _numpy: https://docs.scipy.org/doc/
.. _scipy: https://docs.scipy.org/doc/
.. _pandas>=0.20: http://pandas.pydata.org/
.. _networkx=1.11: https://networkx.github.io/
.. _python-igraph: http://igraph.org/python/
.. _conda: https://anaconda.org/conda-forge/python-igraph
.. _pip: https://pypi.python.org/pypi/python-igraph/0.7
.. _ndex-dev: https://github.com/ndexbio/ndex-python
.. _tuplip-python: https://pypi.python.org/pypi/tulip-python
.. _examples: examples
