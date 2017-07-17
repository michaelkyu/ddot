# The Data-Driven Ontology Toolkit (DDOT)

The Data-Driven Ontology Toolkit (DDOT) implements four recurring analytical procedures:

1. Expand Gene Set: Identify genes whose functions are highly similar to a seed set of genes.
2. Build Data-driven Hierarchy: Run the CLiXO algorithm to derive a hierarchy relating 
3. Hierarchical Alignment: Identify terms in a hierarchy that highly match terms in a reference hierarchy, such as the Gene Ontology.
4. Hierarchical Viewer: Interactively browse the structure of a hierarchy, including the data that supports the existence of each term in the hierarchy.

Note that the focus of DDOT is not on semantic tasks (the typical focus of other tools for analyzing biomedical ontologies). Instead, DDOT is focused on the hierarchical structures that reflect the physical reality of a cell, disease, or greater biological system. DDOT focuses on creating, manipulating, and analyzing the hierarchical structure of an ontology for downstream tasks in machine learning.

# Dependencies

## Python Packages
* Python (major version 2)
* NumPy
* SciPy
* Pandas
* python-igraph

# Installation

TODO: installation through PyPI.

`pip install ddot`

# DDOT as Python package
DDOT can be invoked in two ways. The first way is as a Python package. DDOT takes advantage of an expansive data analysis ecosystem, including numerical and statistical libraries (NumPy, SciPy, and pandas), machine learning libraries (scikit-learn), general network analysis (NetworkX, igraph), plotting libraries (matplotlib), and iPython notebooks. The core of this package is an “Ontology” class that mediates analyses with a minimal number of lines of code in an object-oriented fashion. For example, the function Ontology.run_CLiXO() will execute the CLiXO algorithm and return an Ontology object. Further analysis of this ontology can be easily done using other methods, such Ontology.align_hierarchies() for comparison with a reference hierarchy. These two functions provide convenient Python wrappers of C++ implementations that we have developed previously. Some of our code has existed as C++ libraries and others as iPython notebooks, requiring a lot of extra boilerplate code to interface between C++ and Python.

# DDOT as a REST service

The second manner that DDOT can be used is a set of four RESTful services, each of which exactly implements one of the above procedures.

# Jupyter Notebook Examples

The <examples> subdirectory contains Jupyter notebooks for example analyses that can be done with DDOT.
