Quick Tutorial
==============

Instantiating an Ontology object
---------------------------------

Ontology objects can be created in several ways

1. Through the __init__ constructor::
     
     ont = Ontology(hierarchy, mapping)

2. Loading from NDEx::

     ont = Ontology.from_ndex('8bfa8318-55ed-11e7-a2e2-0660b7976219')

3. Loading from text-based tables::

     ont = Ontology.from_table('example.txt')

4. Assembly from a similarity network using the CLIXO algorithm::
    
     # load table
     similarity = blah
     ont = Ontology.run_clixo(ont)

Inspecting and manipulating the contents of an ontology
-------------------------------------------------------

An Ontology stores four fundamental types of information

1. genes : List of genes 
2. terms : List of term
3. gene_2_term or term_2_gene : dictionary mapping genes to terms, or vice versa
4. child_2_parent or parent_2_child  : dictionary mapping child to parent terms, or vice versa

Alternatively, the hierarchical connections can be viewed as a matrix::

  Ontology.get_connectivity_matrix

A summary of an Ontology's object, i.e. the number of genes, terms, and connections, can be printed by `Ontology.summary()`

Other functions for manipulating the contents of the Ontology

1. rename
2. collapse_node
3. delete_genes
4. delete_terms

Preprocessing ontologies
------------------------

DDOT provides several convenience functions for processing Ontologies into a desirable structure

1. Ontology.collapse

2. Ontology.mutual_collapse

3. Ontology.propagate_annotations

4. Ontology.parse_obo, Ontology.parse_gaf

Analytical functions
---------------------

1. Ontology.align_hierarchy : align two hierarchies

2. Ontology.get_features : generate ontotypes from Yu et al.

Sharing and visualizing ontologies with NDEx and HiView
--------------------------------------------------------

1. to_ndex

http://hiview.ucsd.edu

Interfaces with other libraries
-------------------------------

1. from_pandas
2. from_pandas
3. to_igraph
4. to_networkx
5. to_NdexGraph

