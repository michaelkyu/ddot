Tutorial
========

An ontology is a hierarchical arrangement of two types of nodes: (1)
genes at the leaves of the hierarchy and (2) terms at intermediate
levels of the hierarchy. The hierarchy can be thought of as directed
acyclic graph (DAG), in which each node can have multiple children or
multiple parent nodes. In contrast, nodes in a tree (a.k.a. dendogram)
have at most one parent.

The DDOT Python library provides many functions for assembling,
analyzing, and visualizing ontologies.  The main functionalities are
implemented by an "Ontology" class. DDOT can handle ontologies that
are assembled in a data-driven manner as well as manually curated
ontologies like the Gene Ontology.

Creating an Ontology object
---------------------------

Instantiations of the Ontology class can be created in several ways

1. Through the __init__ constructor::
     
     ont = Ontology(hierarchy, mapping, ...)

   The two main parameters are <hierarchy>, a list of (child term,
   parent term) pairs, and <mapping>, a list of (gene, term) pairs.

   **Example:**::

     hierarchy = [('A', 'B'),
                  ('A', 'C'),
		  ('B', 'D'),
		  ('B', 'E'),
		  ('C', 'D'),
		  ('C', 'E'),
		  ('D', 'F'),
		  ('D', 'G'),
		  ('D', 'H'),
		  ('E', 'G'),
		  ('E', 'H'),
		  ]
     mapping = [('F', 'G1'),
                ('G', 'G1'),
		('H', 'G1')]
     ont = Ontology(hierarchy, mapping, parent_child=True)

1. Assembly from a similarity network using the CLIXO algorithm::
    
     ont = Ontology.run_clixo(similarity, ...)

2. Loading from a tab-delimited table or pandas DataFrame::

     ont = Ontology.from_table('example.txt')

4. Loading from the Network DataBase Exchange (NDEx)::

     ont = Ontology.from_ndex(uuid)

   where `uuid` is a specific the unique identifier (UUID) of a
   network on an NDEX server.

     
     
Inspecting the Hierarchical Structure of an Ontology
-------------------------------------------------------

An Ontology stores four fundamental types of information

1. genes : List of genes 
2. terms : List of term
3. gene_2_term or term_2_gene : dictionary mapping genes to terms, or vice versa
4. child_2_parent or parent_2_child  : dictionary mapping child to parent terms, or vice versa

Alternatively, the hierarchical connections can be viewed as a matrix::

  Ontology.connected()

A summary of an Ontology's object, i.e. the number of genes, terms, and connections, can be printed::

  print ont

Direct manipulation of the Ontology
-----------------------------------

DDOT provides several convenience functions for processing Ontologies into a desirable structure

1. Renaming genes and terms.

   .. class:: ddot.Ontology.rename
	      
2. Removing genes and terms.

   .. class:: ddot.Ontology.delete

Currently, there are no functions for adding genes and terms. If this
is needed, then we recommend creating a new Ontology or manipulating
the contents in a different library, such as NetworkX or igraph, and
transforming the results into Ontology.

Collapsing the Ontology
------------------------
1. Ontology.collapse

2. Ontology.mutual_collapse


Expansion of reduction of transitive connnections in the hierarchy
------------------------------------------------------------------

Ontology.propagate_annotations


Alignment of Ontologies
-----------------------

Ontology.align_hierarchy : align two hierarchies


Creating Ontotypes
------------------

Ontology.get_features : generate ontotypes from Yu et al.


Exporting an Ontology
---------------------

1. Ontology.to_ndex : http://hiview.ucsd.edu
2. from_pandas, from_table
3. to_pandas, to_table
   
Visualizing an Ontology
-----------------------

Because an Ontology is a general hierarchical structure, it is
difficult to visualize it in a 2D-layout without having edge
crossings. As one solution, DDOT allows you to show only a subset of
the edges that form a spanning tree of the DAG. In turn, this tree can
be visualized much more easily using several existing algorithms.


1. Ontology.get_tree()

2. Ontology.unfold()


CX file documentation: <link>


Interfaces with other libraries
-------------------------------

1. to_igraph

2. to_networkx

3. to_NdexGraph

