Tutorial
============

An ontology is a hierarchical arrangement of two types of nodes: (1)
genes at the leaves of the hierarchy and (2) terms at intermediate
levels of the hierarchy. The hierarchy can be thought of as directed
acyclic graph (DAG), in which each node can have multiple children or
multiple parent nodes. DAGs are a generalization of trees (a.k.a.
dendogram), where each node has at most one parent.

The DDOT Python library provides many functions for assembling,
analyzing, and visualizing ontologies. The main functionalities are
implemented in an object-oriented manner by an "Ontology" class. This
class can handle both ontologies that are data-driven as well as those
that are manually curated like the Gene Ontology.

Creating an Ontology object
---------------------------

An object of the Ontology class can be created in several ways. To
demonstrate this, we will build the following ontology

Through the init constructor
------------------------------------

.. code:: 

    # Connections from child terms to parent terms
    hierarchy - [('S3', 'S1'),
		 ('S4', 'S1'),
		 ('S5', 'S1'),
		 ('S5', 'S2'),
		 ('S6', 'S2'),
		 ('S1', 'S0'),
		 ('S2', 'S0')]

    # Connections from genes to terms
    mapping - [('A', 'S3'),
	       ('B', 'S3'),
	       ('C', 'S3'),
	       ('C', 'S4'),
	       ('D', 'S4'),
	       ('E', 'S5'),
	       ('F', 'S5'),
	       ('G', 'S6'),
	       ('H', 'S6')]

    # Construct ontology
    ont - Ontology(hierarchy, mapping)

To and from a tab-delimited table or Pandas dataframe
-----------------------------------------------------

.. code:: 

    ont.to_table('toy_ontology.txt')




.. raw:: html

    <div>
    <style>
	.dataframe thead tr:only-child th {
	    text-align: right;
	}

	.dataframe thead th {
	    text-align: left;
	}

	.dataframe tbody tr th {
	    vertical-align: top;
	}
    </style>
    <table border-"1" class-"dataframe">
      <thead>
	<tr style-"text-align: right;">
	  <th></th>
	  <th>Parent</th>
	  <th>Child</th>
	  <th>EdgeType</th>
	</tr>
      </thead>
      <tbody>
	<tr>
	  <th>0</th>
	  <td>S2</td>
	  <td>S5</td>
	  <td>Child-Parent</td>
	</tr>
	<tr>
	  <th>1</th>
	  <td>S2</td>
	  <td>S6</td>
	  <td>Child-Parent</td>
	</tr>
	<tr>
	  <th>2</th>
	  <td>S1</td>
	  <td>S3</td>
	  <td>Child-Parent</td>
	</tr>
	<tr>
	  <th>3</th>
	  <td>S1</td>
	  <td>S4</td>
	  <td>Child-Parent</td>
	</tr>
	<tr>
	  <th>4</th>
	  <td>S1</td>
	  <td>S5</td>
	  <td>Child-Parent</td>
	</tr>
	<tr>
	  <th>5</th>
	  <td>S0</td>
	  <td>S1</td>
	  <td>Child-Parent</td>
	</tr>
	<tr>
	  <th>6</th>
	  <td>S0</td>
	  <td>S2</td>
	  <td>Child-Parent</td>
	</tr>
	<tr>
	  <th>7</th>
	  <td>S3</td>
	  <td>A</td>
	  <td>Gene-Term</td>
	</tr>
	<tr>
	  <th>8</th>
	  <td>S3</td>
	  <td>C</td>
	  <td>Gene-Term</td>
	</tr>
	<tr>
	  <th>9</th>
	  <td>S4</td>
	  <td>C</td>
	  <td>Gene-Term</td>
	</tr>
	<tr>
	  <th>10</th>
	  <td>S3</td>
	  <td>B</td>
	  <td>Gene-Term</td>
	</tr>
	<tr>
	  <th>11</th>
	  <td>S5</td>
	  <td>E</td>
	  <td>Gene-Term</td>
	</tr>
	<tr>
	  <th>12</th>
	  <td>S4</td>
	  <td>D</td>
	  <td>Gene-Term</td>
	</tr>
	<tr>
	  <th>13</th>
	  <td>S6</td>
	  <td>G</td>
	  <td>Gene-Term</td>
	</tr>
	<tr>
	  <th>14</th>
	  <td>S5</td>
	  <td>F</td>
	  <td>Gene-Term</td>
	</tr>
	<tr>
	  <th>15</th>
	  <td>S6</td>
	  <td>H</td>
	  <td>Gene-Term</td>
	</tr>
      </tbody>
    </table>
    </div>



.. code:: 

    ont - Ontology.from_table('toy_ontology.txt')

From the Network Data Exchange (NDEx). Requires a free user account at http://ndexbio.org/
------------------------------------------------------------------------------------------

.. code:: 

    # Replace with your own NDEx user account
    ndex_server, ndex_user, ndex_pass - 'http://test.ndexbio.org', 'scratch', 'scratch'
    # ndex_user, ndex_pass - 'ddot_test', 'ddot_test'

    url, _ - ont.to_ndex(ndex_server-ndex_server, ndex_user-ndex_user, ndex_pass-ndex_pass)
    print(url)


.. parsed-literal::

    http://dev2.ndexbio.org/v2/network/79ccfa8b-369a-11e8-929a-0660b7976219


.. code:: 

    ont2 - Ontology.from_ndex("http://dev2.ndexbio.org/v2/network/3a103097-35f1-11e8-84e4-0660b7976219")
    print(ont2)


.. parsed-literal::


    8 genes, 7 terms, 9 gene-term relations, 7 term-term relations
    node_attributes: [u'NodeType', 'name', u'x_pos', u'isRoot', u'Vis:Shape', u'y_pos', u'Label', u'Vis:Border Paint', u'Vis:Size', u'Vis:Fill Color', u'Size']
    edge_attributes: [u'Is_Tree_Edge', u'2', u'Vis:Visible', u'EdgeType']


Inspecting the structure of an ontology
---------------------------------------

An Ontology object contains seven attributes:

-  ``genes`` : List of gene names
-  ``terms`` : List of term names
-  ``gene_2_term`` : dictionary mapping a gene name to a list of terms
   connected to that gene. Terms are represented as their 0-based index
   in ``terms``.
-  ``term_2_gene`` : dictionary mapping a term name to a list or genes
   connected to that term. Genes are represented as their 0-based index
   in ``genes``.
-  ``child_2_parent`` : dictionary mapping a child term to its parent
   terms.
-  ``parent_2_child`` : dictionary mapping a parent term to its children
   terms.
-  ``term_sizes`` : A list of each term's size, i.e. the number of
   unique genes contained within this term and its descendants. The
   order of this list is the same as ``terms``. For every ``i``, it
   holds that ``term_sizes[i] - len(self.term_2_gene[self.terms[i]])``

.. code:: 

    ont.genes




.. parsed-literal::

    ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']



.. code:: 

    ont.terms




.. parsed-literal::

    ['S0', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6']



.. code:: 

    ont.gene_2_term




.. parsed-literal::

    {'A': [3],
     'B': [3],
     'C': [3, 4],
     'D': [4],
     'E': [5],
     'F': [5],
     'G': [6],
     'H': [6]}



.. code:: 

    ont.term_2_gene




.. parsed-literal::

    {'S0': [],
     'S1': [],
     'S2': [],
     'S3': [0, 1, 2],
     'S4': [2, 3],
     'S5': [4, 5],
     'S6': [6, 7]}



.. code:: 

    ont.child_2_parent




.. parsed-literal::

    {'S0': [],
     'S1': ('S0',),
     'S2': ('S0',),
     'S3': ('S1',),
     'S4': ('S1',),
     'S5': ('S2', 'S1'),
     'S6': ('S2',)}



Alternatively, the hierarchical connections can be viewed as a matrix,
using ``connected()``

.. code:: 

    ont.connected()




.. parsed-literal::

    array([[ True, False, False, False, False, False, False, False,  True,  True, False,  True, False, False, False],
	   [False,  True, False, False, False, False, False, False,  True,  True, False,  True, False, False, False],
	   [False, False,  True, False, False, False, False, False,  True,  True, False,  True,  True, False, False],
	   [False, False, False,  True, False, False, False, False,  True,  True, False, False,  True, False, False],
	   [False, False, False, False,  True, False, False, False,  True,  True,  True, False, False,  True, False],
	   [False, False, False, False, False,  True, False, False,  True,  True,  True, False, False,  True, False],
	   [False, False, False, False, False, False,  True, False,  True, False,  True, False, False, False,  True],
	   [False, False, False, False, False, False, False,  True,  True, False,  True, False, False, False,  True],
	   [False, False, False, False, False, False, False, False,  True, False, False, False, False, False, False],
	   [False, False, False, False, False, False, False, False,  True,  True, False, False, False, False, False],
	   [False, False, False, False, False, False, False, False,  True, False,  True, False, False, False, False],
	   [False, False, False, False, False, False, False, False,  True,  True, False,  True, False, False, False],
	   [False, False, False, False, False, False, False, False,  True,  True, False, False,  True, False, False],
	   [False, False, False, False, False, False, False, False,  True,  True,  True, False, False,  True, False],
	   [False, False, False, False, False, False, False, False,  True, False,  True, False, False, False,  True]], dtype-bool)



A summary of an Ontology’s object, i.e. the number of genes, terms, and
connections, can be printed ``print(ont)``

.. code:: 

    print(ont)


.. parsed-literal::

    8 genes, 7 terms, 9 gene-term relations, 7 term-term relations
    node_attributes: []
    edge_attributes: [2]


Manipulating the structure of an ontology
-----------------------------------------

DDOT provides several convenience functions for processing Ontologies
into a desirable structure. Currently, there are no functions for adding
genes and terms. If this is needed, then we recommend creating a new
Ontology or manipulating the contents in a different library, such as
NetworkX or igraph, and transforming the results into Ontology.

.. code:: 

    # Renaming genes and terms.
    ont2 - ont.rename(genes-{'A' : 'A_alias'}, terms-{'S0':'S0_alias'})
    ont2.to_table()




.. raw:: html

    <div>
    <style>
	.dataframe thead tr:only-child th {
	    text-align: right;
	}

	.dataframe thead th {
	    text-align: left;
	}

	.dataframe tbody tr th {
	    vertical-align: top;
	}
    </style>
    <table border-"1" class-"dataframe">
      <thead>
	<tr style-"text-align: right;">
	  <th></th>
	  <th>Parent</th>
	  <th>Child</th>
	  <th>EdgeType</th>
	</tr>
      </thead>
      <tbody>
	<tr>
	  <th>0</th>
	  <td>S2</td>
	  <td>S5</td>
	  <td>Child-Parent</td>
	</tr>
	<tr>
	  <th>1</th>
	  <td>S2</td>
	  <td>S6</td>
	  <td>Child-Parent</td>
	</tr>
	<tr>
	  <th>2</th>
	  <td>S1</td>
	  <td>S3</td>
	  <td>Child-Parent</td>
	</tr>
	<tr>
	  <th>3</th>
	  <td>S1</td>
	  <td>S4</td>
	  <td>Child-Parent</td>
	</tr>
	<tr>
	  <th>4</th>
	  <td>S1</td>
	  <td>S5</td>
	  <td>Child-Parent</td>
	</tr>
	<tr>
	  <th>5</th>
	  <td>S0_alias</td>
	  <td>S1</td>
	  <td>Child-Parent</td>
	</tr>
	<tr>
	  <th>6</th>
	  <td>S0_alias</td>
	  <td>S2</td>
	  <td>Child-Parent</td>
	</tr>
	<tr>
	  <th>7</th>
	  <td>S3</td>
	  <td>C</td>
	  <td>Gene-Term</td>
	</tr>
	<tr>
	  <th>8</th>
	  <td>S4</td>
	  <td>C</td>
	  <td>Gene-Term</td>
	</tr>
	<tr>
	  <th>9</th>
	  <td>S3</td>
	  <td>B</td>
	  <td>Gene-Term</td>
	</tr>
	<tr>
	  <th>10</th>
	  <td>S5</td>
	  <td>E</td>
	  <td>Gene-Term</td>
	</tr>
	<tr>
	  <th>11</th>
	  <td>S4</td>
	  <td>D</td>
	  <td>Gene-Term</td>
	</tr>
	<tr>
	  <th>12</th>
	  <td>S6</td>
	  <td>G</td>
	  <td>Gene-Term</td>
	</tr>
	<tr>
	  <th>13</th>
	  <td>S5</td>
	  <td>F</td>
	  <td>Gene-Term</td>
	</tr>
	<tr>
	  <th>14</th>
	  <td>S6</td>
	  <td>H</td>
	  <td>Gene-Term</td>
	</tr>
	<tr>
	  <th>15</th>
	  <td>S3</td>
	  <td>A_alias</td>
	  <td>Gene-Term</td>
	</tr>
      </tbody>
    </table>
    </div>



Delete S1 and G while preserving transitive connections
-------------------------------------------------------

.. code:: 

    ont2 - ont.delete(to_delete-['S1', 'G'])
    print(ont2)


.. parsed-literal::

    7 genes, 6 terms, 8 gene-term relations, 6 term-term relations
    node_attributes: []
    edge_attributes: []


Delete S1 and G (don't preserve transitive connections)
-------------------------------------------------------

.. code:: 

    ont2 - ont.delete(to_delete-['S1', 'G'], preserve_transitivity-False)
    print(ont2)


.. parsed-literal::

    7 genes, 6 terms, 8 gene-term relations, 3 term-term relations
    node_attributes: []
    edge_attributes: []


Propagate gene-term connections
-------------------------------

.. code:: 

    ont2 - ont.propagate(direction-'forward', gene_term-True, term_term-False)
    print(ont2)

    # Remove all transitive connections, and maintain only a parsimonious set of connections
    ont3 - ont2.propagate(direction-'reverse', gene_term-True, term_term-False)


.. parsed-literal::

    8 genes, 7 terms, 27 gene-term relations, 7 term-term relations
    node_attributes: []
    edge_attributes: []


Propagate term-term connections
-------------------------------

.. code:: 

    ont2 - ont.propagate(direction-'forward', gene_term-False, term_term-True)
    print(ont2)

    # Remove all transitive connections, and maintain only a parsimonious set of connections
    ont3 - ont2.propagate(direction-'reverse', gene_term-False, term_term-True)


.. parsed-literal::

    8 genes, 7 terms, 9 gene-term relations, 11 term-term relations
    node_attributes: []
    edge_attributes: []


Take the subbranch consisting of all term and genes under S1
------------------------------------------------------------

.. code:: 

    ont2 - ont.focus(branches-['S1'])
    print(ont2)


.. parsed-literal::

    Genes and Terms to keep: 10
    6 genes, 4 terms, 7 gene-term relations, 3 term-term relations
    node_attributes: ['Original_Size']
    edge_attributes: []


Inferring a data-driven ontology
--------------------------------

An ontology can also be inferred in a data-driven manner based on an
input set of node-node similarities.

.. code:: 

    sim, genes - ont.flatten()
    print(genes)
    print(sim)


.. parsed-literal::

    ['A' 'B' 'C' 'D' 'E' 'F' 'G' 'H']
    [[ 1.41503751  1.41503751  1.41503751  0.41503751  0.41503751  0.41503751 -0.         -0.        ]
     [ 1.41503751  1.41503751  1.41503751  0.41503751  0.41503751  0.41503751 -0.         -0.        ]
     [ 1.41503751  1.41503751  2.          2.          0.41503751  0.41503751 -0.         -0.        ]
     [ 0.41503751  0.41503751  2.          2.          0.41503751  0.41503751 -0.         -0.        ]
     [ 0.41503751  0.41503751  0.41503751  0.41503751  2.          2.          1.          1.        ]
     [ 0.41503751  0.41503751  0.41503751  0.41503751  2.          2.          1.          1.        ]
     [-0.         -0.         -0.         -0.          1.          1.          2.          2.        ]
     [-0.         -0.         -0.         -0.          1.          1.          2.          2.        ]]


.. code:: 

    ont2 - Ontology.run_clixo(sim, 0.0, 1.0, square-True, square_names-genes)

.. code:: 

    print(ont2)


.. parsed-literal::

    8 genes, 7 terms, 9 gene-term relations, 7 term-term relations
    node_attributes: []
    edge_attributes: ['CLIXO_score']


Ontology alignment
------------------

.. code:: 

    ## Make a second ontology

    # Connections from child terms to parent terms
    hierarchy - [('T3', 'T1'),
		 ('T4', 'T1'),
		 ('T1', 'T0'),
		 ('T5', 'T0')]

    # Connections from genes to terms
    mapping - [('A', 'T3'),
	       ('B', 'T3'),
	       ('C', 'T3'),
	       ('D', 'T4'),
	       ('E', 'T4'),
	       ('F', 'T4'),
	       ('G', 'T5'),
	       ('H', 'T5')]

    # Construct ontology
    ont_B - Ontology(hierarchy, mapping)

.. code:: 

    ont.align(ont_B)


.. parsed-literal::

    collapse command: /cellar/users/mikeyu/DeepTranslate/ddot/ddot/alignOntology/collapseRedundantNodes /tmp/tmp7vURpx
    collapse command: /cellar/users/mikeyu/DeepTranslate/ddot/ddot/alignOntology/collapseRedundantNodes /tmp/tmpL2rvmk
    Alignment command: /cellar/users/mikeyu/DeepTranslate/ddot/ddot/alignOntology/calculateFDRs /tmp/tmpbJwfYR /tmp/tmpmzw76c 0.05 criss_cross /tmp/tmpFprCRw 100 40 gene




.. raw:: html

    <div>
    <style>
	.dataframe thead tr:only-child th {
	    text-align: right;
	}

	.dataframe thead th {
	    text-align: left;
	}

	.dataframe tbody tr th {
	    vertical-align: top;
	}
    </style>
    <table border-"1" class-"dataframe">
      <thead>
	<tr style-"text-align: right;">
	  <th></th>
	  <th>Term</th>
	  <th>Similarity</th>
	  <th>FDR</th>
	</tr>
	<tr>
	  <th>Term</th>
	  <th></th>
	  <th></th>
	  <th></th>
	</tr>
      </thead>
      <tbody>
	<tr>
	  <th>S3</th>
	  <td>T3</td>
	  <td>0.985294</td>
	  <td>0.000</td>
	</tr>
	<tr>
	  <th>S1</th>
	  <td>T1</td>
	  <td>0.913608</td>
	  <td>0.000</td>
	</tr>
	<tr>
	  <th>S6</th>
	  <td>T5</td>
	  <td>0.910000</td>
	  <td>0.040</td>
	</tr>
	<tr>
	  <th>S0</th>
	  <td>T0</td>
	  <td>0.892982</td>
	  <td>0.005</td>
	</tr>
      </tbody>
    </table>
    </div>



Construct ontotypes
-------------------

.. code:: 

    # Genotypes can be represented as tuples of mutated genes
    genotypes - [('A', 'B'),
		 ('A', 'E'),
		 ('A', 'H'),
		 ('B', 'E'),
		 ('B', 'H'),
		 ('C', 'F'),
		 ('D', 'E'),
		 ('D', 'H'),
		 ('E', 'H'),
		 ('G', 'H')]

    ontotypes - ont.get_ontotype(genotypes)
    print(ontotypes)


.. parsed-literal::

       S0  S1  S2  S3  S4  S5  S6
    0   0   0   0   2   0   0   0
    1   0   0   0   1   0   1   0
    2   0   0   0   1   0   0   1
    3   0   0   0   1   0   1   0
    4   0   0   0   1   0   0   1
    5   0   0   0   1   1   1   0
    6   0   0   0   0   1   1   0
    7   0   0   0   0   1   0   1
    8   0   0   0   0   0   1   1
    9   0   0   0   0   0   0   2


.. code:: 

    # Genotypes can also be represented a genotype-by-gene matrix
    import pandas as pd, numpy as np
    genotypes_df - pd.DataFrame(np.zeros((len(genotypes), len(ont.genes)), np.float64),
				index-['Genotype%s' % i for i in range(len(genotypes))],
				columns-ont.genes)
    for i, (g1, g2) in enumerate(genotypes):
	genotypes_df.loc['Genotype%s' % i, g1] - 1.0
	genotypes_df.loc['Genotype%s' % i, g2] - 1.0
    print(genotypes_df)

    ontotypes - ont.get_ontotype(genotypes_df, input_format-'matrix')
    print(ontotypes)


.. parsed-literal::

		 A    B    C    D    E    F    G    H
    Genotype0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0
    Genotype1  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0
    Genotype2  1.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0
    Genotype3  0.0  1.0  0.0  0.0  1.0  0.0  0.0  0.0
    Genotype4  0.0  1.0  0.0  0.0  0.0  0.0  0.0  1.0
    Genotype5  0.0  0.0  1.0  0.0  0.0  1.0  0.0  0.0
    Genotype6  0.0  0.0  0.0  1.0  1.0  0.0  0.0  0.0
    Genotype7  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0
    Genotype8  0.0  0.0  0.0  0.0  1.0  0.0  0.0  1.0
    Genotype9  0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0
		S0   S1   S2   S3   S4   S5   S6
    Genotype0  0.0  0.0  0.0  2.0  0.0  0.0  0.0
    Genotype1  0.0  0.0  0.0  1.0  0.0  1.0  0.0
    Genotype2  0.0  0.0  0.0  1.0  0.0  0.0  1.0
    Genotype3  0.0  0.0  0.0  1.0  0.0  1.0  0.0
    Genotype4  0.0  0.0  0.0  1.0  0.0  0.0  1.0
    Genotype5  0.0  0.0  0.0  1.0  1.0  1.0  0.0
    Genotype6  0.0  0.0  0.0  0.0  1.0  1.0  0.0
    Genotype7  0.0  0.0  0.0  0.0  1.0  0.0  1.0
    Genotype8  0.0  0.0  0.0  0.0  0.0  1.0  1.0
    Genotype9  0.0  0.0  0.0  0.0  0.0  0.0  2.0


Conversions to NetworkX and igraph
----------------------------------

.. code:: 

    G - ont.to_igraph()
    print(G)


.. parsed-literal::

    IGRAPH DN-- 15 16 --
    + attr: NodeType (v), name (v), EdgeType (e)
    + edges (vertex names):
    A->S3, B->S3, C->S3, C->S4, D->S4, E->S5, F->S5, G->S6, H->S6, S5->S2, S6->S2,
    S3->S1, S4->S1, S5->S1, S1->S0, S2->S0


.. code:: 

    G - ont.to_networkx()
    print(G.nodes())
    print(G.edges())


.. parsed-literal::

    ['A', 'C', 'B', 'E', 'D', 'G', 'F', 'S3', 'H', 'S1', 'S0', 'S6', 'S5', 'S4', 'S2']
    [('A', 'S3'), ('C', 'S3'), ('C', 'S4'), ('B', 'S3'), ('E', 'S5'), ('D', 'S4'), ('G', 'S6'), ('F', 'S5'), ('S3', 'S1'), ('H', 'S6'), ('S1', 'S0'), ('S6', 'S2'), ('S5', 'S2'), ('S5', 'S1'), ('S4', 'S1'), ('S2', 'S0')]


Visualization in HiView (http://hiview.ucsd.edu)
------------------------------------------------

.. code:: 

    url, _ - ont.to_ndex(ndex_server-ndex_server, ndex_user-ndex_user, ndex_pass-ndex_pass, layout-'bubble-collect')
    print(url)


.. parsed-literal::

    http://dev2.ndexbio.org/v2/network/7c8fc40e-369a-11e8-929a-0660b7976219


After reading this tutorial
---------------------------

You should check out the list of functions of the `Ontology
class <http://the-data-driven-ontology-toolkit-ddot.readthedocs.io/en/latest/ontology.html>`__
and a list of `utility
functions <http://the-data-driven-ontology-toolkit-ddot.readthedocs.io/en/latest/utils.html>`__
that may help you build more concise pipelines
