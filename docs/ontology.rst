Ontology Class
==============

.. autoclass:: ddot.Ontology
	       
Read/write Ontology objects
--------------------------------

.. automethod:: ddot.Ontology.from_table
.. automethod:: ddot.Ontology.to_table
.. automethod:: ddot.Ontology.read_pickle
.. automethod:: ddot.Ontology.to_pickle
.. automethod:: ddot.Ontology.to_ndex
.. automethod:: ddot.Ontology.from_ndex
.. automethod:: ddot.Ontology.to_cx
.. automethod:: ddot.Ontology.to_graphml

NetworkX and igraph
-------------------

.. automethod:: ddot.Ontology.to_networkx
.. automethod:: ddot.Ontology.from_networkx
.. automethod:: ddot.Ontology.from_igraph
.. automethod:: ddot.Ontology.to_igraph	
		
Inspecting structure
-------------------

.. automethod:: ddot.Ontology.connected
.. automethod:: ddot.Ontology.get_best_ancestors
.. automethod:: ddot.Ontology.topological_sorting
		
Manipulating structure
-----------------------

.. automethod:: ddot.Ontology.unfold
.. automethod:: ddot.Ontology.delete
.. automethod:: ddot.Ontology.focus
.. automethod:: ddot.Ontology.propagate

Inferring data-driven ontology
------------------------------

.. automethod:: ddot.Ontology.flatten
.. automethod:: ddot.Ontology.run_clixo
		
Aligning ontologies
-------------------

.. automethod:: ddot.Ontology.align
