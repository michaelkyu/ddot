all:
	cd ddot/alignOntology && $(MAKE)
	cd ddot/clixo_0.3 && $(MAKE)

clean:
	cd ddot/alignOntology && $(MAKE) clean
	cd ddot/clixo_0.3 && $(MAKE) clean


