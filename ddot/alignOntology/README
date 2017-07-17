Instructions for Ontology Alignment
Michael Kramer
Ideker Lab - UCSD
Original created: 9/12/11
Last updated: 4/14/14

C++ programs in main directory (alignOntology_0.1) which need to be compiled for whichever machine you are using.  Simply type make from within the main alignOntology_0.1 folder.

alignOntology:

This program is the main program for aligning ontologies.  It takes two ontologies and a couple of settings parameters as arguments and outputs to stdout an alignment file (described later).  It takes the following arguments:
./alignOntology computedOntology referenceOntology minimum_similarity_value_accepted semantic_verification_version allow_multiple_mappings genes_as_terms terminal_category_name

computedOntology and referenceOntology are simply the paths to the files containing the two ontologies to align.  The ontologies must be in source-target format (three columns separated by tabs where the first column is the source node/concept, the second column is the target node/concept/gene, and the third column is the type of relationship between the two items).  All lines where the target is a terminal node (i.e. gene) must have a third column value that equals terminal_category_name (the 7th parameter, default is "gene").  As far as this program is concerned, these ontologies are symmetric.  The only difference is the order of the output columns (this difference is important for some of the wrapper scripts).

minimum_similarity_value_accepted is the mimimum similarity score that will be considered for a possible mapping between two ontologies.  This value can be any real number from 0 to 1.  In practice, 0.05 works well.  Lower values will cause the program to take longer to converge and will give more insignificant mappings.  Higher values will potentially miss some significant mappings.

semantic_verification_version is a setting that will switch the way the program looks for inconsistencies between pairs of mappings.  "criss_cross" will look for parent-child criss-cross mappings (i.e. cases where e1 < e2 and e1' < e2' but we have the mappings (e1,e2') and (e2,e1').  "strict_hierarchy" will consider it an inconsistency if we have the mappings (e1,e1') and (e2,e2') where and either e1 < e2 or e1' < e2' but not both - this is a stricter requirement than criss_cross (note: this does not perform the check if any of e1,e2,e1', or e2' is a gene, as this would lead to an overly strict requirement that two concepts contain identical genes to be mapped to each other).  "none" will cause the program to do neither of the previous methods (it will still look for double mapping inconsistencies). 'sib_sib' will not allow two sibling terms in one ontology to be mapped to parent-child terms in the second ontology.

allow_multiple_mappings will default to 0 (false) if it is not included.  To allow the same term to be mapped twice (only in the case that two possible mappings are mathematically equivalent), simply set this argument to 1.

The output file for this program will contain several lines that detail the state of the alignment at various iterations of aligning and will give the time that various stages took to complete.  After a line which contains only the word "Matched", the final alignment is printed.  Each line of the file after "Matched" contains one mapping between the two alignments.  Column 1 is the node/gene/concept from computedOntology.  Column 2 is the node/gene/concept from referenceOntology.  Column 3 is the similarity value calculated between these two nodes. Only mappings between internal nodes are reported (mappings between terminal nodes are done simply on identity of names and are trivial).

An example output can be seen in test/test_alignment_results_correct/alignment.out.  This is the result of running (from inside the alignOntology_0.1 folder):

make alignOntology
./alignOntology test/computed_ont.txt test/ref_ont.txt 0.05 criss_cross > test/test_alignment_results_correct/alignment.out

addFDRsToAlignment:

This program will add the false discovery rates to the alignment file after an alignment and several random alignments have been run.  It is best to just let this get called from inside the calculateFDRs script.  This program will need to be compiled for the machine you are running it on before running either of the aforementioned scripts.

SCRIPTS

calculateFDRs:

This program is usually the main way to run the alignment.  It will automatically both align two ontologies and calculate false discovery rates for each individual mapping in the resulting alignment.  It takes as input the computed ontology, the reference ontology, the minimum similarity value accepted, the semantic verification mode, a path to the directory where the results will be placed (directory must already exist and be empty), the number of iterations of random alignments to perform for calculating FDRs, the number of threads that the program should use on the machine, the false discovery rate cutoff (used only for filtering at the end, not in actually performing alignments or calculating FDR), and the terminal category identifier (usually "gene").

An example run of this script from inside the test folder is:

./calculateFDRs computed_ont.txt ref_ont.txt 0.01 criss_cross test_alignment_results 100 30 gene &

The output of this program will be several files, contained in the results directory folder (test_alignment_results in the example above). The raw alignment file will be called alignment_without_descendents (format described above in alignOntology description). The file FDRs will contain the calculated false discovery rates for each of these alignments.  The non-commented lines in this file will have 5 columns.  The first three columns are identical to the raw alignment.  The 4th column is the estimated false discovery rate (FDR) for this size node at the alignment score displayed in column 3. The 5th column is the size of the computed ontology term (i.e. number of reachable terminal nodes).

The expected results of running the above script are contained in test/test_alignment_results_correct.
