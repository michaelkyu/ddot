from __future__ import absolute_import

import itertools, multiprocessing, logging, os, collections, random, math, sys, time
from itertools import groupby, combinations
from operator import *
from collections import Counter
import tempfile
from subprocess import Popen, PIPE, STDOUT
import inspect
import shlex
import shutil
from StringIO import StringIO
            
import numpy as np
import pandas as pd
import pandas.io.pickle
import networkx as nx
import igraph
import scipy, scipy.sparse
from scipy.sparse import csr_matrix
from scipy.stats import hypergeom

import ndex.client as nc
from ndex.networkn import NdexGraph
        
import ddot
import ddot.config
from ddot.utils import time_print, set_node_attributes_from_pandas, set_edge_attributes_from_pandas, nx_to_NdexGraph, NdexGraph_to_nx, parse_ndex_uuid, parse_ndex_server, make_index, update_nx_with_alignment, bubble_layout_nx, split_indices_chunk, invert_dict, make_network_public, nx_edges_to_pandas, nx_nodes_to_pandas

def _collapse_node(g,
                   v,
                   edge_filter=None,
                   use_v_name=False,
                   combine_attrs=None,
                   default_attr=None,
                   verbose=True,
                   fast_collapse=False,
                   delete=True):
    """Collapses a node in a Graph (igraph package) while preserving
    long-range hierarchical relations between descendants and
    ancestral nodes.

    """

    if use_v_name:
        assert isinstance(v, (unicode, str))
        v = g.vs.find(name_eq=v).index

    if fast_collapse:
        parents = g.neighbors(v, mode='out')
        children = g.neighbors(v, mode='in')

        if len(parents) > 0 and len(children) > 0:
            # A faster collapse that adds all new edges
            # simultaneously. Ignores edge attributes
            new_edges = [(c, p) for p in parents for c in children]
            new_edges = [x for x, y in zip(new_edges, g.get_eids(new_edges, error=False)) if y == -1]
            g.add_edges(new_edges)
    else:
        in_edges = g.es[g.incident(v, mode='in')]
        out_edges = g.es[g.incident(v, mode='out')]

        if edge_filter is not None:
            in_edges = [e for e in in_edges if edge_filter(e)]
            out_edges = [e for e in out_edges if edge_filter(e)]

        for e_in in in_edges:
            for e_out in out_edges:

                in_neigh, out_neigh = e_in.source, e_out.target

                # Only add an edge if it doesn't already exist                                                                                                                   
                if g[in_neigh, out_neigh] == 0:
                    g.add_edge(in_neigh, out_neigh)
                    e = g.es[g.get_eid(in_neigh, out_neigh)]
                    if combine_attrs is not None:
                        # Set default value of edge attributes to 0                                                                                                              
                        for key in combine_attrs:  e[key] = None

                e = g.es[g.get_eid(in_neigh, out_neigh)]

                # Update attributes                                                                                                                                              
                if combine_attrs is not None:
                    for key in combine_attrs:
                        e[key] = combine_attrs[key](e_in, e_out, e)
                        if verbose and key=='triangle_edge_priority':
                            print('Setting',
                                  key,
                                  g.vs[in_neigh]['name'],
                                  g.vs[out_neigh]['name'],
                                  'to',
                                  combine_attrs[key](e_in, e_out, e),
                                  (e_in[key], e_out[key]))

                e['collapsed_length'] = e_in['collapsed_length'] + e_out['collapsed_length']

    e['collapsed_terms'] = e_in['collapsed_terms'] + [g.vs[v]['name']] + e_out['collapsed_terms']

    if delete:
        g.delete_vertices(v)

    return g

def read_alignment_file(f, source='Term_1'):
    """Parses an alignment file created from alignOntology's calculateFDRs script

    Parameters
    -----------
    f : str
       Filename of alignment file

    source : str 
       Indicates which ontology will be the index of the
       returned pandas.DataFrame. Value must be either 'Term_1' (first
       ontology) or 'Term_2' (second ontology)

    Returns
    --------
    : pandas.DataFrame
       DataFrame with four columns: 'Term', 'Similarity', 'FDR', and 'Size'.
       The index of the DataFrame are the names of terms in the "source" ontology. 
    
    """

    # Five columns in the input file
    # 1) Term from first "computed" ontology
    # 2) Term from second "reference" ontology
    # 3) Similarity value
    # 4) FDR
    # 5) Size of the term in the first ontology

    df = pd.read_table(f,
                       names=['Term_1', 'Term_2', 'Similarity', 'FDR', 'Size'],
                       dtype={'Term_1':str,
                              'Term_2':str,
                              'Similarity':np.float64,
                              'FDR':np.float64,
                              'Size':np.int64},
                       header=None)
    target = 'Term_2' if source=='Term_1' else 'Term_1'
    df.rename(columns={target : 'Term'}, inplace=True)
    df.set_index(source, inplace=True)
    df.index.rename('Term', inplace=True)
    return df

def align_hierarchies(hier1,
                      hier2,
                      iterations,
                      threads,
                      update_hier1=False,
                      update_hier2=False,
                      calculateFDRs=None,
                      mutual_collapse=True,
                      output=None):
    if output is None:
        with tempfile.NamedTemporaryFile('w', delete=True) as output_file:
            return align_hierarchies(hier1, hier2, iterations, threads,
                                     update_hier1=update_hier1, update_hier2=update_hier2,
                                     mutual_collapse=mutual_collapse,
                                     output=output_file.name,
                                     calculateFDRs=calculateFDRs)

    hier1_orig, hier2_orig = hier1, hier2
    if mutual_collapse:
        hier1, hier2 = Ontology.mutual_collapse(hier1, hier2)

    def to_file(hier):
        if isinstance(hier, Ontology):
            with tempfile.NamedTemporaryFile('w', delete=False) as f:
                hier.to_3col_table(f, parent_child=True)
            hier = f.name
        else:
            assert isinstance(hier, file) or os.path.exists(hier)
        return hier

    hier1 = to_file(hier1)
    hier2 = to_file(hier2)

    if calculateFDRs is None:
        top_level = os.path.dirname(os.path.abspath(inspect.getfile(ddot)))
        calculateFDRs = os.path.join(top_level, 'alignOntology', 'calculateFDRs')

        #assert os.path.isdir(ddot.config.alignOntology)
        #calculateFDRs = os.path.join(ddot.config.alignOntology, 'calculateFDRs')
    assert os.path.isfile(calculateFDRs)

    output_dir = tempfile.mkdtemp(prefix='tmp')
    cmd = '{5} {0} {1} 0.05 criss_cross {2} {3} {4} gene'.format(
              hier1, hier2, output_dir, iterations, threads, calculateFDRs)
    print 'Alignment command:', cmd

    p = Popen(shlex.split(cmd), shell=False)
    
    try:        
        p.wait()
        shutil.copy(os.path.join(output_dir, 'alignments_FDR_0.1_t_0.1'), output)
    finally:
        if os.path.isdir(output_dir):
            shutil.rmtree(output_dir)

        if p.poll() is None:
            if verbose: time_print('Killing alignment process %s. Output: %s' % (p.pid, output))
            p.kill()  # Kill the process

        align1 = read_alignment_file(output)[['Term', 'Similarity', 'FDR']]
        align2 = align1.copy()
        align2.index, align2['Term'] = align2['Term'].values.copy(), align2.index.values.copy()
        
        append_prefix = lambda x: 'Aligned_%s' % x

        if update_hier1:
            try:
                node_attr = hier2_orig.node_attr[update_hier1]
            except KeyError:
                node_attr = hier2_orig.node_attr

            hier2_import = pd.merge(pd.DataFrame(index=align2.index), node_attr, left_index=True, right_index=True, how='left')
            assert (hier2_import.index == align2.index).all()
            # Change index to terms in hier1
            hier2_import.index = align2['Term'].copy()
            hier2_import.rename(columns=append_prefix, inplace=True)
            
        if update_hier2:
            try:
                node_attr = hier1_orig.node_attr[update_hier2]
            except KeyError:
                node_attr = hier1_orig.node_attr
                
            hier1_import = pd.merge(pd.DataFrame(index=align1.index), node_attr, left_index=True, right_index=True, how='left')
            assert (hier1_import.index == align1.index).all()
            # Change index to terms in hier2
            hier1_import.index = align1['Term'].copy()
            hier1_import.rename(columns=append_prefix, inplace=True)
            
        if update_hier1:
            hier1_orig.update_node_attr(align1.rename(columns=append_prefix))
            hier1_orig.update_node_attr(hier2_import)
        if update_hier2:
            hier2_orig.update_node_attr(align2.rename(columns=append_prefix))
            hier2_orig.update_node_attr(hier1_import)

        return align1
        
def read_term_descriptions(description_table):
    """Read mapping of GO term ID to their term descriptions

    Parameters
    ----------
    description_table : str
       Filename of a tab-delimited table with two columns: GO term ID and term description

    Returns
    -------
    : pandas.Series
       The index is the GO ID and the values are the term descriptions
    
    """

    go_descriptions = pd.read_table(description_table,
                                    header=None,
                                    names=['Term', 'Term_Description'],
                                    index_col=0)
    go_descriptions = go_descriptions['Term_Description']
    assert go_descriptions.index.duplicated().sum() == 0

    return go_descriptions

def parse_obo(obo,
              output_file=None,
              id2name_file=None,
              id2namespace_file=None,
              alt_id_file=None):
    """Parses an OBO file and writes the results into several tables.

    Parameters
    ----------
    obo : str

        Filename of OBO file

    output_file : str

        Filename to write table that describes the ontology's
        hierarchical structure. The table has four columns: (1) parent
        term, (2) child term, (3) relation type (e.g. "is_a" or
        "part_of"), (4) namespace of relation
        (e.g. "biological_process" or "cellular component")

    id2name_file : str

        Filename to write table of term descriptions.  The table has
        two columns: (1) Ontology term (e.g. "GO:0000030"), (2)
        description (e.g. "mannosyltransferase activity")

    id2namespace_file : str
    
        Filename to write table of term namespaces.  The table has two
        columns: (1) Ontology term (e.g. "GO:0000030"), (2) namespace
        of the term (e.g. "biological_process")

    alt_id_file : str
    
        Filename to write table of alternative Term IDs that are
        synonyms and refer to the same term. The table has two
        columns: (1) Primary Term ID, (2) Alternative Term ID

    """

    ## Keywords that screw up parsing:
    # import, is_anonymous, intersection_of, union_of

    ## Relations
    # 'is_a:'
    # 'relationship: has_part'  # Not in filtered GO
    # 'relationship: occurs_in' # Not in filtered GO
    # 'relationship: part_of'   
    # 'relationship: positively_regulates' 
    # 'relationship: negatively_regulates'
    # 'relationship: regulates'
    # 'relationship: results_in' # Not in filtered GO

    stanza, edges = [], []
    id2name = dict()
    id2namespace = dict()
    alt_id = dict()
    in_term_stanza = False
    default_namespace_exists = False
    for line in open(obo).read().splitlines():

        line = line.split('!')[0].strip()  # Remove comments

        if len(line)>0 and line[0]=='[' and line[-1]==']':
            # Add last stanza if it was a term stanza.  Include namespace.
            if in_term_stanza:
                edges.extend(x+(namespace, ) for x in stanza)

            # Start new term stanza
            stanza = []
            
            # Set the default namespace, if it exists
            if default_namespace_exists:
                namespace = default_namespace
            
            # In a term stanzo or not
            in_term_stanza = line =='[Term]'

            name = None
                
        #if 'alt_id:' in line: assert False

        if 'id:' == line[:3]:
            curr_term = line.split('id:')[1].strip()
        elif 'alt_id:' in line:
            alt_term = line.split('alt_id:')[1].strip()
            if alt_id.has_key(curr_term):  alt_id[curr_term].append(alt_term)
            else:                          alt_id[curr_term] = [alt_term]
            id2name[alt_term] = name
        elif 'name:' in line:
            name = line.split('name:')[1].strip()
            assert not id2name.has_key(curr_term)
            id2name[curr_term] = name
        elif 'is_a:' in line:
            parent = line.split('is_a:')[1].strip()
            stanza.append((parent, curr_term, 'is_a'))
        elif 'relationship:' in line:
            line = line.split('relationship:')[1].strip().split()
            if len(line)!=2: print line
            assert len(line)==2
            relation, parent = line
            stanza.append((parent, curr_term, relation))
        elif 'namespace:' == line[:10]:
            namespace = line.split('namespace:')[1].strip()
            assert not id2namespace.has_key(curr_term)
            id2namespace[curr_term] = namespace
        elif 'default-namespace:' == line[:18]:
            namespace = line.split('default-namespace:')[1].strip()
            default_namespace_exists = True
            default_namespace = namespace

    pd.DataFrame(edges).to_csv(output_file, header=False, index=False, sep='\t')
    pd.Series(id2name).to_csv(id2name_file, sep='\t')
    pd.Series(id2namespace).to_csv(id2namespace_file, sep='\t')
    pd.Series(dict([(a, c) for a, b in alt_id.items() for c in b])).to_csv(alt_id_file, sep='\t')

def parse_gaf(gaf):
    """
    Read gene-term annotations from GAF file format:

    http://geneontology.org/page/go-annotation-file-gaf-format-21

    Parameters
    ----------
    gaf : str
        Filename of GAF file

    Returns
    --------
    A list of 2-tuples (gene, GO term)
    """

    gaf_columns = ['DB', 'DB Object ID', 'DB Object Symbol',
                   'Qualifier', 'GO ID', 'DB:Reference',
                   'Evidence Code', 'With (or) From', 'Aspect',
                   'DB Object Name', 'DB Object Synonym',
                   'DB Object Type', 'Taxon', 'Date',
                   'Assigned By', 'Annotation Extension',
                   'Gene Product Form ID']
    df = pd.read_table(gaf, header=None, comment='!', names=gaf_columns)

    # Check that all annotations are to UniProtKB protein IDs
    # assert df['DB'].unique().size == 1 and df['DB'].unique()[0]=='UniProtKB'    

#    return df.loc[:, ['DB Object ID', 'GO ID']].values.tolist()
    return df

class Ontology:
    """A Python representation for constructing, analyzing, and
    manipulating the hierarchical structure of ontologies.

    An Ontology object contains the following attributes for
    inspecting the hierarchical structure. Do not directly modify
    these attributes.

    Parameters
    ----------
    genes : list
        Names of genes
    
    terms : list
        Names of terms

    gene_2_term : dict
 
        gene_2_term[<gene>] --> list of terms that <gene> is mapped
        to. Terms are represented as their 0-based index in self.terms

    term_2_gene : dict

        term_2_gene[<term>] --> list of genes that <term>
        contains. Genes are represented as their 0-based index in self.genes

    child_2_parent : dict

        child_2_parent[<child>] --> list of the parent terms of <child>

    parent_2_child : dict

        parent_2_child[<parent>] --> list of the children terms of <parent>

    term_sizes : list

        A list of every term's size, i.e. the number of genes it
        contains. This list has the same order as self.terms.

        Equivalently, `term_sizes[i] = len(self.term_2_gene[self.terms[i]])`.
    """

    NODETYPE_ATTR = 'NodeType'
    GENE_NODETYPE = 'Gene'
    TERM_NODETYPE = 'Term'

    EDGETYPE_ATTR = 'EdgeType'
    GENE_TERM_EDGTYPE = 'Gene-Term'
    CHILD_PARENT_EDGTYPE = 'Child-Parent'

    def __init__(self,
                 hierarchy,
                 mapping,
                 edge_attr=None,
                 node_attr=None,
                 parent_child=False,
                 add_root_name=None,
                 propagate=None,
                 verbose=True):
        """Construct an Ontology object.

        Parameters
        ----------
        hierarchy : list, tuple

            Iterable of (child term, parent term). E.g. list of 2-tuples

        mapping : list, tuple

            Iterable of (gene, term) pairs. E.g. list of 2-tuples

        edge_attr : pandas.DataFrame

            Meta-data describing (child_term, parent_term)
            pairs. Suggestion: The index of the DataFrame must be a
            pandas.MultiIndex, where the first level is the child term
            and the second level is the parent term.

        parent_child : bool

            If True, then the definitions of <hierarchy> and <mapping>
            are reversed so that they iterate over (parent term, child
            term) and (term, gene) pairs.
        
        propagate : None, str

            The direction ('forward' or 'reverse') to propagate
            gene-term annotations up the hierarchy with
            Ontology.propagate_annotations(). If None, then don't
            propagate annotations.

        add_root_name : bool

            The name of an artificial root. If there are multiple
            roots in the ontology, then they are joined into one root
            with this name. Default: Don't create this root.

        """

        if parent_child:
            hierarchy = [(x[1],x[0]) for x in hierarchy]
            mapping = [(x[1],x[0]) for x in mapping]
            
        ## Read term-to-term edges        
        # parent_2_child[<term_name>] --> list of <term_name>'s children terms
        self.parent_2_child = {r: tuple(p[0] for p in q) for r, q in \
                           itertools.groupby(sorted(hierarchy,
                                                    key=lambda a:a[1]),
                                             key=lambda a:a[1])}

        if add_root_name is not None:
            ## Check if there is a single unifying root term of the
            ## ontology. If not, then identify the multiple roots and
            ## join them under an artificial root
            root_list = self.get_roots()
            print 'Unifying %s roots into one super-root' % len(root_list)
            if len(root_list) > 1:
                root_name = add_root_name
                self.parent_2_child[root_name] = root_list

        ## Read gene-to-term edges
        # self.gene_2_term[<gene_name>] --> list of terms that <gene_name> is mapped to
        self.gene_2_term = {key: tuple(set([a[1] for a in group])) for key, group in \
                            itertools.groupby(sorted(mapping,
                                                     key=lambda a:a[0]),
                                              key=lambda a:a[0])}

        ## Check that the set of terms is the same according to
        ## parent_2_child and self.gene_2_term
        terms_A = set.union(set(self.parent_2_child.keys()),
                            *[set(x) for x in self.parent_2_child.values()])
        terms_B = set.union(*[set(x) for x in self.gene_2_term.values()])

        if verbose and len(terms_B - terms_A)>0:
            print 'WARNING: There are {} terms that are annotated to genes but not connected to the rest of the ontology'.format(len(terms_B - terms_A))
        if verbose and len(terms_A - terms_B)>0:
            print 'WARNING: There are {} terms that have no direct gene annotations'.format(len(terms_A - terms_B))
            
        self.terms = sorted(list(terms_A | terms_B))
        self.genes = sorted(self.gene_2_term.keys())
        
        ## terms_index[<term_name>] --> index in self.terms
        self.terms_index = make_index(self.terms)

        ## self.genes_index[<gene_name>] --> index in self.genes        
        self.genes_index = make_index(self.genes)

        ## Convert self.gene_2_term to list term indices rather than term names
        for k, v in self.gene_2_term.items():
            self.gene_2_term[k] = [self.terms_index[x] for x in self.gene_2_term[k]]

        if node_attr is None:
            self.clear_node_attr()
        else:
            assert node_attr.index.nlevels == 1
            if node_attr.index.name != 'Node':
                if verbose:
                    print "Changing node_attr index name from %s to 'Node'" % node_attr.index.name
                node_attr.index.name = 'Node'                
            self.node_attr = node_attr

        if edge_attr is None:
            self.clear_edge_attr()
        else:
            assert edge_attr.index.nlevels == 2
            if edge_attr.index.names != ['Child', 'Parent']:
                if verbose:
                    print "Changing edge_attr index names from %s to ['Child', 'Parent']" % edge_attr.index.names
                edge_attr.index.names = ['Child', 'Parent']
            self.edge_attr = edge_attr

        self._update_fields()

        if propagate:
            self.propagate_annotations(direction=propagate, inplace=True)
            self._update_fields()

    def _update_fields(self):
        self.child_2_parent, self.child_2_parent_indices = self._get_child_2_parent()
        self.term_2_gene = self._get_term_2_gene()
        self.term_sizes = self.get_term_sizes()

    def _get_child_2_parent(self):
        """
        Read term-to-term edges
        # child_2_parent[<term_name>] --> list of <term_name>'s parent term names
        # child_2_parent_indices[<term_name>] --> list of indices of <term_name>'s parent terms
        """
    
        cp_pairs = []
        for p, c_list in self.parent_2_child.items():
            for c in c_list:
                cp_pairs.append((c,p))
        first = lambda a: a[0]
        cp_pairs.sort(key=first)

        child_2_parent = {
            r: tuple(p[1] for p in q) for r, q in 
            itertools.groupby(cp_pairs, key=first)
        }

        child_2_parent_indices = {
            r : [self.terms_index[x] for x in v] 
            for r, v in child_2_parent.items()
        }

        return child_2_parent, child_2_parent_indices

    def clear_node_attr(self):
        """Resets the node attributes to be empty.""" 

        self.node_attr = pd.DataFrame()
        self.node_attr.index.name = 'Node'

    def clear_edge_attr(self):
        """Resets the edge attributes to be empty."""
        
        self.edge_attr = pd.DataFrame()
        self.edge_attr.index = pd.MultiIndex(levels=[[],[]],
                                             labels=[[],[]],
                                             names=['Child', 'Parent'])

    def update_node_attr(self, node_attr):
        """Update existing node attributes or add new node attributes.

        Parameters
        ----------
        node_attr : pandas.DataFrame

            Dataframe where index are the names of genes or terms and
            where the columns are the names of node attributes.

        """
        
        ####
        # TODO : make sure that renaming/deleting/collapsing of genes and columns respect the node_attr and edge_attr

        # Filter for genes and terms in the ontology
        nodes = set(self.genes) | set(self.terms)
        node_attr = node_attr.loc[[x for x in node_attr.index if x in nodes], :]

        assert node_attr.index.duplicated().sum() == 0
        
        # Update index
        self.node_attr = self.node_attr.reindex(self.node_attr.index.union(node_attr.index))

        # Update columns
        for col in node_attr.columns:
            self.node_attr.loc[node_attr.index, col] = node_attr[col].values

    def update_edge_attr(self, edge_attr):
        """Update existing edge attributes or add new edge attributes.

        Parameters
        ----------
        edge_attr : pandas.DataFrame

            Dataframe where the index is a MultiIndex represents edges
            in the Ontology, such that the first level is the name of
            a gene or child term, and the second level is the name of
            a parent term. Columns are the names of edge attributes.

        """

        # Filter for genes and terms in the ontology
        edges = []
        for child, parent_list in self.child_2_parent.items():
            for parent in parent_list:
                edges.append((child, parent))
        for gene, term_list in self.gene_2_term.items():
            for term in term_list:
                edges.append((gene, self.terms[term]))
        edges = set(edges)
        edge_attr = edge_attr.loc[[x for x in edge_attr.index if x in edges], :]

        assert edge_attr.index.duplicated().sum() == 0
        
        # Update index
        self.edge_attr = self.edge_attr.reindex(self.edge_attr.index.union(edge_attr.index))

        # Update values for overlapping columns
        for col in edge_attr.columns:
            self.edge_attr.loc[edge_attr.index, col] = edge_attr[col].values
        
    def get_roots(self):
        """Returns a list of the root term(s).

        Returns
        -------
        : list

        """

        return sorted(set(self.parent_2_child.keys()) - set([y for x in self.parent_2_child.values() for y in x]))

    def align(self,
              hier,
              iterations,
              threads,
              update_self=False,
              update_ref=False,
              calculateFDRs=None,
              mutual_collapse=True,
              output=None):
        """Identifies one-to-one matches between terms in this ontology with
        highly similar terms in another ontology.


        This function wraps around the C++ code in the alignOntology
        package by Michael Kramer at
        https://github.com/mhk7/alignOntology

        Reference:

        Dutkowski, J., Kramer, M., Surma, M.A., Balakrishnan, R.,
        Cherry, J.M., Krogan, N.J. and Ideker, T., 2013. "A gene
        ontology inferred from molecular networks." *Nature
        biotechnology*, 31(1).

        Parameters
        ----------
        hier : ddot.Ontology.Ontology

            The ontology to align against.

        iterations : int
        
            The number of null model randomizations to create FDR score.

        threads : int

            Number of CPU processes to run simultaneously. Used to
            parallelize the the null model randomizations.

        update_self : bool

            If True, then import the node attributes from the
            reference hierarchy as attributes in this hierarchy
        
        update_ref : bool

            If True, then import the node attributes from the this
            hierarchy as attributes in the reference hierarchy

        mutual_collapse : bool
        
            If True, then remove genes that are unique to either
            ontology, and then remove redundant terms in both
            ontologies using Ontology.collapse_ontology().

        calculate_FDRs : str

            Filename of the 'calculateFDRs' scripts in the
            alignOntology C++package at
            https://github.com/mhk7/alignOntology. Default: use the
            'calculateFDRs' script that comes built-in with ddot.

        output : str

            Filename to write the results of the alignment as a
            tab-delimited file. Default: don't write to a file

        Returns
        -------
        : pandas.DataFrame

            Dataframe where index are names of terms in this
            ontology. There are three columns: 'Term' (name of the
            aligned term), 'Similarity' (the similarity score for the
            alignment), 'FDR' (the FDR of this alignment given the
            null models).

        """
        
        return align_hierarchies(
            self,
            hier,
            iterations,
            threads,
            update_hier1=update_self,
            update_hier2=update_ref,
            calculateFDRs=calculateFDRs,
            mutual_collapse=mutual_collapse,
            output=output)

    def _make_dummy(self):
        ont = self.propagate_annotations(direction='reverse', inplace=False)

        new_gene_2_term = []
        new_child_2_parent = []
        for t in ont.terms:
            if len(ont.term_2_gene[t]) > 0:
                dummy_term = 'dummy_%s' % t
                for g in ont.term_2_gene[t]:
                    new_gene_2_term.append([ont.genes[g], dummy_term])
                new_child_2_parent.append([dummy_term, t])
            else:
                for g in ont.term_2_gene[t]:            
                    new_gene_2_term.append([ont.genes[g], term])
            for p in ont.child_2_parent.get(t, []):
                new_child_2_parent.append([t, p])

        ont_dummy = Ontology(new_child_2_parent, new_gene_2_term)
        return ont_dummy

    def to_networkx(self,
                    layout='bubble',
                    spanning_tree=True):

        """Converts Ontology into a NetworkX object.

        Parameters
        ----------
        
        node_attr : pandas.DataFrame
            
            Meta-data about genes and terms that will be included as node
            attributes in the NetworkX object.

        edge_attr : pandas.DataFrame
            
            Meta-data about connections among genes and terms that
            will be included as edge attributes in the NetworkX
            object.

        spanning_tree : bool
        
            If True, then identify a spanning tree of the DAG. include
            an edge attribute "Is_Tree_Edge" that indicates

        layout : str

            The name of the layout algorithm for laying out the
            Ontology as a graph. Node positions are astored in the
            node attributes 'x_pos' and 'y_pos'. If None, then do not
            perform a layout.

        Returns:
        --------
        : nx.DiGraph

        """

        G = nx.DiGraph()

        #################################
        ### Add nodes and node attributes

        G.add_nodes_from([g for g in self.genes])
        G.add_nodes_from([t for t in self.terms])

        term_sizes = dict(zip(self.terms, self.get_term_sizes(propagate=True)))

        for t in self.terms:
            G.node[t][self.NODETYPE_ATTR] = self.TERM_NODETYPE
            G.node[t]['Size'] = term_sizes[t]
            G.node[t]['isRoot'] = False
        for g in self.genes:
            G.node[g][self.NODETYPE_ATTR] = self.GENE_NODETYPE
            G.node[g]['Size'] = 1
            G.node[g]['isRoot'] = False

        # Identify the root
        root = self.get_roots()[0]
        G.node[root]['isRoot'] = True

        set_node_attributes_from_pandas(G, self.node_attr)
        set_edge_attributes_from_pandas(G, self.edge_attr)
            
        #################################
        ### Add edges and edge attributes

        G.add_edges_from([(g, self.terms[t],
                           {self.EDGETYPE_ATTR : self.GENE_TERM_EDGETYPE}) \
                          for g in self.genes for t in self.gene_2_term[g]])

        G.add_edges_from([(c, p,  
                           {self.EDGETYPE_ATTR : self.CHILD_PARENT_EDGETYPE}) \
                          for p in self.terms for c in self.parent_2_child.get(p, [])])

        if spanning_tree:
            # Identify a spanning tree
            tree_edges = self.get_tree_edges()
            nx.set_edge_attributes(G,
                                   'Is_Tree_Edge',
                                   {(s,t) : 'Tree' if ((s,t) in tree_edges) else 'Not_Tree'
                                    for s, t in G.edges_iter(data=False)})

            if layout:
                if layout=='bubble':
                    ont_dummy = self._make_dummy()
                    G_tree = ont_dummy.to_networkx(layout=None, spanning_tree=False)
                    for u, v, data in G.edges(data=True):
                        if data['Is_Tree_Edge']=='Not_Tree':
                            if data[self.EDGETYPE_ATTR]==self.GENE_TERM_EDGETYPE:
                                G_tree.remove_edge(u, 'dummy_' + v)
                            else:
                                G_tree.remove_edge(u, v)

                    G.pos = bubble_layout_nx(G_tree)
                    G.pos = {n : p for n, p in G.pos.items() if 'dummy_' not in n}
                    nx.set_node_attributes(G, 'x_pos', {n : x for n, (x,y) in G.pos.items()})
                    nx.set_node_attributes(G, 'y_pos', {n : y for n, (x,y) in G.pos.items()})
                else:
                    raise Exception('Unsupported layout: %s', layout)

        return G

    @classmethod
    def from_pandas(cls,
                    table,
                    parent_col='Parent',
                    child_col='Child',
                    is_mapping=None,
                    propagate=False,
                    verbose=False):
        """Constructs an Ontology from a pandas Dataframe describing the edges.

        TODO: consider deleting this method because it is redundant
        with from_table. The only non-redundant aspect is that the
        default is_mapping function is different from from_table.
        
        Parameters
        -----------
        table : pandas.Dataframe

        parent : int or str

            Column for parent terms (index or name of column)
        
        child : int or str

            Column for child terms and genes (index or name of column)

        is_mapping : function
        
            Function applied on each row to determine if it represents
            a gene-term annotation. If the function returns True, then
            row represents a (gene, term) pair. If False, it
            represents a (child term, parent term) pair. The function
            is evaluated using df.apply(is_mapping, axis=1).

        propagate : str

            The direction ('forward' or 'reverse') for propagating
            gene-term annotations up the hierarchy with
            Ontology.propagate_annotations(). If None, then don't
            propagate annotations.

        Returns
        -------
        : ddot.Ontology.Ontology

        """

        if is_mapping is None:
            is_mapping = lambda x: x[self.EDGETYPE_ATTR]==self.GENE_TERM_EDGETYPE
            
        return cls.from_table(df,
                              is_mapping=is_mapping,
                              parent_col=parent_col,
                              child_col=child_col,
                              mapping_file=None,
                              propagate=propagate,
                              verbose=verbose)

    @classmethod
    def from_table(cls,
                   table,
                   is_mapping=None,
                   parent=0,
                   child=1,                   
                   mapping=None,
                   mapping_parent=0,
                   mapping_child=0,
                   propagate=False,
                   verbose=False,
                   **kwargs):
        """Create Ontology from a tab-delimited table.

        Parameters
        ----------
        table : pandas.DataFrame, file-like object, or filename

            A table that lists (child term, parent term) pairs. If
            mapping is None, then this table should also list (gene,
            term) pairs.

        is_mapping : function

            Function applied on each row to determine if it represents
            a gene-term annotation. If the function returns True, then
            row represents a (gene, term) pair. If False, it
            represents a (child term, parent term) pair. Only applied
            when mapping=None.

        parent : int or str

            Column for parent terms (index or name of column)
        
        child : int or str

            Column for child terms and genes (index or name of column)

        propagate : str

            The direction ('forward' or 'reverse') for propagating
            gene-term annotations up the hierarchy with
            Ontology.propagate_annotations(). If None, then don't
            propagate annotations.

        mapping : pandas.DataFrame, file-like object, or filename (optional)

            A separate table listing only gene-term annotations

        parent : int or str

            Column for parent terms (index or name of column)
        
        child : int or str

            Column for child terms and genes (index or name of column)
        
        Returns
        -------
        : ddot.Ontology.Ontology

        """

        if is_mapping is None:
            is_mapping = lambda x: x[2]==self.GENE_TERM_EDGETYPE,
            
        # Read table
        try:
            table = pd.read_table(table, comment='#', header=None)
        except:
            assert isinstance(table, pd.DataFrame)

        for col in [child, parent]:
            table.loc[:,col] = table.loc[:,col].astype(str)

        edge_attr = table.set_index([child, parent])
        edge_attr.index.rename(['Child', 'Parent'], inplace=True)

        if mapping is None:
            # Extract gene-term annotations from table
            mask = table.apply(is_mapping, axis=1)
            mapping = table.loc[mask, :].loc[:,[child, parent]].values.tolist()
            hierarchy = table.loc[~mask, :].loc[:,[child, parent]].values.tolist()
        else:
            # Read separate table of gene-term annotations
            try:
                mapping = pd.read_table(mapping, comment='#', header=None)
            except:
                assert isinstance(mapping, pd.DataFrame)

            for col in [mapping_child, mapping_parent]:
                mapping.loc[:,col] = mapping.loc[:,col].astype(str)
                
            mapping_attr = mapping.set_index([mapping_child, mapping_parent])
            mapping_attr.index.rename(['Child', 'Parent'], inplace=True)

            edge_attr = pd.concat([edge_attr, mapping_attr])
            mapping = mapping.loc[:,[mapping_child, mapping_parent]].values.tolist()
            hierarchy = table.loc[:,[child, parent]].values.tolist()

        return cls(hierarchy,
                   mapping,
                   parent_child=False,
                   edge_attr=edge_attr,
                   propagate=propagate,
                   verbose=verbose,
                   **kwargs)

    @classmethod
    def from_ndex(cls,
                  ndex_uuid,
                  ndex_server=None,
                  ndex_user=None,
                  ndex_pass=None,
                  edgetype_attr=None,
                  edgetype_value=None):
        """Reads an Ontology stored on NDEx. Gene and terms are distinguished
        according by an edge attribute.

        Parameters
        ----------
        ndex_uuid : str
            NDEx UUID of ontology

        edgetype_attr : str

            Name of the edge attribute that distinguishes a (gene,
            term) pair from a (child term, parent term) pair

        gene_value : str
        
            Value of the edge attribute for (gene, term) pairs

        Returns
        -------
        : ddot.Ontology.Ontology

        """

        if ndex_server is None:
            ndex_server = ddot.config.ndex_server
        if ndex_user is None:
            ndex_pass = ddot.config.ndex_user
        if ndex_pass is None:
            ndex_pass = ddot.config.ndex_pass
            
        if '/' in ndex_uuid:
            ndex_server = parse_ndex_server(ndex_uuid)
            ndex_uuid = parse_ndex_uuid(ndex_uuid)

        G = NdexGraph(
            server=ndex_server, 
            username=ndex_user,
            password=ndex_pass,
            uuid=ndex_uuid)
        
        return cls.from_NdexGraph(
            G,
            edgetype_attr=edgetype_attr,
            edgetype_value=edgetype_value)        

    @classmethod
    def from_NdexGraph(cls,
                       G,
                       edgetype_attr=None,
                       edgetype_value=None):
        """Converts a NdexGraph object to an Ontology object. Gene and terms
        are distinguished by an edge attribute.

        Parameters
        ----------
        G : NdexGraph

        edgetype_attr : str

            Name of the edge attribute that distinguishes a (gene,
            term) pair from a (child term, parent term) pair

        gene_value : str
        
            Value of the edge attribute for (gene, term) pairs

        Returns
        -------
        : ddot.Ontology.Ontology

        """
        
        return cls.from_networkx(
            NdexGraph_to_nx(G),
            edgetype_attr=edgetype_attr,
            edgetype_value=edgetype_value)

    @classmethod
    def from_networkx(cls,
                      G,
                      edgetype_attr=None,
                      edgetype_value=None):
        """Converts a NetworkX object to an Ontology object. Gene and terms
        are distinguished by an edge attribute.

        Parameters
        ----------
        G : nx.DiGraph

        edgetype_attr : str

            Name of the edge attribute that distinguishes a (gene,
            term) pair from a (child term, parent term) pair

        gene_value : str
        
            Value of the edge attribute for (gene, term) pairs

        Returns
        -------
        : ddot.Ontology.Ontology

        """
        
        if edgetype_attr is None:
            edgetype_attr=self.EDGETYPE_ATTR
        if edgetype_value is None:
            edgetype_value=self.GENE_TERM_EDGETYPE

        hierarchy = []
        mapping = []
        for u, v, attr in G.edges_iter(data=True):
            if attr[edgetype_attr] == edgetype_value:
                mapping.append((u, v))
            else:
                hierarchy.append((u, v))

        edge_attr = nx_edges_to_pandas(G)
        node_attr = nx_nodes_to_pandas(G)
        
        return cls(hierarchy,
                   mapping,
                   node_attr=node_attr,
                   edge_attr=edge_attr)
    
    def collapse_ontology(self,                          
                          method='mhkramer',                          
                          min_term_size=2,
                          default_relation='default',
                          collapseRedundantNodes=None,
                          verbose=True):
        """Remove redundant and empty terms. When a term T is removed,
        hierarchical relations are preserved by connecting every child
        of T with every parent of T. This removal operation has the
        nice property of being commutative, i.e. the order of removal
        does not matter.

        Parameters
        -----------
        method : str

            If "mhkramer", then use the collapseRedundantNodes script
            in the alignOntology package. If "python", then use an
            internal Python script.

        min_term_size : int
        
            Remove terms that are below this size. TODO: not yet supported

        Returns
        -------
        : ddot.ddot.Ontology

            A new Ontology object

        """

        if method=='mhkramer':

            # Propagate forward and then reverse
            ont = self.copy()
            ont = self.propagate_annotations(direction='forward', inplace=False)
            ont.propagate_annotations(direction='reverse', inplace=True)

            if collapseRedundantNodes is None:
                assert os.path.isdir(ddot.config.alignOntology)
                collapseRedundantNodes = os.path.join(ddot.config.alignOntology, 'collapseRedundantNodes')
            assert os.path.isfile(collapseRedundantNodes)

            with tempfile.NamedTemporaryFile('w', delete=False) as f:
                ont.to_3col_table(f, parent_child=True, encoding=None, default_relation=u'default')
            try:                
                cmd = '%s %s' % (collapseRedundantNodes, f.name)
                print 'collapse command:', cmd
                p = Popen(shlex.split(cmd), shell=False, stdout=PIPE, stderr=PIPE)
                collapsed, err = p.communicate()
            finally:
                os.remove(f.name)

            ont = Ontology.from_table(
                StringIO(collapsed),
                is_mapping=lambda x: x[2]=='gene'
            )
            ont.clear_edge_attr()

#            ont.update_node_attr(self.node_attr)
#            ont.update_edge_attr(self.edge_attr)

            return ont

        elif method=='python':

            raise Exception('Unsupported')
        
            g = self.to_igraph()
            if verbose: print len(g.vs), 'total nodes'

            if verbose: print 'Propagating annotations'
            self.propagate_annotations()

            parity = True
            while True:
                
                # Calculate a unique hash for every term based on its set of gense
                names_2_idx = make_index(g.vs['name'])
                term_hash = {names_2_idx[t] : (len(g_list), hash(tuple(g_list)))
                             for t, g_list in self.term_2_gene.items()
                             if names_2_idx.has_key(t)}

                if verbose: time_print('Identify nodes to collapse')
                node_order = g.topological_sorting(mode='out')

                small_terms = [v for v in node_order if term_hash[v][0]<min_term_size]

                # Terms that have the same set of genes as all of its parents

                same_as_parents = []
                same_as_children = []
                for v in node_order:
                    if (len(g.neighbors(v, mode='out')) > 0
                        and all(term_hash[v]==term_hash[y] for y in g.neighbors(v, mode='out'))):
                        same_as_parents.append(v)
                        
                # same_as_all_parents = [v for v in node_order
                #                        if (len(g.neighbors(v, mode='out'))>0 and all(term_hash[v]==term_hash[y] for y in g.neighbors(v, mode='out')))]
                # same_as_all_children = [v for v in node_order
                #                         if (len(g.neighbors(v, mode='in'))>0 and all(term_hash[v]==term_hash[y] for y in g.neighbors(v, mode='in')))]

                if verbose:
                    time_print('%s empty terms, %s (%s) terms that are redundant with all their parents (children)' %
                               (len(small_terms), len(same_as_parents), len(same_as_children)))

                to_delete = list(set(small_terms) | set(same_as_children if parity else same_as_parents))

                if verbose:
                    time_print('Collapsing %s empty terms and %s terms that redundant with its %s' % \
                               (len(small_terms),
                                len(same_as_all_children) if parity else len(same_as_all_parents),
                                'children' if parity else 'parents'))
                parity = not parity

                if len(to_delete)==0:
                    break
                else:
                    for v in to_delete:
                        g = _collapse_node(g,
                                           v,
                                           use_v_name=False,
                                           verbose=False,
                                           fast_collapse=True,
                                           delete=False)
                    g.delete_vertices(to_delete)

            g.es(relation_eq=None)['relation'] = default_relation

            remaining_terms = set([self.terms_index[x] for x in g.vs['name']])
            return Ontology(g,
                            mapping_file=[(gene, self.terms[t]) 
                                          for gene, t_list in self.gene_2_term.items() 
                                          for t in t_list if t in remaining_terms],
                            combined_file=False, parent_child=False)
    
    @classmethod
    def mutual_collapse(cls, ont1, ont2, verbose=False):
        """Collapses two ontologies to the common set of genes.

        Parameters
        -----------
        ont1 : ddot.Ontology.Ontology

        ont2 : ddot.Ontology.Ontology

        Returns
        -------
        ont1_collapsed : ddot.Ontology.Ontology

        ont2_collapsed : ddot.Ontology.Ontology

        """

        common_genes = set(ont1.genes) & set(ont2.genes)

        if verbose:
            print('Common genes:', len(common_genes))

        ont1 = ont1.delete_genes(set(ont1.genes) - common_genes, inplace=False)
        ont1_collapsed = ont1.collapse_ontology(method='mhkramer')
        ont2 = ont2.delete_genes(set(ont2.genes) - common_genes, inplace=False)
        ont2_collapsed = ont2.collapse_ontology(method='mhkramer')

        if verbose:
            print 'ont1_collapsed:', ont1_collapsed.summary()
            print 'ont2_collapsed:', ont2_collapsed.summary()

        return ont1_collapsed, ont2_collapsed

    def delete_terms(self, terms, inplace=False):
        """Delete genes or terms from the ontology.

        Note that if a gene is directly connected to a term T but not
        to an ancestor of T, then deleting T will lose the information
        about the indirect connection to the ancestor. For this
        reason, it is recommended that propagate_annotations() is
        called first.

        Parameters
        ----------
        terms : iterable of str
           
            Terms to delete

        genes : iterable of str

            Genes to delete

        inplace : bool
        
            If True, then modify the ontology. If False, then create and modify a copy.

        Returns
        -------
        : ddot.Ontology.Ontology

        """
        
        if inplace:
            ont = self
        else:
            ont = self.copy()
        
        terms = set(terms)
        tmp_gene_2_term = {g : [ont.terms[t] for t in t_list]
                           for g, t_list in ont.gene_2_term.items()}
        ont.terms = [t for t in ont.terms if t not in terms]
        ont.terms_index = make_index(ont.terms)
        ont.gene_2_term = {g : [ont.terms_index[t] for t in t_list if t not in terms]
                           for g, t_list in tmp_gene_2_term.items()}
        ont.parent_2_child = {p : [c for c in c_list if c not in terms]
                              for p, c_list in ont.parent_2_child.items()
                              if p not in terms}
        
        genes = set(genes)
        ont.genes = [g for g in ont.genes if g not in genes]
        ont.genes_index = make_index(ont.genes)
        ont.gene_2_term = {g : t for g, t in ont.gene_2_term.items() 
                           if g not in genes}
        ont._update_fields()

        return ont
        
    def rename(self,
               genes={},
               terms={},
               inplace=False):
        """Rename gene and/or term names. 

        TODO: update node and edge attr

        Parameters
        ----------
        genes : dict or function

           If dictionary, then it maps current gene names to new
           names. Genes not in dictionary are deleted.

           If function, then genes(name) returns the new name.

        terms : dict or function

           If dictionary, then it maps current term names to new
           names. Terms not in dictionary are deleted.

           If function, then terms(name) returns the new name.

        inplace : bool

           If True, then modify the ontology. If False, then create
           and modify a copy.

        Returns
        -------
        : ddot.Ontology.Ontology

        """

        try:
            terms = {t : terms(t) for t in self.terms}
        except:
            pass

        try:
            genes = {g : genes(g) for g in self.genes}
        except:
            pass

        if inplace:
            ont = self
        else:
            ont = self.copy()

        if genes:
            ont.delete_genes([g for g in ont.genes if not genes.has_key(g)], inplace=True)

            new_genes = []
            new_gene_2_term = {}
            for g in ont.genes:
                new_g = genes[g]
                if isinstance(new_g, (unicode, str)):
                    new_genes.append(new_g)
                    new_gene_2_term[new_g] = ont.gene_2_term[g]
                elif hasattr(new_g, '__iter__'):
                    for new_gg in new_g:
                        new_genes.append(new_gg)
                        new_gene_2_term[new_gg] = ont.gene_2_term[g]
            ont.genes = new_genes
            ont.gene_2_term = new_gene_2_term
            ont.genes_index = make_index(ont.genes)
            ont._update_fields()
        if terms:
            ont.delete_terms([t for t in ont.terms if not terms.has_key(t)], inplace=True)
            ont.terms = [terms[t] for t in ont.terms]
            ont.terms_index = make_index(ont.terms)
            ont.parent_2_child = {terms[p] : [terms[c] for c in c_list]
                                  for p, c_list in ont.parent_2_child.items()}
            ont._update_fields()

        conversions = genes.copy()
        conversions.update(terms)
        
        # Update node and edge attributes
        index = ont.node_attr.index
        ont.node_attr.index = pd.Series(index).replace(conversions)

        index = ont.edge_attr.index
        ont.edge_attr.index = index.to_frame().replace({'Child' : conversions, 'Parent' : conversions}).set_index(['Child', 'Parent']).index
        
        return ont

    def to_pandas(self,
                  term_2_term=True,
                  gene_2_term=True,
                  default_relation=u'default'):
        """Convert Ontology to pandas.DataFrame
        
        Parameters
        ----------
        term_2_term : bool

            Include (child term, parent term) pairs

        gene_2_term : bool

            Include (gene, term) pairs

        default_relation : str
        
            The relation type assigned to all (child term, parent term) pairs

        Returns
        -------
        : pandas.DataFrame

            Three columns: (1) Parent (2) Child -- a gene or child
            term, (3) EdgeType type.  The relation type of all (gene,
            term) pairs is "gene".

        """

        df = pd.DataFrame(columns=['Parent','Child','EdgeType'])
        if term_2_term:
            df = df.append(self._hierarchy_to_pandas(default_relation),
                           ignore_index=True)
        if gene_2_term:
            tmp = self._mapping_to_pandas()
            tmp.rename(columns={'Gene':'Child', 'Term':'Parent'},
                       inplace=True)
            tmp['EdgeType'] = 'gene'
            df = df.append(tmp, ignore_index=True)
        return df

    def _hierarchy_to_pandas(self, default_relation=u'default'):

        triples = [(p,c, relation_dict.get((c, p), default_relation)) \
                   for p, c_list in self.parent_2_child.items() for c in c_list]
        df = pd.DataFrame(triples, columns=['Parent', 'Child', 'EdgeType'])
        return df

    def _mapping_to_pandas(self):
        pairs = [(g, self.terms[t]) for g, t_list in self.gene_2_term.items() for t in t_list]
        df = pd.DataFrame(pairs, columns=['Gene', 'Term'])        
        return df

    def to_3col_table(self,
                      output,
                      header=False,
                      parent_child=True,
                      encoding=None,
                      default_relation=u'default'):
        """Write Ontology to tab-delimited table with 3 columns and 2 types of rows:

        Row Type A. Parent term, child term, ontology relation 
           e.g. "GO:0007005      GO:0007006      is_a"
        Row Type B. Term, gene, "gene"
           e.g. "GO:0007005      ATG7    gene"

        Parameters
        ----------
        output : str
            Filename of table

        header : bool
            If True, include a header line in the table

        parent_child : bool
            If False, then the first and second columns switch places

        encoding
            Character encoding

        default_relation : str
            
            EdgeType type for all (child term, parent term)
            pairs. Note that the relation type for (gene, term) pairs
            is "gene".

        """

        df = self.to_pandas(default_relation=default_relation)
        if parent_child:
            df = df[['Parent','Child','EdgeType']]
        else:
            df = df[['Child','Parent','EdgeType']]
        return df.to_csv(output, header=header, index=False, sep='\t')

    def copy(self):
        """Create a deep copy of the Ontology object"""

        hierarchy = [(c, p) for p, c_list in self.parent_2_child.items() for c in c_list]
        mapping = [(g, self.terms[t]) for g, t_list in self.gene_2_term.items() for t in t_list]
        return Ontology(hierarchy,
                        mapping,
                        edge_attr=None if (self.edge_attr is None) else self.edge_attr.copy(),
                        node_attr=None if (self.node_attr is None) else self.node_attr.copy(),
                        parent_child=False,
                        verbose=False)

    def semantic_similarity(self,
                            genes_subset=None,
                            term_sizes='subset',
                            between_terms=False,
                            output='Resnik'):
        """Computes the semantic similarity between pair of genes in
        <genes_subset>. Similarity s(g1,g2) is defined as
        :math:`-log_2(|T_sca| / |T_root|)` where :math:`|T|` is the number of genes in
        <genes_subset> that are under term T. :math:`T_sca` is the "smallest
        common ancestor", the common ancestral term with the smallest
        term size. :math:`T_root` is the root term of the ontology.

        Parameters
        -----------
        genes_subset : iterable

            The set of genfes, over which pairs the similarity will be
            calculated. If <term_sizes>='subset', then term sizes will
            be recalculated according to only these genes, rather than
            all genes in the ontology

        between_terms : bool
            if True, then output similarity between all terms

        output : str
            type of semantic similarity

        """

        # If no genes are specified, then compute the similarity between all pairs of genes
        if genes_subset is None: genes_subset = self.genes

        time_print('Calculating term sizes')
        if term_sizes=='all_genes':
            term_sizes = np.array(self.get_term_sizes())
        elif term_sizes=='subset':
            ### Reimplement this by deleting genes

            # Recompute term sizes, with respect to the intersection of genes
            term_2_gene = self.term_2_gene
            genes_subset_set = set(genes_subset)
            term_sizes = np.array([len(set([self.genes[x] for x in term_2_gene[t]]) & genes_subset_set) for t in self.terms])
        else:
            raise Exception()

        if output=='Resnik':
            graph = self.to_igraph()

            sca = get_smallest_ancestor(graph, term_sizes)
            ss_terms = -1 * np.log2(term_sizes[sca] / float(term_sizes.max()))

            if between_terms:
                return ss_terms
            else:
                # For convenience below, set the similarity between a term and itself to be 0
                ss_terms[[np.arange(ss_terms.shape[0]), np.arange(ss_terms.shape[0])]] = 0                
                
                idx_list = [self.gene_2_term[g] for g in genes_subset]
                ss = [ss_terms[idx1, :][:, idx2].max() for idx1, idx2 in combinations(idx_list, 2)]
                return ss

            # # Add nodes in the igraph object to represent genes
            # graph.add_vertices(self.genes)
            # graph.add_edges([(g, t) for g, t_list in self.gene_2_term.items() for t in t_list])
            # assert graph.vs[-len(self.genes):]['name'] == self.genes

            # sca = get_smallest_ancestor(graph, term_sizes)
            # # Offset the indices
            # idx = [len(self.terms) + self.genes_index[g] for g in genes_subset]
            # ss = (-1 * np.log2(term_sizes / float(term_sizes.max())))[sca[idx, :][:, idx][np.triu_indices(len(idx), k=1)]]

            # return ss
        elif output=='sca_list':
            ## For each pair of gene, return a list of the smallest
            ## common ancestors (sca). There may be more than one sca with the same size.

            gene_2_term_numpy = {g : np.array(t_list) for g, t_list in self.gene_2_term.items()}            
            common_ancestors = [np.intersect1d(gene_2_term_numpy[g1], gene_2_term_numpy[g2], assume_unique=True) \
                                for g1, g2 in combinations(genes_subset, 2)]
            assert all(x.size > 0 for x in common_ancestors)
            min_size = [term_sizes[x].min() for x in common_ancestors]
            sca_list = [x[term_sizes[x]==m] for x, m in zip(common_ancestors, min_size)]

            # Dict: (g1,g2) gene pairs --> list of term indices
            return {(g1,g2) : x for (g1, g2), x in zip(combinations(genes_subset, 2), sca_list)}
        
    def _get_term_2_gene(self, verbose=False): 
        if verbose: print 'Calculating term_2_gene'

        term_2_gene = invert_dict(
            self.gene_2_term,
            keymap=make_index(self.genes),
            valmap=dict(enumerate(self.terms)))

        for t in self.terms:
            if not term_2_gene.has_key(t):
                term_2_gene[t] = []
        return term_2_gene

    def get_term_sizes(self, propagate=False):
        """Returns an array of term sizes in the same order as self.terms"""

        if propagate:
            ont = self.propagate_annotations(verbose=False, inplace=False)
        else:
            ont = self

        tmp = Counter([x for y in ont.gene_2_term.values() for x in y])
        term_sizes = [tmp[x] for x in range(len(ont.terms))]
        return term_sizes

    def get_information_gain(self):
        for p in terms:
            self.parent_2_children[p]

    def shuffle_genes(self, inplace=False):
        """Shuffle the names of genes"""

        new_order = self.genes.copy()
        random.shuffle(new_order)
        rename = dict(zip(self.genes, new_order))

        return self.rename(rename, inplace=False)

    def get_tree_edges(self):
        """Identify a spanning tree of the DAG (including genes as part of the
        DAG), and return a list of (u, v) edges in the tree.

        """

        tree = self.to_igraph(include_genes=True, spanning_tree=True)

        tree_edges = set([(tree.vs[e.source]['name'],
                           tree.vs[e.target]['name']) 
                          for e in tree.es if e['Is_Tree_Edge']=='Tree'])
        return tree_edges

    def is_dag(self):
        """Return True if the Ontology is a valid directed acyclic graph,
        False otherwise.

        """

        return self.to_igraph().is_dag()

    def to_igraph(self, include_genes=False, spanning_tree=True):
        """Convert Ontology to an igraph.Graph object. Gene and term names are
           stored in the 'name' vertex attribute of the igraph object.

        Parameters
        ----------
        include_genes : bool

           Include genes as vertices in the igraph object.

        spanning_tree : bool

            If True, then identify a spanning tree of the DAG. include
            an edge attribute "Is_Tree_Edge" that indicates

        Returns
        -------
        : igraph.Graph

        """

        if include_genes:
            terms_index_offset = {t : v + len(self.genes) for t, v in self.terms_index.items()}
            edges = ([(self.genes_index[g], terms_index_offset[self.terms[t]])
                      for g in self.genes
                      for t in self.gene_2_term[g]] +
                     [(terms_index_offset[c], terms_index_offset[p]) 
                      for p, children in self.parent_2_child.items()
                      for c in children])
            graph = igraph.Graph(n=len(self.genes) + len(self.terms),
                                 edges=edges,
                                 directed=True,
                                 vertex_attrs={
                                     'name':self.genes + self.terms,
                                     self.NODETYPE_ATTR:[self.GENE_NODETYPE for x in self.genes] + [self.TERM_NODETYPE for x in self.terms]
                                 })
        else:
            edges = [(self.terms_index[c], self.terms_index[p]) for p, children in self.parent_2_child.items() for c in children]
            graph = igraph.Graph(n=len(self.terms),
                                 edges=edges,
                                 directed=True,
                                 vertex_attrs={'name':self.terms})
        if spanning_tree:
            tmp = self.get_term_sizes(propagate=True)
            parent_priority = [tmp[self.terms_index[v['name']]] if self.terms_index.has_key(v['name']) else 1 for v in graph.vs]

            # Identify spanning tree
            graph = self._make_tree_igraph(
                graph,
                parent_priority=parent_priority,
                optim=min,
                edge_name='Is_Tree_Edge')
            graph.es['Is_Tree_Edge'] = ['Tree' if x else 'Not_Tree' for x in graph.es['Is_Tree_Edge']]

        return graph

    def get_shortest_paths(self, sparse=False, chunk_size=500):
        """Calculate the length of the shortest paths between all pairs of
        terms.

        Parameters
        ----------
        sparse

            If True, return a scipy.sparse matrix. If False, return a
            NumPy array

        chunk_size

            Computational optimization: shortest paths are calculated in batches.

        Returns
        -------
        d : np.ndarray or scipy.sparse.matrix

            d[x,y] is the length of the shortest directed path from a
            descendant term with index x to an ancestral term with
            index y. Term indices are defined by
            self.terms_index. d[x,y]=0 if no directed path exists.

        """
        
        graph = self.to_igraph()

        tmp = [graph.shortest_paths(
                  graph.vs[x[0]:x[1]],
                  graph.vs,
                  mode='out')
               for x in split_indices_chunk(len(graph.vs), chunk_size)]

        if sparse:
            return scipy.sparse.vstack([scipy.sparse.csr_matrix(x) for x in tmp])
        else:
            return np.vstack(tmp, order='C')
        
    def get_longest_paths(self):
        """Computes the lengths of the longest directed paths between all pairs
        of terms.

        Returns
        -------
        : np.ndarray

           NumPy array d where d[x,y] is length of the longest
           directed path from a descendant term x to an ancestral term y

        """

        graph = self.to_igraph()

        return -1 * np.array(graph.shortest_paths(graph.vs, graph.vs, weights=[-1 for x in graph.es], mode='out'), order='C')

    def get_connectivity_matrix(self, sparse=False):
        """Calculate which terms are descendants/ancestors of other terms

        Parameters
        -----------
        sparse : bool

            If True, return a scipy.sparse matrix. If False, return a
            NumPy array

        Creates a term-by-term matrix d where

        Returns
        -------
        d : np.ndarray or scipy.sparse.matrix

            d[i,j] is 1 if term i is an ancestor of term j, and 0
            otherwise. Note that d[i,i] == 1 and d[root,i] == 0, for
            every i.

        """

        time_print('Calculating connectivity matrix')
        if sparse:
            paths = self.get_shortest_paths(sparse=True)
            d = scipy.sparse.coo_matrix((np.isfinite(paths[paths.nonzero()]),
                                              (paths.nonzero()[0], paths.nonzero()[1])),
                                             dtype=np.int8)
        else:
            d = np.int8(np.isfinite(self.get_shortest_paths()))
        return d

    def get_connectivity_matrix_nodiag(self):
        """Returns the same matrix as Ontology.get_connectivity_matrix(),
        except the diagonal of the matrix is set to 0.

        Note: 
            d[a, a] == 0 instead of 1

        Returns
        -------
        : np.ndarray

        """
        
        d = self.get_connectivity_matrix(sparse=False)
        
        d[np.diag_indices(d.shape[0])] = 0
        assert not np.isfortran(d)
        return d

    def get_leaves(self, terms_list, children_list=None):
        """Returns terms in <terms_list> that are not ancestors of any term in
        <children_list>
        
        If <children_list> is None, then select the terms in
        <terms_list> that are not ancestors of any of the other terms
        in <terms_list>.

        Parameters
        ----------
        terms_list : list

        children : list
        
        """
        
        connectivity_matrix_nodiag = self.get_connectivity_matrix_nodiag()
        
        terms_list = np.array(terms_list)
        if children_list is None:
            children_list = terms_list
        else:
            children_list = np.array(children_list)

        return terms_list[~ np.any(connectivity_matrix_nodiag[children_list, :][:, terms_list], axis=0)]

    def propagate_annotations(self,
                              direction='forward',
                              verbose=False,
                              inplace=True):
        """Propagate gene-term annotations through the ontology.

        Example) Consider an ontology with one gene g and three terms t1, t2, t3. The connections are

        Child Term - Parent Term relations:
        t1 --> t2
        t2 --> t3

        Gene - Term relations:
        g --> t1
        g --> t2

        In 'forward' propagation, a new relation g --> t3 is added. In
        'reverse' propagation, the relation g --> t2 is deleted
        because it is an indirect relation inferred from g --> t1 and
        t1 --> t2.

        TODO: consider 'reverse' instead of 'reverse'
        
        Parameters
        ----------
        direction : str

            The direction of propgation. Either 'forward' or 'reverse'

        inplace : bool

            If True, then modify the ontology. If False, then create
            and modify a copy.

        Returns
        -------
        : ddot.Ontology.Ontology

        """
            
        if inplace:
            ont = self
        else:
            ont = self.copy()

        if direction=='forward':
            child_2_parent_idx = {
                ont.terms_index[c] : [ont.terms_index[p] for p in p_list]
                for c, p_list in ont.child_2_parent.items()
            }
            gene_2_term_set = {
                g : set(t_list)
                for g, t_list in ont.gene_2_term.items()
            }

            genes_to_update = set(ont.gene_2_term.keys())
            count = 0
            while len(genes_to_update) > 0:
                # Iterate over a copy of genes_to_update
                for g in genes_to_update.copy():
                    curr_terms = gene_2_term_set[g]
                    num_old = len(curr_terms)
                    curr_terms.update(
                        set([p for t in curr_terms 
                             for p in child_2_parent_idx.get(t, [])]))
                    if len(curr_terms) == num_old:
                        genes_to_update.remove(g)                        

                if verbose: print count,
                count +=1
                if count == 1000:
                    raise Exception('ERROR: Ontology depth >1000. Stopping in case of bug in code')
            if verbose: print
            ont.gene_2_term = {g : sorted(t_set) for g, t_set in gene_2_term_set.items()}
            
            self._update_fields()

        elif direction=='reverse':
            ont.propagate_annotations(direction='forward', inplace=True)

            term_2_gene_set = {}
            for t, g in ont.term_2_gene.items():
                term_2_gene_set[t] = set(g)

            graph = ont.to_igraph(spanning_tree=False)

            for c_idx in graph.topological_sorting(mode='in'):
                child = graph.vs[c_idx]['name']
                for parent in ont.child_2_parent.get(child, []):
                    term_2_gene_set[parent] -= term_2_gene_set[child]

            ont.gene_2_term = invert_dict(term_2_gene_set,
                                          keymap=make_index(ont.terms),
                                          valmap=dict(enumerate(ont.genes)))
            ont.term_2_gene = {a : list(b) for a, b in term_2_gene_set.items()}
        else:
            raise Exception('Unsupported propagation direction: %s' % direction)

        return ont

    def propagate_ontotypes(self, ontotypes, prop, ontotype_size, max_ontotype, method='fixed_size'):
        """Propagates a list of base ontotypes.
        
        Parameters
        -----------
        ontotypes

            An array of ontotypes in scipy.sparse.csr_matrix format.
            Each row is a separate ontotype.

        method

            If 'fixed_size', then the sum of values in every ontotype
            is exactly the same, namely <ontotype_size>

        ontotype_size

            If method=='fixed_size', then this is the sum of values in
            every ontotype

        """

        # Just return the array
        if prop=='genes':
            if method=='fixed_size':
                if ontotype_size == 0:
                    return ontotypes

                assert ontotype_size in [1,2]        

                ## Each non-zero term in the ontotype is assumed to be hit
                ## by exactly one gene and that term is the only term the
                ## gene hits
                graph = self.to_igraph()

                time_print('Calculating node ancestry')
                d = self.get_connectivity_matrix()

                assert d.dtype == ontotypes.dtype

                # The usage of np.int8 for <d> requires ontotype size be <= 127
                assert ontotype_size <= 127

                time_print('Formatting ontotype array')
                ontotype_arr = scipy.sparse.csr_matrix(ontotypes)
                nnz = ontotype_arr.nonzero()

                time_print('Sanity checks')

                # # Number of starting perturbations in each ontotype
                # # Currently supporting only 2
                # assert ontotype_size==2

                # Check that the number of non-zero elements in every row is equal ontotype_size
                assert( ((ontotype_arr != 0).sum(axis=1) != ontotype_size).sum() == 0 )

                # Check that the every consecutive block of <ontotype_size> elements is located in a new row a new row
                assert( (nnz[0][0:len(nnz[0]):ontotype_size,] != range(ontotype_arr.shape[0])).sum() == 0 )

                # Check that the non-zero elements were returned in the order of traversing row-by-row
                assert( (np.argsort(np.int64(nnz[0]) * ontotype_arr.shape[1] + np.int64(nnz[1])) != np.array(range(len(nnz[0])))).sum() == 0)

                # Check that the number of non-zero elements is a multiple of ontotype_size
                assert( len(nnz[1]) % ontotype_size == 0 )

                # import pdb
                # pdb.set_trace()

                time_print('Doing propagation')
                # For all ontotypes, calculate the propagation of the
                # first term perturbed as an array, and likewise for the
                # second term.  Then, add the arrays.

                if ontotype_size == 2:
                    # Matrix for propagating the first perturbed term
                    a = np.take(d, nnz[1][0:len(nnz[1]):2], axis=0) * \
                        np.reshape(np.array(ontotype_arr[(nnz[0][0:len(nnz[0]):2],
                                                          nnz[1][0:len(nnz[1]):2])]), (len(nnz[0])/2, 1))

                    # Matrix for propagating the second perturbed term
                    b = np.take(d, nnz[1][1:len(nnz[1]):2], axis=0) * \
                        np.reshape(np.array(ontotype_arr[(nnz[0][1:len(nnz[0]):2],
                                                          nnz[1][1:len(nnz[1]):2])]), (len(nnz[0])/2, 1))

                    # Since the we're propagating genes, we need to limit the ontotype value
                    ontotype_arr = np.minimum(a + b, max_ontotype)

                elif ontotype_size == 1:
                    ontotype_arr = np.take(d, nnz[1], axis=0) * \
                        np.reshape(np.array(ontotype_arr[(nnz[0], nnz[1])]), (len(nnz[0]), 1))

                    assert (ontotype_arr > max_ontotype).sum().sum() == 0
                    ontotype_arr = np.minimum(ontotype_arr, max_ontotype)

                time_print('Done')
                return ontotype_arr

            elif method=='topological':

                assert ontotypes.dtype == np.int8
                if ontotypes.dtype == np.int8: assert max_ontotype <= 127
                
                if scipy.sparse.issparse(ontotypes):
                    ontotypes = ontotypes.toarray(order='F')
                else:
                    ontotypes = np.array(ontotypes, order='F')

                # Topological sort on terms, going bottom-up the ontology
                for i in self.to_igraph().topological_sorting(mode='out'):
                    
                    # Add the values from child to all its parents.
                    # Cap the value at <max_ontotype>.  When
                    # propagating the gene-term annotations in an
                    # ontology, you should set max_ontotype=1.
                    if self.terms[i] in self.child_2_parent_indices.keys():

                        parents = np.array(self.child_2_parent_indices[self.terms[i]])

                        ontotypes[:, parents] = np.minimum(max_ontotype,
                                                           ontotypes[:, parents] + ontotypes[:, i].reshape(ontotypes.shape[0], 1))

                        if not np.all(ontotypes <= max_ontotype):
                            print i
                            assert False
                assert np.all(ontotypes <= max_ontotype)

                return ontotypes

    def get_ontotype(self,
                     genotypes,
                     input_format='gene_list',
                     output_format='sparse',
                     matrix_columns=None):
        """Transform genotypes to ontotypes.

        Parameters
        ----------
        genotypes : list, np.ndarray, scipy.sparse.spmatrix, pd.DataFrame

        input_format : str

            If 'gene_list', then genotypes is interpreted as a list,
            where each element is an individual genotype represented
            as a list of genes. Thus, genotypes is a list of list of
            genes.

            If 'matrix', then genotypes is interpreted as a
            genotype-by-gene matrix, where the value at position (i,j)
            represents the mutation value of gene j in genotype i.

        output_format : str

            If 'sparse', then return a sparse matrix as a
            scipy.sparse.csr_matrix object. 

            If 'dataframe', then return a pandas.DataFrame

        matrix_columns : list
            
            Must be set when input_format is 'matrix' and genotypes is
            a NumPy array or SciPy sparse matrix. In this case,
            matrix_columns represents a list of the genes that are
            represented by the columns in genotypes.

            Not used when input_format is 'gene_list', or when
            input_format is 'matrix' and genotypes is a
            pandas.DataFrame.

        Returns
        -------
        : scipy.sparse.csr_matrix, pandas.DataFrame

            genotype-by-term matrix. The ordering of genotypes is the
            same as the input. The ordering of terms is the same as
            self.terms.

        """

        if input_format=='gene_list':
            gene_2_term = {k: np.array(v) for k, v in self.gene_2_term.items()}
            genotypes_x = [np.concatenate([gene_2_term[g] for g in gset]) if len(gset)>0 else np.array([]) for gset in genotypes]
            indices = np.concatenate(genotypes_x)
            indptr = np.append(0,
                               np.cumsum([gset.size for gset in genotypes_x]))
            data = np.ones((indices.size, ), dtype=np.int64)

            ontotypes = scipy.sparse.csr_matrix(
                (data, indices, indptr),
                (len(genotypes), len(self.terms))
            )
            ontotypes.sum_duplicates()
        
        elif input_format=='matrix':
            if isinstance(genotypes, pd.DataFrame):
                matrix_columns = genotypes.columns
                genotypes = genotypes.values
            elif isinstance(genotypes, np.ndarray) or scipy.sparse.issparse(genotypes):
                assert matrix_columns is not None
            else:
                raise Exception('<genotypes> must be a genotype-by-gene matrix or pd.DataFrame')
            
            contained = np.array([self.genes_index.has_key(g) for g in matrix_columns])
            
            genotypes = scipy.sparse.csc_matrix(genotypes)[:,contained]
            
            annotation_matrix = scipy.sparse.csr_matrix(self.get_annotation_matrix())
            annotation_matrix = scipy.sparse.csr_matrix(annotation_matrix)[contained,:]

            ontotypes = genotypes.dot(annotation_matrix)
        else:
            raise Exception('Invalid input format')

        if output_format=='dataframe':
            ontotypes = pd.DataFrame(ontotypes.toarray(), columns=self.terms)
        elif output_format=='sparse':
            pass
        else:
            raise Exception('Invalid output format')
            
        return ontotypes

    def get_annotation_matrix(self):
        """Returns a gene-by-term matrix stored as a scipy.sparse.coo_matrix
        
        Returns
        -------
        : scipy.sparse.coo_matrix

        """

        # Convert gene names to indices
        gene_2_term = [(self.genes_index[g], t_list) 
                       for g, t_list in self.gene_2_term.items()]

        annotation_matrix = scipy.sparse.coo_matrix(
            ([1 for g, t_list in gene_2_term for t in t_list],
             ([g for g, t_list in gene_2_term for t in t_list],
              [t for g, t_list in gene_2_term for t in t_list])),
            shape=(len(self.genes), len(self.terms)))
        
        return annotation_matrix

    def summary(self):
        """Summarize the Ontology's contents with respect to number of genes,
        terms, and connections.

        Returns
        --------
        : str

        """

        if self.node_attr is None:
            node_attr_names = None
        else:
            node_attr_names = ', '.join(map(str, self.node_attr.columns))

        if self.edge_attr is None:
            edge_attr_names = None
        else:
            edge_attr_names = ', '.join(map(str, self.edge_attr.columns))
            
        summary = ('%s genes, '
                    '%s terms, '
                    '%s gene-term relations, '
                    '%s term-term relations'
                    '\nnode_attributes: %s'
                    '\nedge_attributes: %s') % (
                        len(self.genes),
                        len(self.terms),
                        sum([len(x) for x in self.gene_2_term.values()]),
                        sum([len(x) for x in self.parent_2_child.values()]),
                        node_attr_names,
                        edge_attr_names)
        return summary

    @classmethod
    def run_clixo(cls,
                  graph,
                  alpha,
                  beta,
                  min_dt=-10000000,
                  timeout=100000000,
                  output=None,
                  output_log=None,
                  verbose=True):
        """Runs the CLIXO algorithm and returns the result as an Ontology object.

        Acts as a wrapper for the C++ package at
        https://mhk7.github.io/clixo_0.3/.
        
        Parameters
        ----------

        graph : pandas.DataFrame

           3-column pandas.DataFrame (long format) or square
           pandas.DataFrame (wide format)

        alpha : float

           CLIXO alpha parameter

        beta : float

           CLIXO beta parameter

        min_dt : float

           Minimum similarity score

        timeout : int

           Maximum time (in seconds) allowed to run CLIXO

        output : str

           Filename to write the resulting Ontology as a
           table. Default: don't write to file

        output_log : str

           Filename to write log information from CLIXO. Default:
           don't write to file

        Returns
        --------
        : ddot.Ontology.Ontology

        """

        rerun = False

        if output is None:
            output_file = tempfile.NamedTemporaryFile('w', delete=False)
            output = output_file.name
            if verbose:
                print 'temp output:', output
            rerun, delete_output = True, True

        if not (isinstance(graph, str) and os.path.exists(graph)):

            # if wide is None:
            #     wide = graph.shape[1] == graph.shape[0]
            # if wide:        
            #     df_sq = pd.DataFrame(graph, index=expand, columns=expand)
            #     df = melt_square(df_sq)

            # Write graph into a temporary file.
            # Assumes that <graph> is a list of 3-tuples (parent, child, score)
            with tempfile.NamedTemporaryFile('w', delete=False) as graph_file:
                try:
                    graph.to_csv(graph_file, sep='\t', header=False, index=False)
                except:
                    graph_file.write('\n'.join(['\t'.join([str(x[0]),str(x[1]),str(x[2])]) for x in graph]) + '\n')
            graph = graph_file.name
            if verbose:
                print 'temp graph:', graph
            rerun, delete_graph = True, True
                                         
        if not (isinstance(output_log, str) and os.path.exists(output_log)):
            output_log_file = tempfile.NamedTemporaryFile('w', delete=False)
            output_log = output_log_file.name
            if verbose:
                print 'temp output log:', output_log
            rerun, delete_output_log = True, True

        if rerun:            
            try:
                return cls.run_clixo(
                    graph, alpha, beta,
                    min_dt=min_dt,
                    timeout=timeout,
                    output=output,
                    output_log=output_log,
                    verbose=verbose)
            finally:
                if delete_output:
                    os.remove(output)
                if delete_output_log:
                    os.remove(output_log)
                if delete_graph:
                    os.remove(graph)

                pass

        if verbose:
            time_print('\t'.join(map(str, [graph, alpha, beta, min_dt])))
        
        top_level = os.path.dirname(os.path.abspath(inspect.getfile(ddot)))
        clixo_cmd = os.path.join(top_level, 'mhk7-clixo_0.3-cec3674', 'clixo')

        # For timestamping everyline: awk '{ print strftime("%Y-%m-%d %H:%M:%S"), $0; fflush(); }'
        cmd = ("""{0} {1} {2} {3} | awk""".format(clixo_cmd, graph, alpha, beta) + 
               """ '{if ( $1 ~ /^#/ ) {print "\#", strftime("%Y-%m-%d %H:%M:%S"), $0 ; fflush() } else {print $0}}'""" + 
               """ | tee {}""".format(output_log))
        if verbose:
            print 'CLIXO command:', cmd

        p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, bufsize=1)

        # p = Popen(shlex.split(cmd), shell=False, stdout=PIPE, stderr=STDOUT, bufsize=1)

        # p_cmd = "{0} {1} {2} {3} | awk""".format(clixo_cmd, graph, alpha, beta)
        # p2_cmd = """awk '{if ( $1 ~ /^#/ ) {print "\#", strftime("%Y-%m-%d %H:%M:%S"), $0 ; fflush() } else {print $0}}'"""
        # print 'CLIXO command:', cmd

        curr_dt = None
        start = time.time()

        # Asynchronous readout of the command's output
        while p.poll() is None:

            # Break if passed the maximum processing time
            if time.time() - start > timeout:
                if verbose:
                    time_print(
                        ('Killing process %s (OUT OF TIME). '
                         'Current dt: %s: Output: %s') % (p.pid, curr_dt, output_log))
                break

            line = p.stdout.readline()

            # Remove newline character
            line = line[:-1]

            # Break if the min_dt has been met
            if '# dt: ' in line:
                curr_dt = float(line.split('# dt: ')[1])
                if curr_dt < min_dt:
                    if verbose:
                        time_print(
                            ('Killing process %s (BEYOND MIN THRESHOLD). '
                             'Current dt: %s, min_dt: %s') % (p.pid, curr_dt, min_dt))
                    break

            # If line was empty, then sleep a bit
            if line=='':  time.sleep(0.1)

        if p.poll() is None:
            if verbose: time_print('Killing process %s. Output: %s' % (p.pid, output_log))
            p.kill()  # Kill the process

            # Extract ontology with extractOnt
            p_ext = Popen('%s %s 0 0 %s' % (extract_cmd, output_log, output),
                          shell=True, stdout=PIPE, stderr=STDOUT)
            p_ext.communicate()
        else:
            if verbose: time_print('Extracting by grep -v #')

            # Extract ontology with grep -v '#'
            p_ext = Popen("grep -v '#' %s | grep -v '@' > %s" % (output_log, output), shell=True)
            p_ext.communicate()

        if verbose:
            time_print('Elapsed time (sec): %s' % (time.time() - start))
        
        ont = cls.from_table(
            output,
            is_mapping=lambda x: x[2]=='gene',
            verbose=verbose)
        ont.rename(terms=lambda x: 'CLIXO:%s' % x, inplace=True)

        ont.edge_attr.columns = map(str, ont.edge_attr.columns)
        ont.edge_attr.drop('2', inplace=True, axis=1)
        ont.edge_attr.rename(columns={'3':'CLIXO_score'}, inplace=True)

        if verbose:
            print 'Ontology:', ont.summary()

        return ont

    def to_ndex(self,
                name=None,
                description=None,
                ndex_server=None,
                ndex_user=None,
                ndex_pass=None,
                network=None,
                features=[],
                subnet_max_term_size=None,
                layout='bubble',
                propagate='reverse',
                style=ddot.config.ontology_style,
                public=False,
                verbose=False):
        """Upload an Ontology to NDEx. The Ontology can be preformatted in
        several ways including

        1. Set a name and description of the Ontology
        2. Upload a supporting gene-gene subnetwork for every term in the Ontology
        3. Propagate gene-term annotations
        4. Layout the nodes.
        5. Apply a visual style, e.g. specifying node and edge colors


        Parameters
        ----------
        name : str
        
            Name of Ontology

        description : str

            Description of Ontology

        layout : str

            The name of the layout algorithm for laying out the
            Ontology as a graph. Node positions are astored in the
            node attributes 'x_pos' and 'y_pos'. If None, then do not
            perform a layout.

        style : ndex.networkn.NdexGraph

            The Cytoscape.js visual style on NDEx. Represented using
            CX and stored in an NdexGraph.

        network : pandas.Dataframe

            Dataframe describing gene-gene network from which to
            create subnetworks for every term. To be passed to
            Ontology.upload_subnets_ndex().

        features : list of str
        
            Columns in the gene-gene network to upload. To be passed
            to Ontology.upload_subnets_ndex().

        ndex_server : str
        
             URL of NDEx server

        ndex_user : str

             NDEx username

        ndex_pass : str
        
             NDEx password

        public : bool

            Whether to make the Ontology public on NDEx

        Returns
        -------
        : ndex.networkn.NdexGraph

        """

        if ndex_server is None:
            ndex_server = ddot.config.ndex_server
        if ndex_user is None:
            ndex_pass = ddot.config.ndex_user
        if ndex_pass is None:
            ndex_pass = ddot.config.ndex_pass
            
        if propagate is not None:
            ont = self.propagate_annotations(direction=propagate)
        else:
            ont = self

        if (sim is not None) and (len(features) > 0):
            if subnet_max_term_size is not None:
                terms = [t for t,s in zip(ont.terms, ont.term_sizes) if s <= subnet_max_term_size]
            else:
                terms = ont.terms

            term_2_uuid = ont.upload_subnets_ndex(
                sim,
                features,
                name,
                ndex_server=ndex_server,
                ndex_user=ndex_user,
                ndex_pass=ndex_pass,
                terms=terms,
                public=public,
                verbose=verbose
            )
        else:
            term_2_uuid = {}

        G = ont.to_NdexGraph(
                name=name,
                description=description,
                term_2_uuid=term_2_uuid,
                layout=layout)

        if style:
            import ndex.beta.toolbox as toolbox
            toolbox.apply_network_as_template(G, style)

        ont_url = G.upload_to(ndex_server, ndex_user, ndex_pass)

        if public:
            ont_uuid = parse_ndex_uuid(ont_url)
            if verbose: print 'Making ontology at %s public on NDEX:' % ont_uuid
            make_network_public(ont_uuid,
                                ndex_server=ndex_server, ndex_user=ndex_user, ndex_pass=ndex_pass,
                                timeout=100)
        
        return ont_url, G
            
    def to_NdexGraph(self,
                     name=None,
                     description=None,
                     term_2_uuid=False,
                     layout='bubble'):
        """Formats an Ontology object into a NetworkX object with extra node
        attributes that are accessed by the hierarchical viewer.

        Parameters
        -----------
        name : str
        
            Name of Ontology, as would appear if uploaded to NDEx.

        description : str

            Description of Ontology, as would appear if uploaded to NDEx.

        term_2_uuid : list

            A dictionary mapping a term to a NDEx UUID of a gene-gene
            subnetwork of genes in that term. the UUID will be stored
            in the node attribute 'ndex:internallink'. If uploaded to
            NDEx, then this attribute will provide a hyperlink to the
            gene-gene subnetwork when the term is clicked upon on the
            NDEx page for this ontology.

            This dictionary can be created using
            Ontology.upload_subnets_ndex(). Default: no dictionary.

        layout : str

            Layout the genes and terms in this Ontology. Stored in the
            node attributes 'x_pos' and 'y_pos'. If None, then do not
            perform a layout.

        Returns
        -------
        : ndex.networkn.NdexGraph

        """

        # Convert to NetworkX
        G = self.to_networkx(layout=layout, spanning_tree=True)

        # Set node attribute 'Label'
        for t in self.terms:
            if not G.node[t].has_key('Label'):
                G.node[t]['Label'] = t
        for g in self.genes:
            if not G.node[g].has_key('Label'):
                G.node[g]['Label'] = g

                # Set links to subnetworks supporting each term
        if term_2_uuid:
            for t in self.terms:        
                if term_2_uuid.has_key(t):
                    G.node[t]['ndex:internalLink'] = '[%s](%s)' % (G.node[t]['Label'], term_2_uuid[t])

        G = nx_to_NdexGraph(G)
        if name is not None:
            G.set_name(name)
        if description is not None:
            G.set_network_attribute('Description', description)

        return G

    def _force_directed_layout(self, G):
        """Force-directed layout on only the terms"""

        sub_nx = G.copy()
        sub_nx.remove_edges_from([(u,v) for u,v,attr in sub_nx.edges_iter(data=True) if attr['Is_Tree_Edge']=='Not_Tree'])
        pos = nx.spring_layout(sub_nx, dim=2, k=None,
                               pos=None,
                               fixed=None,
                               iterations=50,
                               weight=None,
                               scale=1.0)

        tmp = np.array([x[0] for x in pos.values()])
        x_min, x_max = tmp.min(), tmp.max()
        tmp = np.array([x[1] for x in pos.values()])
        y_min, y_max = tmp.min(), tmp.max()

        x_scale = 500. / (y_max - y_min)
        y_scale = 500. / (x_max - x_min)
        pos = {a : [b[0] * x_scale, b[1] * y_scale] for a, b in pos.items()}
        return pos

    def upload_subnets_ndex(self,
                            network,
                            features,
                            name,
                            ndex_server=None,
                            ndex_user=None,
                            ndex_pass=None,
                            terms=None,
                            gene_columns=['Gene1', 'Gene2'],
                            propagate='forward',
                            public=False,
                            node_attr=None,
                            verbose=False):
        """For each term in the ontology, upload a subnetwork of interactions
        between the genes in that term to NDEx.

        TODO: instead of specifying gene_columns, add another
        parameter use_index to specify that genes are the network's
        index

        Parameters
        ----------
        network : pandas.Dataframe

            Dataframe describing network

        features : list of str
        
            Columns in network to upload

        name : str
        
            Prefix for the names of all subnetworks

        ndex_server : str
        
             URL of NDEx server

        ndex_user : str

             NDEx username

        ndex_pass : str
        
             NDEx password

        terms : list

            List of terms to upload a subnetwork. Default: upload for
            all terms.

        gene_columns : list

            Columns in network that represent the two genes.

        propagate : str

            The direction ('forward' or 'reverse') to propagate
            gene-term annotations up the hierarchy with
            Ontology.propagate_annotations(). If None, then don't
            propagate annotations.
        
        public : bool

            Whether to make networks public on NDEx
        
        node_attr : pandas.DataFrame

        """

        if ndex_server is None:
            ndex_server = ddot.config.ndex_server
        if ndex_user is None:
            ndex_pass = ddot.config.ndex_user
        if ndex_pass is None:
            ndex_pass = ddot.config.ndex_pass
        
        if propagate:
            ontology = self.copy()
            ontology.propagate_annotations(direction=propagate)
        else:
            ontology = self

        ndex = nc.Ndex(ndex_server, ndex_user, ndex_pass)
        term_2_uuid = {}

        start = time.time()
        g1, g2 = gene_columns[0] + '_lex', gene_columns[1] + '_lex'

        if verbose:
            print 'features:', features
            print 'gene_columns:', gene_columns

        network = network[features + gene_columns].copy()

        # Filter dataframe for gene pairs within the ontology
        genes_set = set(ontology.genes)
        tmp = [x in genes_set and y in genes_set
               for x, y in zip(network[gene_columns[0]], network[gene_columns[1]])]
        network = network.loc[tmp, :]

        for feat in features:
            network[feat] = network[feat].astype(np.float64)

        # Lexicographically sort gene1 and gene2 so that gene1 < gene2
        network[g1] = network[gene_columns].min(axis=1)
        network[g2] = network[gene_columns].max(axis=1)
        network_idx = {x : i for i, x in enumerate(zip(network[g1], network[g2]))}

        if verbose:
            print 'Setup time:', time.time() - start

        # Normalize features into z-scores
        tmp = network[features]
        network[features] = (tmp - tmp.mean()) / tmp.std()

        # Calculate the min/max range of features
        feature_mins = network[features].min().astype(np.str)
        feature_maxs = network[features].max().astype(np.str)

        if terms is None:
            terms = ontology.terms
            
        for t in terms:
            start = time.time()

            genes = [ontology.genes[g] for g in ontology.term_2_gene[t]]
            genes.sort()
            gene_pairs_idx = [network_idx[gp] for gp in itertools.combinations(genes, 2) \
                              if network_idx.has_key(gp)]

            G_nx = nx.from_pandas_dataframe(network.iloc[gene_pairs_idx, :], g1, g2,
                                            edge_attr=features)
            if node_attr is not None:
                set_node_attributes_from_pandas(G_nx, node_attr)

            G = nx_to_NdexGraph(G_nx)
            G.set_name('%s supporting network for %s' % (name, t))
            G.set_network_attribute('Description', '%s supporting network for %s' % (name, t))
            for f in features:
                G.set_network_attribute('%s min' % f, feature_mins[f])
                G.set_network_attribute('%s max' % f, feature_maxs[f])

            start_upload = time.time()
            ndex_url = G.upload_to(ndex_server, ndex_user, ndex_pass)
            term_2_uuid[t] = parse_ndex_uuid(ndex_url)
            upload_time = time.time() - start_upload

            if verbose:
                print(ontology.terms_index[t],
                      'Term:', t,
                      'Gene pairs:', len(gene_pairs_idx),
                      'Genes:', len(genes),
                      'Time:', round(time.time() - start, 4),
                      'Upload time:', round(upload_time, 4))

        if public:            
            to_upload = set(term_2_uuid.values())
            while len(to_upload) > 0:
                completed = set()
                for ndex_uuid in to_upload:
                    try:
                        ndex.make_network_public(ndex_uuid)
                        completed.add(ndex_uuid)
                    except:
                        pass
                to_upload = to_upload - completed
                time.sleep(0.5)

        return term_2_uuid

    def get_best_ancestor_matrix(self, node_order=None, verbose=False):
        """Compute the 'best' ancestor for every pair of terms. 'Best' is
        specified by a ranking of terms. For example, if terms are
        ranked by size, from smallest to largest, then the smallest
        common ancestor is calculated.

        Parameters
        ----------
        node_order : list

           A list of terms, ordered by their rank with the 'best' term at the beginning.

        Returns
        --------
        ancestors : np.ndarray

            ancestors[a,b] = the best common ancestor of terms a and
            b, represented as a 0-based index of self.terms

        """

        graph = self.to_igraph()

        if node_order is None:
            # By default, sort from smallest to largest terms
            node_order = [a[0] for a in sorted(enumerate(self.get_term_sizes()), key=lambda a: a[1])]

        d = np.int8(np.isfinite(np.array(graph.shortest_paths(graph.vs, graph.vs, mode='out'), order='C')))

        ancestor_matrix = np.zeros(d.shape, dtype=np.int32)
        ancestor_matrix.fill(-1)

        if verbose: time_print('Iterating:')
        for idx, i in enumerate(node_order):

            # Note: includes self as a child
            children_list = np.where(d[:,i] == 1)[0]

            # For those descendants without a computed LCA yet, set their LCA to this term
            lca_sub = ancestor_matrix[children_list.reshape(-1,1), children_list]
            lca_sub[lca_sub == -1] = i
            ancestor_matrix[children_list.reshape(-1,1), children_list] = lca_sub

        # Check symmetry
        assert (ancestor_matrix.T == ancestor_matrix).all()
        assert (-1 == ancestor_matrix).sum() == 0, 'The ontology may have more than one root'

        return ancestor_matrix

    @classmethod
    def _make_tree_igraph(self,
                          graph=None,
                          method='priority',
                          edge_name='smallest_parent',
                          parent_priority=None, edge_priority=None, default_priority=None, optim='max'):
        """Returns copy of graph with new edge attribute marking spanning
        tree"""

        if graph is None:
            graph = self.to_igraph()

        if method=='priority':
            assert 1 == (parent_priority is not None) + (edge_priority is not None)            
            if edge_priority is not None: assert default_priority is not None

            if optim=='min': optim=min
            if optim=='max': optim=max

            graph.es[edge_name] = False

            for v in graph.vs:
                parents = graph.neighbors(v.index, mode='out')

                if len(parents) > 0:

                    """Choose the parent with the highest valued priority"""
                    if parent_priority is not None:
                        small_parent = optim(parents, key=lambda p: parent_priority[p])
                    elif edge_priority is not None:
                        small_parent = optim(parents, key=lambda p: edge_priority.get(graph.get_eid(v.index, p), default_priority))

                    graph.es[graph.get_eid(v.index, small_parent)][edge_name] = True

        else:
            raise Exception('Method not supported')

        return graph

    def to_pickle(self, file, compression='infer'):
        pandas.io.pickle.to_pickle(self, file, compression=compression)

    @classmethod
    def read_pickle(cls, file, compression='infer'):
        return pandas.io.pickle.read_pickle(file, compression=compression)

    # def transitive_closure(self):

        # """Computes the transitive closure on (child term, parent term)
        # relations. Transitivity rules are defined at

        # http://www.geneontology.org/page/ontology-relations

        # """
        
        # relations = ['is_a', 'regulates', 'positively_regulates', 'negatively_regulates', 'has_part', 'part_of', 'gene']
        # relations_index = make_index(relations)
        # go_reasoning = dict()
        # for r in relations:
        #     go_reasoning[('is_a', r)] = r
        #     go_reasoning[(r, 'is_a')] = r
        #     if r != 'has_part':
        #         go_reasoning[(r, 'part_of')] = 'part_of'
        # for r in ['regulates', 'positively_regulates', 'negatively_regulates']:
        #     go_reasoning[(r, 'part_of')] = 'regulates'
        # go_reasoning[('has_part', 'has_part')] = 'has_part'

        # g = self.to_igraph()

        # # Sort from leaves up
        # for i in g.topological_sorting(mode='out'):
        
        #     # Get parents
        #     for j in g.neighbors(i, mode='out'):
                
        #         # Infer GO relations of new edges
        #         base_relation = g.es[g.get_eid(i, j)]['relation']

        #         # Iterate over grandparents
        #         for p in g.neighbors(j, mode='out'):
        #             r = go_reasoning.get((base_relation, g.es[g.get_eid(j, p)]['relation']), None)

        #             # If a relation can't be inferred, then don't add a new edge
        #             if r is not None:
        #                 if -1 != g.get_eid(i, p, error=False):
        #                     # Edge already exists, so take the higher-ranked relation
        #                     e = g.es[g.get_eid(i, p)]
        #                     e['relation'] = r if relations_index[r] > relations_index[e['relation']] else e['relation']
        #                 else:
        #                     # Add new edge with relation
        #                     g.add_edge(i, p, relation=r)                                        

        # return g
        
        # ## Update parent_2_child , child_2_parent, child_2_parent_indices
        
