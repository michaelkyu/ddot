from __future__ import absolute_import, print_function, division

import itertools, multiprocessing, logging, os, collections, random, math, sys, time
from itertools import groupby, combinations
from operator import *
from collections import Counter
import tempfile
from subprocess import Popen, PIPE, STDOUT
import inspect
import shlex
import shutil
import io
from io import StringIO
import json
import datetime

import numpy as np
import pandas as pd
import pandas.io.pickle
import networkx as nx
import igraph
import scipy, scipy.sparse
from scipy.sparse import csr_matrix, coo_matrix
from scipy.stats import hypergeom

import ndex.client as nc
from ndex.networkn import NdexGraph
import ndex.beta.layouts as layouts

import ddot
import ddot.config
from ddot.utils import time_print, set_node_attributes_from_pandas, set_edge_attributes_from_pandas, nx_to_NdexGraph, NdexGraph_to_nx, parse_ndex_uuid, parse_ndex_server, make_index, update_nx_with_alignment, bubble_layout_nx, split_indices_chunk, invert_dict, make_network_public, nx_edges_to_pandas, nx_nodes_to_pandas, ig_edges_to_pandas, ig_nodes_to_pandas, melt_square, nx_set_tree_edges, gridify

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
        assert isinstance(v, str)
        v = g.vs.find(name_eq=v).index

    try:
        g.vs[v]
    except:
        raise Exception("Can't find vertex %s in graph. Consider setting use_v_name=True" % v)
    
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
        g.es['collapsed_length'] = 0
        g.es['collapsed_terms'] = [[] for x in g.es]
        
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
                      output=None,
                      verbose=False):
    
    if output is None:
        with tempfile.NamedTemporaryFile('w', delete=True) as output_file:
            return align_hierarchies(hier1, hier2, iterations, threads,
                                     update_hier1=update_hier1, update_hier2=update_hier2,
                                     mutual_collapse=mutual_collapse,
                                     output=output_file.name,
                                     calculateFDRs=calculateFDRs,
                                     verbose=verbose)

    common_genes = set(hier1.genes) & set(hier2.genes)

    hier1_orig, hier2_orig = hier1, hier2
        
    if len(common_genes) > 0:
        if mutual_collapse:
            hier1, hier2 = Ontology.mutual_collapse(hier1, hier2, verbose=verbose)
            hier1.clear_node_attr()
            hier1.clear_edge_attr()
            hier2.clear_node_attr()
            hier2.clear_edge_attr()
            hier1.propagate('reverse', inplace=True)
            hier2.propagate('reverse', inplace=True)

        def to_file(hier):
            if isinstance(hier, Ontology):
                with tempfile.NamedTemporaryFile('w', delete=False) as f:
                    hier.to_table(f, clixo_format=True)
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

        if threads is None:
            import multiprocessing
            threads = multiprocessing.cpu_count()

        output_dir = tempfile.mkdtemp(prefix='tmp')
        cmd = '{5} {0} {1} 0.05 criss_cross {2} {3} {4} gene'.format(
                  hier1, hier2, output_dir, iterations, threads, calculateFDRs)
        print('Alignment command:', cmd)

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
    else:    
        align1 = pd.DataFrame(columns=['Term', 'Similarity', 'FDR'])
    
    align2 = align1.copy()
    align2.index, align2['Term'] = align2['Term'].values.copy(), align2.index.values.copy()

    append_prefix = lambda x: 'Aligned_%s' % x

    if update_hier1:
        if hasattr(update_hier1, '__iter__'):
            node_attr = hier2_orig.node_attr[update_hier1]
        else:
            node_attr = hier2_orig.node_attr

        hier2_import = pd.merge(pd.DataFrame(index=align2.index), node_attr, left_index=True, right_index=True, how='left')
        assert (hier2_import.index == align2.index).all()
        # Change index to terms in hier1
        hier2_import.index = align2['Term'].copy()
        hier2_import.rename(columns=append_prefix, inplace=True)

    if update_hier2:
        if hasattr(update_hier2, '__iter__'):
            node_attr = hier1_orig.node_attr[update_hier2]
        else:
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
    for line in io.open(obo).read().splitlines():

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
            if curr_term in alt_id:  alt_id[curr_term].append(alt_term)
            else:                          alt_id[curr_term] = [alt_term]
            id2name[alt_term] = name
        elif 'name:' in line:
            name = line.split('name:')[1].strip()
            assert not curr_term in id2name
            id2name[curr_term] = name
        elif 'is_a:' in line:
            parent = line.split('is_a:')[1].strip()
            stanza.append((parent, curr_term, 'is_a'))
        elif 'relationship:' in line:
            line = line.split('relationship:')[1].strip().split()
            if len(line)!=2: print(line)
            assert len(line)==2
            relation, parent = line
            stanza.append((parent, curr_term, relation))
        elif 'namespace:' == line[:10]:
            namespace = line.split('namespace:')[1].strip()
            assert not curr_term in id2namespace
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

    # Remove annotations that have a NOT qualifier
    df = df.loc[df['Qualifier']!='NOT', :]

#    return df.loc[:, ['DB Object ID', 'GO ID']].values.tolist()
    return df

class Ontology(object):
    """A Python representation for constructing, analyzing, and
    manipulating the hierarchical structure of ontologies.

    An Ontology object contains the following attributes for
    representing the hierarchical structure. Do not directly modify
    these attributes.

    Parameters
    ----------
    genes : list
        Names of genes
    
    terms : list
        Names of terms

    gene_2_term : dict
 
        gene_2_term[<gene>] --> list of terms connected to
        <gene>. Terms are represented as their 0-based index in
        self.terms.

    term_2_gene : dict

        term_2_gene[<term>] --> list of genes connected to
        <term>. Genes are represented as their 0-based index in
        self.genes.

    child_2_parent : dict

        child_2_parent[<child>] --> list of the parent terms of <child>

    parent_2_child : dict

        parent_2_child[<parent>] --> list of the children terms of <parent>

    term_sizes : list

        A list of every term's size, i.e. the number of unique genes
        that it and its descendant terms contain. This list has the
        same order as self.terms. It holds that for every i,

        `term_sizes[i] = len(self.term_2_gene[self.terms[i]])`

    """

    NODETYPE_ATTR = 'NodeType'
    GENE_NODETYPE = 'Gene'
    TERM_NODETYPE = 'Term'

    EDGETYPE_ATTR = 'EdgeType'
    GENE_TERM_EDGETYPE = 'Gene-Term'
    CHILD_PARENT_EDGETYPE = 'Child-Parent'

    def __init__(self,
                 hierarchy,
                 mapping,
                 edge_attr=None,
                 node_attr=None,
                 parent_child=False,
                 add_root_name=None,
                 propagate=None,
                 ignore_orphan_terms=False,                 
                 verbose=True,
                 **kwargs):
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
            Ontology.propagate(). If None, then don't
            propagate annotations.

        add_root_name : bool

            The name of an artificial root. If there are multiple
            roots in the ontology, then they are joined into one root
            with this name. Default: Don't create this root.

        ignore_orphan_terms : bool

        """

        if 'empty' in kwargs and kwargs['empty'] is True:
            return
        
        if parent_child:
            hierarchy = [(x[1],x[0]) for x in hierarchy]
            mapping = [(x[1],x[0]) for x in mapping]

        # Cast all node names to strings
        hierarchy = [(str(x[0]),str(x[1])) for x in hierarchy]
        mapping = [(str(x[0]),str(x[1])) for x in mapping]
            
        ## Read term-to-term edges        
        # parent_2_child[<term_name>] --> list of <term_name>'s children terms
        self.parent_2_child = {r: [p[0] for p in q] for r, q in \
                           itertools.groupby(sorted(hierarchy,
                                                    key=lambda a:a[1]),
                                             key=lambda a:a[1])}

        ## Read gene-to-term edges
        # self.gene_2_term[<gene_name>] --> list of terms that <gene_name> is mapped to
        self.gene_2_term = {key: set([a[1] for a in group]) for key, group in \
                            itertools.groupby(sorted(mapping,
                                                     key=lambda a:a[0]),
                                              key=lambda a:a[0])}

        ## Check that the set of terms is the same according to
        ## parent_2_child and self.gene_2_term
        terms_A = set.union(set(self.parent_2_child.keys()),
                            *[set(x) for x in self.parent_2_child.values()])
        if len(self.gene_2_term) > 0:
            terms_B = set.union(*self.gene_2_term.values())
        else:
            terms_B = set([])
            
        if verbose and ignore_orphan_terms and len(terms_B - terms_A)>0:
            print('WARNING: Ignoring {} terms are connected to genes but not to other terms'.format(len(terms_B - terms_A)))
        # if verbose and len(terms_A - terms_B)>0:
        #     print 'WARNING: {} terms connected to other terms but not to genes'.format(len(terms_A - terms_B))

        if ignore_orphan_terms:
            self.terms = sorted(terms_A)
        else:
            self.terms = sorted(terms_A | terms_B)
        self.genes = sorted(self.gene_2_term.keys())

        if add_root_name is not None:
            root_list = self.get_roots()
            if len(root_list) > 1:
                print('Unifying %s roots into one super-root' % len(root_list))
                self.parent_2_child[add_root_name] = root_list
                self.terms.append(add_root_name)
        
        ## terms_index[<term_name>] --> index in self.terms
        self.terms_index = make_index(self.terms)

        ## self.genes_index[<gene_name>] --> index in self.genes        
        self.genes_index = make_index(self.genes)

        ## Convert self.gene_2_term to list term indices rather than term names
        for k, v in self.gene_2_term.items():            
            self.gene_2_term[k] = [self.terms_index[x] for x in self.gene_2_term[k] if x in self.terms_index]

        if node_attr is None:
            self.clear_node_attr()
        else:        
            assert node_attr.index.nlevels == 1
            if node_attr.index.name != 'Node':
                # if verbose:
                #     print("Changing node_attr index name from %s to 'Node'" % node_attr.index.name)
                #     # import traceback
                #     # print traceback.print_stack()
                    
                node_attr.index.name = 'Node'                
            self.node_attr = node_attr
            
        if edge_attr is None:
            self.clear_edge_attr()
        else:
            assert edge_attr.index.nlevels == 2
            edge_attr.index.names = ['Child', 'Parent']
            # if 'Child' in edge_attr.index.names and 'Parent' in edge_attr.index.names:
            #     edge_attr.index = edge_attr.index[['Child', 'Parent']]
            # else:
            #     edge_attr.index.names = ['Child', 'Parent']
            # if edge_attr.index.names != ['Child', 'Parent']:
            #     if verbose:
            #         print("Changing edge_attr index names from %s to ['Child', 'Parent']" % edge_attr.index.names)
            #     edge_attr.index.names = ['Child', 'Parent']
            self.edge_attr = edge_attr

        self._update_fields()

        if propagate:
            self.propagate(direction=propagate, inplace=True)
            self._update_fields()

        self._check_valid()

        # ## Not necessary and requires extra start-up time (perhaps set as a __init__ parameter to precalculate many things)
        
        # empty_terms = sum([x==0 for x in self.term_sizes])
        # if verbose and empty_terms > 0:
        #     print 'WARNING: {} terms are connected to other terms but not to genes'.format(empty_terms), [t for t, x in zip(self.terms, self.term_sizes) if x==0][:5]
        #     # import traceback
        #     # print traceback.print_stack()
                
    def _update_fields(self, reset_term_sizes=True):
        self.child_2_parent = self._get_child_2_parent()
        self.term_2_gene = self._get_term_2_gene()
        if reset_term_sizes:
            self._term_sizes = None

        for t in self.terms:
            if t not in self.parent_2_child:
                self.parent_2_child[t] = []
            if t not in self.child_2_parent:
                self.child_2_parent[t] = []
        
    def add_root(self, root_name, inplace=False):
        """Check if there is a single unifying root term of the ontology. If
        not, then identify the multiple roots and join them under an
        artificial root."""

        if inplace:
            ont = self
        else:
            ont = self.copy()
            
        assert root_name not in ont.terms        
        root_list = ont.get_roots()
        if len(root_list) > 1:
            print('Unifying %s roots into one super-root' % len(root_list))
            ont.parent_2_child[root_name] = root_list
            
        ont.terms.append(root_name)
        ont.terms_index = make_index(sorted(ont.terms))

        for g, t_list in ont.gene_2_term.items():
            ont.gene_2_term[g] = [ont.terms_index[ont.terms[t]] for t in t_list]

        ont.terms.sort()
        ont._update_fields()
        return ont

    def _get_child_2_parent(self):
        """
        Converts self.parent_2_child to child_2_parent

        # child_2_parent[<term_name>] --> list of <term_name>'s parent term names
        """
    
        cp_pairs = []
        for p, c_list in self.parent_2_child.items():
            for c in c_list:
                cp_pairs.append((c,p))
        first = lambda a: a[0]
        cp_pairs.sort(key=first)

        child_2_parent = {
            r: [p[1] for p in q] for r, q in 
            itertools.groupby(cp_pairs, key=first)
        }

        for t in self.terms:
            if t not in child_2_parent:
                child_2_parent[t] = []

        return child_2_parent

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
        
        # Update index to the union of current and new node_attr
        self.node_attr = self.node_attr.reindex(self.node_attr.index.union(node_attr.index))
        
        # Update columns
        for col in node_attr.columns:
            self.node_attr.loc[node_attr.index, col] = node_attr[col]

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

        tmp = set(self.terms) - set([y for x in self.parent_2_child.values() for y in x])
        return sorted(tmp)

    def align(self,
              hier,
              iterations=100,
              threads=None,
              update_self=False,
              update_ref=False,
              align_label=None,
              calculateFDRs=None,
              mutual_collapse=True,
              output=None,
              verbose=False):
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
            parallelize the the null model randomizations. Default:
            The number of CPU cores returned by
            multiprocessing.cpu_count()

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
        
        alignment = align_hierarchies(
            self,
            hier,
            iterations,
            threads,
            update_hier1=update_self,
            update_hier2=update_ref,
            calculateFDRs=calculateFDRs,
            mutual_collapse=mutual_collapse,
            output=output,
            verbose=verbose)

        # Set labels based on ontology alignment
        if align_label and (('Aligned_%s' % align_label) in self.node_attr.columns):
            label_attr = self.node_attr['Aligned_%s' % align_label]            
            label_attr = {k : '%s\n%s' % (k,v) for k, v in label_attr.iteritems() if not pd.isnull(v)}
            label_attr = pd.Series(label_attr, name='Label').to_frame()
            self.update_node_attr(label_attr)

        return alignment
    
    def _make_dummy(self, tree_edges=None):
        """For each term T in the ontology, create a new dummy term that
        indirectly connect T's to T. For example, if g1 and g2 are in
        T, then a new term dummy_T is created so that the new ontology
        consists of

        g1 --> T_dummy
        g2 --> T_dummy
        T_dummy --> T

        Parameters
        ----------
        tree_edges : list

            List of (child, parent) edges that constitute a spanning
            tree of the ontology. If specified, then for each term T,
            only the genes that are connected to T in the spanning
            tree will be re-routed to the dummy node.

            Default: None. This restriction will not apply

        Returns
        -------
        : ddot.Ontology.Ontology

        """
        
        ont = self
        
        new_gene_2_term = []
        new_child_2_parent = []
        for t in ont.terms:
            
            used_dummy = False
            if len(ont.parent_2_child[t]) > 0:
                dummy_term = 'dummy2_%s' % t
            else:
                dummy_term = t
            for g in [ont.genes[g] for g in ont.term_2_gene[t]]:
                if (tree_edges is None) or (g,t) in tree_edges:
                    new_gene_2_term.append((g, dummy_term))
                    used_dummy=True
                    
            if used_dummy and dummy_term != t:
                new_child_2_parent.append([dummy_term, t])
                
            for p in ont.child_2_parent[t]:
                if (tree_edges is None) or (t,p) in tree_edges:
                    new_child_2_parent.append((t, p))

        ont_dummy = Ontology(new_child_2_parent, new_gene_2_term)
        return ont_dummy

    def _collect_transform(self,
                           tree_edges=None,
                           hidden_gene=True,
                           hidden_parent=True,
                           hidden_child=True):
        """
        Creates intermediate duplicate nodes 
        """
        
        ont = self
                
        if tree_edges is None:
            tree_edges = self.get_tree()

        nodes_copy = {v : 1 for v in ont.genes + ont.terms}
        def get_copy(u):
            u_name = '%s.%s' % (u, nodes_copy[u])
            nodes_copy[u] += 1            
            return u_name

        collect_nodes = []
        new_gene_2_term = []
        new_child_2_parent = []
        for t in ont.terms:
            
            ## Gene-term connections
            collect_hidden_gene = 'collect_hidden_gene_%s' % t
            used_hidden_gene = False
            for g in [ont.genes[g] for g in ont.term_2_gene[t]]:
                if (not hidden_gene) or ((g, t) in tree_edges):
                    new_gene_2_term.append((g, collect_hidden_gene))
                    used_hidden_gene = True
                else:
                    new_gene_2_term.append((get_copy(g), collect_hidden_gene))
                    used_hidden_gene = True
            if used_hidden_gene:
                collect_nodes.append(collect_hidden_gene)
                new_child_2_parent.append((collect_hidden_gene, t))

            ## Parent-child term connections            
            collect_hidden_child = 'collect_hidden_child_%s' % t
            collect_hidden_parent = 'collect_hidden_parent_%s' % t
            used_hidden_child, used_hidden_parent = False, False 
            for c in ont.parent_2_child[t]:
                if (not hidden_child) or ((c,t) in tree_edges):
                    new_child_2_parent.append((c,t))
                else:
                    new_child_2_parent.append((get_copy(c), collect_hidden_child))
                    used_hidden_child = True
            for p in ont.child_2_parent[t]:
                if hidden_parent and ((t,p) not in tree_edges):
                    new_child_2_parent.append((get_copy(p), collect_hidden_parent))
                    used_hidden_parent = True                
            if used_hidden_child:
                collect_nodes.append(collect_hidden_child)
                new_child_2_parent.append((collect_hidden_child, t))        
            if used_hidden_parent:
                collect_nodes.append(collect_hidden_parent)
                new_child_2_parent.append((collect_hidden_parent, t))
        
        ont_collect = Ontology(new_child_2_parent,
                               new_gene_2_term,
                               node_attr=ont.node_attr.copy(),
                               edge_attr=ont.edge_attr.copy(),
                               verbose=False)    
                
        ##################################################
        # Set Original_Name and Size for Duplicate Nodes #

        new_and_orig = [('%s.%s' %(v,i), v) for v, copy_num in nodes_copy.items()
                        for i in (range(1, copy_num) if copy_num>1 else [])]
        new_2_orig = dict(new_and_orig)

        df = pd.DataFrame({'orig_tmp' : [x[1] for x in new_and_orig],
                           'Hidden' : True},
                          index=[x[0] for x in new_and_orig])
        df = df.astype({'orig_tmp' : np.str, 'Hidden' : np.bool})

        # For duplicate nodes, set the Original_Name attribute to the name of the original node
        merge = pd.merge(df, ont.node_attr, how='left', left_on=['orig_tmp'], right_index=True)
        if 'Original_Name' in merge:
            unset = pd.isnull(merge['Original_Name']).values
            merge.loc[unset, 'Original_Name'] = df.loc[unset, 'orig_tmp'].values
        else:
            merge['Original_Name'] = df['orig_tmp'].values
        del merge['orig_tmp']

        # Set the 'Size' attribute of duplicate nodes to be the 'Size'
        # of the original node. If the original node is a term with no
        # 'Size' attribute, then set 'Size' to be the number of genes
        # in the term
        in_merge = set(merge.index)
        for node in merge.index:
            if node in new_2_orig:
                orig = new_2_orig[node]
                if orig in in_merge and not pd.isnull(merge.loc[orig, 'Size']):
                    merge.loc[node, 'Size'] = merge.loc[new_2_orig[node], 'Size']
                elif orig in ont.terms_index:
                    merge.loc[node, 'Size'] = ont.term_sizes[ont.terms_index[orig]]
        
        # Append attributes for the new nodes
        try:
            # Used for pandas version >= 0.23
            ont_collect.node_attr = pd.concat([ont.node_attr, merge], axis=0, sort=True)
        except:
            ont_collect.node_attr = pd.concat([ont.node_attr, merge], axis=0)
        
        ########################################
        # Set Label and Size for collect nodes #
        ########################################
        
        def get_label(x):
            if 'hidden_child' in x:
                return 'Linked Children'
            elif 'hidden_parent' in x:
                return 'Linked Parents'
            elif 'hidden_gene' in x:
                return 'Linked Genes'
            elif 'tree_gene' in x:
                return 'Genes'
        
        collect_attr = pd.DataFrame(
            {'Size' : 1,
             'Label' : [get_label(x) for x in collect_nodes],
             'is_collect_node' : True},
             index=collect_nodes)        
        ont_collect.update_node_attr(collect_attr)
                
        return ont_collect

    def unfold(self,
               duplicate=None,
               genes_only=False,
               levels=None,
               tree_edges=None):
        """Traverses the ontology from the root to the leaves while
        duplicating nodes during the traversal to create a tree representation.

        Traverse the ontology from the root nodes to the leaves in a
        breadth-first manner. Each time a node is traversed, then
        create a duplicate of it

        Parameters
        ---------- 

        duplicate : list

            Nodes to duplicate for unfolding. Default: all genes and terms

        genes_only : bool

            If True, then duplicate all of the genes and none of the terms. Default: False

        levels :

        """

        ont = self.propagate(direction='reverse', inplace=False)
                
        hidden_mode = levels is not None
        if hidden_mode:            
            if tree_edges is None:
                tree_edges = self.get_tree()
            hidden_depth = {}
        
        if genes_only:
            duplicate = ont.genes
        elif duplicate is None:
            duplicate = ont.genes + ont.terms
        nodes_copy = {x : 0 for x in duplicate}
            
        def get_name(u):
            if u in nodes_copy:
                u_name = '%s.%s' % (u, nodes_copy[u])
                nodes_copy[u] += 1
            else:
                u_name = u
            return u_name

        to_expand = []
        new_2_orig = {}
        for u in ont.get_roots():
            u_name = get_name(u)
            new_2_orig[u_name] = u
            to_expand.append(u_name)

            if hidden_mode:
                hidden_depth[u_name] = 0
        expanded = set(to_expand)

        hierarchy, mapping = [], []

        # Manual bfs
        curr = 0
        while curr < len(to_expand):
            v_name = to_expand[curr]
            v = new_2_orig[v_name]
                
            for u in [ont.genes[u] for u in ont.term_2_gene[v]]:
                u_name = get_name(u)
                new_2_orig[u_name] = u
                mapping.append((u_name, v_name))

                if hidden_mode:
                    v_depth = hidden_depth[v_name]
                    if v_depth==0:
                        if (u,v) in tree_edges:
                            hidden_depth[u_name] = 0
                        else:
                            hidden_depth[u_name] = 1                        
                    elif v_depth < levels:
                        hidden_depth[u_name] = v_depth + 1
                
            for u in ont.parent_2_child[v]:
                u_name = get_name(u)
                new_2_orig[u_name] = u
                hierarchy.append((u_name, v_name))

                if hidden_mode:
                    v_depth = hidden_depth[v_name]
                    insert = u_name not in expanded
                    if v_depth==0 and ((u,v) in tree_edges):
                        hidden_depth[u_name] = 0                        
                    elif v_depth < levels:
                        hidden_depth[u_name] = v_depth + 1
                    else:
                        insert = False
                else:
                    insert = u_name not in expanded
                        
                if insert:
                    to_expand.append(u_name)
                    expanded.add(u_name)

            curr += 1

        new_nodes, orig_nodes = zip(*new_2_orig.items())
        new_nodes, orig_nodes = list(new_nodes), list(orig_nodes)
        ont.node_attr = ont.node_attr.reindex(list(set(orig_nodes)))
        
        node_attr = ont.node_attr.loc[orig_nodes, :].copy()
        if 'Original_Name' in node_attr:
            unset = pd.isnull(node_attr['Original_Name']).values
            node_attr.loc[unset, 'Original_Name'] = np.array(orig_nodes)[unset]
        else:            
            node_attr['Original_Name'] = orig_nodes

        if hidden_mode:
            node_attr['Level'] = [hidden_depth[v] for v in new_nodes]
        
        node_attr.index = new_nodes
        node_attr.dropna(axis=0, how='all', inplace=True)
        
        new_edges = hierarchy + mapping
        old_edges = [(new_2_orig[u], new_2_orig[v]) for u, v in new_edges]
        in_index = [x in ont.edge_attr.index for x in old_edges]
        if sum(in_index) > 0:
            edge_attr = ont.edge_attr.loc[[x for x, y in zip(old_edges, in_index) if y], :].copy()
            edge_attr.index = pd.MultiIndex.from_tuples([x for x, y in zip(new_edges, in_index) if y])
            edge_attr.dropna(axis=0, how='all', inplace=True)
        else:
            edge_attr = None
            
        ont = Ontology(hierarchy,
                       mapping,
                       edge_attr=edge_attr,
                       node_attr=node_attr,
                       parent_child=False,
                       verbose=False)
        return ont
    
    
    def _to_networkx_no_layout(self):
        G = nx.DiGraph()

        #################################
        ### Add nodes and node attributes

        G.add_nodes_from(self.genes + self.terms)

        set_node_attributes_from_pandas(G, self.node_attr)

        # Ensure that all 'Size' values are the same numeric type
        if 'Size' in self.node_attr.columns:
            dtype = self.node_attr['Size'].dtype
            if dtype in [np.dtype('float16'), np.dtype('float32'), np.dtype('float64')]:
                dtype = float
            else:
                dtype = int
        else:
            dtype = int
                
        for t in self.terms:
            G.node[t][self.NODETYPE_ATTR] = self.TERM_NODETYPE
            if ('Size' not in G.node[t]) or pd.isnull(G.node[t]['Size']):
                G.node[t]['Size'] = dtype(self.term_sizes[self.terms_index[t]])
            G.node[t]['isRoot'] = False
        for g in self.genes:
            G.node[g][self.NODETYPE_ATTR] = self.GENE_NODETYPE
            if ('Size' not in G.node[g]) or pd.isnull(G.node[g]['Size']):
                G.node[g]['Size'] = dtype(1)
            G.node[g]['isRoot'] = False

        # Identify the root
        root = self.get_roots()[0]
        G.node[root]['isRoot'] = True

        # Set the node attribute 'Label'. If the node has a "Original
        # Name" attribute, indicating that it is a duplicate, then use
        # that. Otherwise, use the node's name.
        for x in self.genes + self.terms:
            data = G.node[x]            
            if ('Label' not in data) or pd.isnull(data['Label']):
                if ('Original_Name' in data) and (not pd.isnull(data['Original_Name'])):
                    data['Label'] = data['Original_Name']
                else:
                    data['Label'] = x

        #################################
        ### Add edges and edge attributes

        G.add_edges_from([(g, self.terms[t],
                           {self.EDGETYPE_ATTR : self.GENE_TERM_EDGETYPE}) \
                          for g in self.genes for t in self.gene_2_term[g]])

        G.add_edges_from([(c, p,  
                           {self.EDGETYPE_ATTR : self.CHILD_PARENT_EDGETYPE}) \
                          for p in self.terms for c in self.parent_2_child.get(p, [])])

        set_edge_attributes_from_pandas(G, self.edge_attr)

        return G

    def expand(self, spanning_tree=True):
        if spanning_tree is True:
            ont = self._collect_transform()
        else:
            # Assume a list of tree edges are supplied
            ont = self._collect_transform(spanning_tree)

        G_tree = ont.get_tree(ret='ontology')._to_networkx_no_layout()
        pos = bubble_layout_nx(G_tree)
        tmp = ont.node_attr[['Label', 'is_collect_node']].dropna()
        collect_nodes = tmp[tmp['is_collect_node']].index
        gridify(collect_nodes, pos, G_tree)

        ## Remove collector nodes               

        def decide_delete(v):
            return ((not layout_params['hidden_parent'] and v=='Linked Parents') or
                    (not layout_params['hidden_child'] and v=='Linked Children') or
                    (not layout_params['hidden_gene'] and v=='Linked Genes'))
        tmp = ont.node_attr[['Label', 'is_collect_node']].dropna()
        tmp = tmp[tmp['is_collect_node']]
        tmp = tmp[tmp['Label'].apply(decide_delete)]
        to_delete = tmp.index.tolist()

        ont_red = ont
        if len(to_delete) > 0:                    
            # Need fast special delete
            ont_red = ont_red.delete(to_delete=to_delete, preserve_transitivity=True)

        # Set the original term sizes for the original copy of
        # each term (not the duplicates)
        ont_red.update_node_attr(pd.DataFrame({'Size' : self.term_sizes}, index=self.terms))
        G = ont_red._to_networkx_no_layout()

        nodes_set = set(G.nodes())
        G.pos = {n : (float(scale*p[0]), float(scale*p[1])) for n, p in pos.items() if n in nodes_set}
        nx_set_tree_edges(G, ont_red.get_tree())

        ######################################################
        # TODO: move this visual styling outside of the layout
        # functionality

        nx.set_edge_attributes(G, values='ARROW', name='Vis:EDGE_SOURCE_ARROW_SHAPE')
        nx.set_edge_attributes(G, values='NONE', name='Vis:EDGE_TARGET_ARROW_SHAPE')

        for v, data in G.nodes(data=True):
            # if 'collect_hidden' in v and 'is_collect_node' in data and data['is_collect_node']:
            #     for u in G.predecessors(v):
            #         G.node[u]['Vis:Fill Color'] = '#3182BD'
            try:
                if 'collect_hidden_parent' in v and 'is_collect_node' in data and data['is_collect_node']:
                    for _, _, data in G.in_edges(v, data=True):
                        data["Vis:EDGE_TARGET_ARROW_SHAPE"] = 'ARROW'
                        data["Vis:EDGE_SOURCE_ARROW_SHAPE"] = 'NONE'
            except:
                print(data)
                print('v', v)
                print('collect_hidden_parent' in v)
                print('is_collect_node' in data)
                print(data['is_collect_node'])
                raise
        
    def to_networkx(self,
                    layout='bubble',
                    spanning_tree=True,
                    layout_params=None,
                    verbose=False):

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

        Returns
        -------

        : nx.DiGraph

        """

        default_layout_params = {'hidden_parent' : True,
                                 'hidden_child' : False,
                                 'hidden_gene' : False}
        if layout_params is not None:
            default_layout_params.update(layout_params)
        layout_params = default_layout_params
            
        if spanning_tree:
            scale = 1

            if layout is None or layout=='bubble':
                G = self._to_networkx_no_layout()
                if spanning_tree is True:
                    tree_edges = self.get_tree()
                else:
                    tree_edges = spanning_tree
                nx_set_tree_edges(G, tree_edges)
  
                if layout=='bubble':
                    G_tree = self.propagate('reverse')._make_dummy(tree_edges)._to_networkx_no_layout()
                    pos = bubble_layout_nx(G_tree, verbose=verbose)
                    gridify([v for v in G_tree.nodes() if 'dummy2' in v], pos, G_tree)
                    G.pos = {n : (float(scale*p[0]), float(scale*p[1])) for n, p in pos.items() if 'dummy2' not in n}

            elif layout=='bubble-collect':

                if spanning_tree is True:
                    ont = self._collect_transform()
                else:
                    # Assume a list of tree edges are supplied
                    ont = self._collect_transform(spanning_tree)
                    
                G_tree = ont.get_tree(ret='ontology')._to_networkx_no_layout()
                pos = bubble_layout_nx(G_tree)
                tmp = ont.node_attr[['Label', 'is_collect_node']].dropna()
                collect_nodes = tmp[tmp['is_collect_node']].index
                gridify(collect_nodes, pos, G_tree)

                ## Remove collector nodes               

                def decide_delete(v):
                    return ((not layout_params['hidden_parent'] and v=='Linked Parents') or
                            (not layout_params['hidden_child'] and v=='Linked Children') or
                            (not layout_params['hidden_gene'] and v=='Linked Genes'))
                tmp = ont.node_attr[['Label', 'is_collect_node']].dropna()
                tmp = tmp[tmp['is_collect_node']]
                tmp = tmp[tmp['Label'].apply(decide_delete)]
                to_delete = tmp.index.tolist()
                
                ont_red = ont
                if len(to_delete) > 0:                    
                    # Need fast special delete
                    ont_red = ont_red.delete(to_delete=to_delete, preserve_transitivity=True)
                    
                # Set the original term sizes for the original copy of
                # each term (not the duplicates)
                ont_red.update_node_attr(pd.DataFrame({'Size' : self.term_sizes}, index=self.terms))
                G = ont_red._to_networkx_no_layout()

                nodes_set = set(G.nodes())
                G.pos = {n : (float(scale*p[0]), float(scale*p[1])) for n, p in pos.items() if n in nodes_set}
                nx_set_tree_edges(G, ont_red.get_tree())

                ######################################################
                # TODO: move this visual styling outside of the layout
                # functionality
                
                nx.set_edge_attributes(G, values='ARROW', name='Vis:EDGE_SOURCE_ARROW_SHAPE')
                nx.set_edge_attributes(G, values='NONE', name='Vis:EDGE_TARGET_ARROW_SHAPE')

                for v, data in G.nodes(data=True):
                    # if 'collect_hidden' in v and 'is_collect_node' in data and data['is_collect_node']:
                    #     for u in G.predecessors(v):
                    #         G.node[u]['Vis:Fill Color'] = '#3182BD'
                    try:
                        if 'collect_hidden_parent' in v and 'is_collect_node' in data and data['is_collect_node']:
                            for _, _, data in G.in_edges(v, data=True):
                                data["Vis:EDGE_TARGET_ARROW_SHAPE"] = 'ARROW'
                                data["Vis:EDGE_SOURCE_ARROW_SHAPE"] = 'NONE'
                    except:
                        print(data)
                        print('v', v)
                        print('collect_hidden_parent' in v)
                        print('is_collect_node' in data)
                        print(data['is_collect_node'])
                        raise
            else:
                raise Exception('Unsupported layout: %s', layout)

            if layout is not None:
                nx.set_node_attributes(G, values={n : x for n, (x,y) in G.pos.items()}, name='x_pos')
                nx.set_node_attributes(G, values={n : y for n, (x,y) in G.pos.items()}, name='y_pos')
                
        else:
            G = self._to_networkx_no_layout()
            
        return G

    @classmethod
    def from_table(cls,
                   table,
                   parent=0,
                   child=1,
                   is_mapping=None,
                   mapping=None,
                   mapping_parent=0,
                   mapping_child=1,
                   header=0,
                   propagate=False,
                   verbose=False,
                   clixo_format=False,
                   clear_default_attr=True,
                   **kwargs):
        """Create Ontology from a tab-delimited table or pandas DataFrame.

        Duplicate gene-term or term-term connections in the table are removed.

        Parameters
        ----------

        table : pandas.DataFrame, file-like object, or filename

            A table that lists (child term, parent term) pairs. If
            mapping==None, then this table should also include (gene,
            term) pairs.

        parent : int or str

            Column for parent terms in table (index or name of column)
        
        child : int or str

            Column for child terms and genes in table (index or name of column)

        is_mapping : function

            A function that is applied on each row and returns True if
            the row represents a (gene, term) pair and False
            otherwise. This function is only applied when a separate
            table of (gene, term) pairs is not specified,
            i.e. mapping==None.

            The default function is `lambda row: row[2]=={0}`

            which tests if the third column equals the string "{0}".
            
        mapping : pandas.DataFrame, file-like object, or filename (optional)

            A separate table listing only (gene, term) pairs

        mapping_parent : int or str

            Column for terms in mapping table (index or name of column)
        
        mappping_child : int or str

            Column for genes in mapping table (index or name of column)

        header : int or None

            Row number to use as the column names, which are then
            stored in the resulting Ontology object's `edge_attr`
            field. For example if `header=0` (default), then the first
            row is assumed to be column names. If `header=None`, then
            no column names are assumed.

        propagate : None or str

            The direction ('forward' or 'reverse') for propagating
            gene-term annotations up the hierarchy with
            Ontology.propagate(). If None, then don't
            propagate annotations.

        clixo_format : bool

            If True, The table is assumed to be in the same format
            produced by the CLIXO C++ implementation. In particular,
            table has three columns:

            Column 1) Parent Term
            Column 2) Child Term or Gene
            Column 3) The string "gene" if the row is a
                      gene-term mapping, otherwise the string "default".

            The table is also assumed to have no column headers (i.e. header=False)

        clear_default_attr: bool

            If True (default), then remove the edge attribute
            'EdgeType' created using Ontology.to_table(). This
            attribute was created to make the table be an equivalent
            representation of an Ontology object; however, it is no
            longer necessary after reconstructing the Ontology object.
        
        Returns
        -------

        : ddot.Ontology.Ontology

        """.format(cls.GENE_TERM_EDGETYPE)

        if clixo_format:
            ont = cls.from_table(
                table,
                parent=0,
                child=1,
                is_mapping=lambda x: x[2]=='gene',
                header=None,
                clixo_format=False,
                verbose=verbose)
            ont.edge_attr.columns = map(str, ont.edge_attr.columns)
            del ont.edge_attr['2']

            return ont
        
        if is_mapping is None:
            if mapping is None:
            #     print('WARNING: no gene-term connections '
            #           'were specified by the is_mapping '
            #           'function or separate table. '
            #           'Default: assume a gene-term connection when the 3rd column equals %s' % cls.GENE_TERM_EDGETYPE)
                is_mapping = lambda x: x.iloc[2]==cls.GENE_TERM_EDGETYPE
                
        # Read table
        try:
            table = pd.read_table(table, comment='#', header=header)
        except:
            assert isinstance(table, pd.DataFrame)

        if child not in table.columns:
            child = table.columns[child]
        if parent not in table.columns:
            parent = table.columns[parent]
            
        for col in [child, parent]:
            table.loc[:,col] = table.loc[:,col].astype(str)

        edge_attr = table.set_index([child, parent])
        edge_attr.index.rename(['Child', 'Parent'], inplace=True)

        if mapping is None:
            # Extract gene-term connections from table
            mask = table.apply(is_mapping, axis=1)
            mapping = table.loc[mask, :].loc[:,[child, parent]]
            hierarchy = table.loc[~mask, :].loc[:,[child, parent]]
            mapping_child, mapping_parent = child, parent
        else:
            # Read separate table of gene-term connections
            try:
                mapping = pd.read_table(mapping, comment='#', header=header)
            except:
                assert isinstance(mapping, pd.DataFrame)
                
            if mapping_child not in mapping.columns:
                mapping_child = mapping.columns[mapping_child]
            if mapping_parent not in mapping.columns:
                mapping_parent = mapping.columns[mapping_parent]
        
            for col in [mapping_child, mapping_parent]:
                mapping.loc[:,col] = mapping.loc[:,col].astype(str)
                
            mapping_attr = mapping.set_index([mapping_child, mapping_parent])
            mapping_attr.index.rename(['Child', 'Parent'], inplace=True)

            try:
                # Used for pandas version >= 0.23
                edge_attr = pd.concat([edge_attr, mapping_attr], sort=True)
            except:
                edge_attr = pd.concat([edge_attr, mapping_attr])
            mapping = mapping.loc[:,[mapping_child, mapping_parent]]
            hierarchy = table.loc[:,[child, parent]]

        dups = mapping.duplicated([mapping_child, mapping_parent]).sum()
        if dups > 0:
            print('WARNING: Dropping %s duplicate gene-term connections' % dups)
            mapping.drop_duplicates([mapping_child, mapping_parent], inplace=True)

        dups = hierarchy.duplicated([child, parent]).sum()
        if dups > 0:
            print('WARNING: Dropping %s duplicate term-term connections' % dups)
            hierarchy.drop_duplicates([child, parent], inplace=True)

        edge_attr = edge_attr.loc[~ edge_attr.index.duplicated(), :]
        edge_attr.index.names = ['Child', 'Parent']

        if clear_default_attr:
            if cls.EDGETYPE_ATTR in edge_attr:
                del edge_attr[cls.EDGETYPE_ATTR]
            
        mapping, hierarchy = mapping.values.tolist(), hierarchy.values.tolist()        
        
        return cls(hierarchy,
                   mapping,
                   parent_child=False,
                   edge_attr=edge_attr,
                   propagate=propagate,
                   verbose=verbose,
                   **kwargs)

    @classmethod
    def from_scipy_linkage(cls, Z):
        """Creates an Ontology object from a linkage matrix created by scipy's
        hierarchical/agglomerative clustering. Note that this form of
        clustering produces a binary tree.

        """
        import scipy.cluster.hierarchy
                
        rootnode, nodelist = scipy.cluster.hierarchy.to_tree(Z, rd=True)
        leaves = set(scipy.cluster.hierarchy.leaves_list(Z))
        hierarchy, mapping = [], []
        for v in nodelist:
            v_id = v.get_id()
            if v.get_left():
                child = v.get_left().get_id()
                if child in leaves:
                    mapping.append((v_id, child))
                else:
                    hierarchy.append((v_id, child))
            if v.get_right():
                child = v.get_right().get_id()
                if child in leaves:
                    mapping.append((v_id, child))
                else:
                    hierarchy.append((v_id, child))
        return cls(hierarchy, mapping, parent_child=True)
    
    @classmethod
    def from_ndex(cls,
                  ndex_uuid,
                  ndex_user=None,
                  ndex_pass=None,
                  ndex_server=None,
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
            ndex_user = ddot.config.ndex_user
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
                      edgetype_value=None,
                      clear_default_attr=True):
        """Converts a NetworkX object to an Ontology object. Gene and terms
        are distinguished by an edge attribute.

        Parameters
        ----------
        G : nx.DiGraph

        edgetype_attr : str

            Name of the edge attribute that distinguishes a (gene,
            term) pair from a (child term, parent term) pair

        edgetype_value : str
        
            Value of the edge attribute for (gene, term) pairs

        clear_default_attr : bool

            If True (default), then remove the node and edge
            attributes that are created in a NetworkX graph using
            Ontology.to_networkx() or Ontology.to_ndex(). These
            attributes include 'Label', 'Size', 'NodeType', and
            'EdgeType'. These attributes were created to make the
            NetworkX graph be an equivalent representation of an
            Ontology object; however, they are no longer necessary
            after reconstrcting the Ontology object.

        Returns
        -------
        : ddot.Ontology.Ontology

        """
        
        if edgetype_attr is None:
            edgetype_attr=cls.EDGETYPE_ATTR
        if edgetype_value is None:
            edgetype_value=cls.GENE_TERM_EDGETYPE
            
        hierarchy = []
        mapping = []
        for u, v, attr in G.edges(data=True):
            if attr[edgetype_attr] == edgetype_value:
                mapping.append((u, v))
            else:
                hierarchy.append((u, v))

        edge_attr = nx_edges_to_pandas(G)
        node_attr = nx_nodes_to_pandas(G)
        
        ont = cls(hierarchy,
                  mapping,
                  node_attr=node_attr,
                  edge_attr=edge_attr)

        if clear_default_attr:
            for attr in [Ontology.NODETYPE_ATTR, 'Label', 'Size', 'isRoot', 'x_pos', 'y_pos']:
                if attr in ont.node_attr.columns:
                    del ont.node_attr[attr]

            for attr in [edgetype_attr, 'Is_Tree_Edge']:
                if attr in ont.edge_attr.columns:
                    del ont.edge_attr[attr]            

        return ont
    
    @classmethod
    def from_igraph(cls,
                    G,
                    edgetype_attr=None,
                    edgetype_value=None,
                    verbose=False):
        """Converts a igraph Graph object to an Ontology object. Gene and terms
        are distinguished by an edge attribute.

        Parameters
        ----------
        G : igraph.Graph

        edgetype_attr : str

            Name of the edge attribute that distinguishes a (gene,
            term) pair from a (child term, parent term) pair

        edgetype_value : str
        
            Value of the edge attribute for (gene, term) pairs

        Returns
        -------
        : ddot.Ontology.Ontology

        """
                
        if edgetype_attr is None:
            edgetype_attr=cls.EDGETYPE_ATTR
        if edgetype_value is None:
            edgetype_value=cls.GENE_TERM_EDGETYPE
        
        hierarchy = []
        mapping = []
        for e in G.es:
            u = G.vs[e.source]['name']
            v = G.vs[e.target]['name']
            if e[edgetype_attr] == edgetype_value:
                mapping.append((u, v))
            else:
                hierarchy.append((u, v))

        edge_attr = ig_edges_to_pandas(G)
        node_attr = ig_nodes_to_pandas(G)
        edge_attr.index.names = ['Child', 'Parent']
        node_attr.index.name = 'Node'

        ont = cls(hierarchy,
                  mapping,
                  node_attr=node_attr,
                  edge_attr=edge_attr,
                  verbose=verbose)

        for attr in [Ontology.NODETYPE_ATTR]:
            if attr in ont.node_attr.columns:
                del ont.node_attr[attr]
            
        for attr in [edgetype_attr, 'Is_Tree_Edge']:
            if attr in ont.edge_attr.columns:
                del ont.edge_attr[attr]            

        return ont

    def collapse_ontology(self,                          
                          method='mhkramer',
                          to_keep=None,
                          min_term_size=2,
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
            assert to_keep is None, 'to_keep is only supported for method=="python"'
            
            # Propagate forward and then reverse
            ont = self.copy()
            ont = self.propagate(direction='forward', inplace=False)
            ont.propagate(direction='reverse', inplace=True)

            top_level = os.path.dirname(os.path.abspath(inspect.getfile(ddot)))
            collapseRedundantNodes = os.path.join(top_level, 'alignOntology', 'collapseRedundantNodes')
            # assert os.path.isdir(ddot.config.alignOntology)
            # collapseRedundantNodes = os.path.join(ddot.config.alignOntology, 'collapseRedundantNodes')
            assert os.path.isfile(collapseRedundantNodes)

            with tempfile.NamedTemporaryFile('w', delete=False) as f:
                ont.to_table(f, clixo_format=True)
            try:                
                cmd = '%s %s' % (collapseRedundantNodes, f.name)
                print('collapse command:', cmd)
                p = Popen(shlex.split(cmd), shell=False, stdout=PIPE, stderr=PIPE)
                collapsed, err = p.communicate()
                collapsed = collapsed.decode()
            finally:
                os.remove(f.name)
            
            ont = Ontology.from_table(
                StringIO(collapsed),
                is_mapping=lambda x: x[2]=='gene',
                clixo_format=True
            )
            ont.clear_edge_attr()
            ont.update_node_attr(self.node_attr)
            ont.update_edge_attr(self.edge_attr)

            return ont

        elif method=='python':
            ont = self.propagate('forward', inplace=False)
            term_hash = {t : hash(tuple(g_list)) for t, g_list in ont.term_2_gene.items()}
            to_collapse = set()
            for p in ont.parent_2_child:
                for c in ont.parent_2_child[p]:
                    if term_hash[p] == term_hash[c]:
                        to_collapse.add(p)
            
            if min_term_size is not None:
                to_collapse = to_collapse | set([t for t, s in zip(ont.terms, ont.term_sizes) if s < min_term_size])

            if to_keep is not None:
                to_collapse = to_collapse - set(to_keep)

#            print('to_collapse:', sorted(to_collapse))
            
            ont.propagate('reverse', inplace=True)
            ont_red = ont.delete(to_delete=to_collapse, preserve_transitivity=True)

            return ont_red
    
    @classmethod
    def mutual_collapse(cls,
                        ont1,
                        ont2,                        
                        verbose=False):
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

        if len(common_genes) > 0:
            ont1 = ont1.delete(to_delete=set(ont1.genes) - common_genes, inplace=False)
            ont1_collapsed = ont1.collapse_ontology(method='mhkramer')
            ont2 = ont2.delete(to_delete=set(ont2.genes) - common_genes, inplace=False)
            ont2_collapsed = ont2.collapse_ontology(method='mhkramer')
        else:
            raise Exception('No common genes between ontologies')

        if verbose:
            print('ont1_collapsed:', ont1_collapsed.summary())
            print('ont2_collapsed:', ont2_collapsed.summary())

        return ont1_collapsed, ont2_collapsed

    def focus(self,
              branches=None,
              genes=None,
              collapse=False,
              root=True,
              verbose=True):
        """

        """

        assert (branches is not None) or (genes is not None)
        
        to_keep = np.array(self.genes + self.terms)
        if branches is not None:
            to_keep = to_keep[self.connected(to_keep, branches).sum(1) > 0]
            if verbose:
                print('Genes and Terms to keep:', to_keep.size)
        if genes is not None:
            to_keep = to_keep[self.connected(genes, to_keep).sum(0) > 0]
            if verbose:
                print('Genes and Terms to keep:', to_keep.size)

        if root:
            while True:
                common_root = self.common_ancestors(to_keep, minimal=True)
                if common_root in to_keep or len(common_root)<=1:
                    break            
                else:
                    print('Adding', common_root)
                    to_keep = np.append(to_keep, common_root)
        
        ont = self.delete(to_keep=to_keep, preserve_transitivity=True)

        if collapse:
            ont = ont.collapse_ontology(method='python', to_keep=ont.get_roots())
        
        df = ont.to_table(edge_attr=True)

        new_connections = []
        for t in ont.terms:
            removed_genes = set([self.genes[g] for g in self.term_2_gene[t]]) - set([ont.genes[g] for g in ont.term_2_gene[t]])
            removed_terms = set(self.parent_2_child[t]) - set(ont.parent_2_child[t])
            if len(removed_genes) > 0:
                new_connections.append(('%s_%s_other_genes' % (t, len(removed_genes)), t, self.GENE_TERM_EDGETYPE))
            if len(removed_terms) > 0:
                new_connections.append(('%s_%s_other_terms' % (t, len(removed_terms)), t, self.CHILD_PARENT_EDGETYPE))

        if len(new_connections) > 0:
            new_connections = pd.DataFrame(new_connections)
            new_connections.columns = ['Child', 'Parent', self.EDGETYPE_ATTR]
            new_nodes = new_connections['Child'].values.tolist()
            new_connections['Summary'] = True
            df['Summary'] = False
            try:
                # Used for pandas version >= 0.23
                tmp = pd.concat([df, new_connections], ignore_index=True, sort=True)
            except:
                tmp = pd.concat([df, new_connections], ignore_index=True)
            df = tmp[df.columns]

        ont = Ontology.from_table(df)
        ont.update_node_attr(self.node_attr)
        # orig_sizes = pd.DataFrame({'Original_Size' : self.term_sizes}, index=self.terms)
        # ont.update_node_attr(orig_sizes)        
        # if len(new_connections)>0:
        #     summary_sizes = pd.DataFrame({'Original_Size' : [int(x.split('_')[1]) for x in new_nodes]}, index=new_nodes)
        #     ont.update_node_attr(summary_sizes)

        if len(new_connections) > 0:
            ont.update_node_attr(pd.DataFrame({'Label':['_'.join(x.split('_')[1:]) for x in new_nodes]}, index=new_nodes))
        
        return ont
    
    def delete(self,
               to_delete=None,
               to_keep=None,
               preserve_transitivity=True,
               inplace=False):
        """Delete genes and/or terms from the ontology.

        Parameters
        ----------
        to_delete : array-like (optional)
           
            Names of genes and/or terms to delete. Either to_delete or
            to_keep must be specified.

        to_keep : array-like (optional)
           
            Names of genes and/or terms to keep; all other genes/terms
            are delete. Only used if to_delete is not specified.

        preserve_transitivity : bool

            If True, then maintain transitive relations when deleting
            terms. For example, if the hierarchical structure consists
            of

            geneA --> term1
            term1 --> term2
            term2 --> term3
            term2 --> term4

            then deleting term2 will result in the structure:

            geneA --> term1
            term1 --> term3
            term3 --> term4

            If False, then deleting term2 will result in a
            disconnected structure:

            geneA --> term1            

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

        if to_delete is not None:
            terms = set([x for x in to_delete if x in ont.terms_index])
            genes = set([x for x in to_delete if x in ont.genes_index])
        elif to_keep is not None:
            terms = set(ont.terms) - set([x for x in to_keep if x in ont.terms_index])
            genes = set(ont.genes) - set([x for x in to_keep if x in ont.genes_index])
        else:
            raise Exception('Must specify nodes to delete or to keep')

        if len(genes) > 0:
            ont.genes = [g for g in ont.genes if g not in genes]
            ont.genes_index = make_index(ont.genes)
            ont.gene_2_term = {g : t for g, t in ont.gene_2_term.items()
                               if g not in genes}
            ont._update_fields()
            
        if len(terms) > 0:
            if preserve_transitivity:
                gene_2_term_set = {g : set([ont.terms[s] for s in t]) for g, t in ont.gene_2_term.items()}
                term_2_gene_set = {a : set(b) for a, b in ont.term_2_gene.items()}
                child_2_parent_set = {a : set(b) for a, b in ont.child_2_parent.items()}
                parent_2_child_set = {a : set(b) for a, b in ont.parent_2_child.items()}
                for t in terms:
                    t_parents = child_2_parent_set[t]
                    t_genes = term_2_gene_set[t]
                    t_children = parent_2_child_set[t]
                    for g_i in t_genes:
                        g = ont.genes[g_i]
                        gene_2_term_set[g].update(t_parents)
                        gene_2_term_set[g].remove(t)
                    for p in t_parents:
                        term_2_gene_set[p].update(t_genes)                        
                        parent_2_child_set[p].update(t_children)
                        parent_2_child_set[p].remove(t)
                    for c in t_children:
                        child_2_parent_set[c].update(t_parents)
                        child_2_parent_set[c].remove(t)
                    del child_2_parent_set[t]
                    del parent_2_child_set[t]
                    del term_2_gene_set[t]

                ont.terms = [t for t in ont.terms if t not in terms]
                terms_index = make_index(ont.terms)
                ont.terms_index = terms_index
                ont.gene_2_term = {g : sorted([terms_index[s] for s in t]) for g, t in gene_2_term_set.items()}
                ont.child_2_parent = {c : sorted(p) for c, p in child_2_parent_set.items()}
                ont.parent_2_child = invert_dict(ont.child_2_parent)
                ont._update_fields()                
            else:
                tmp_gene_2_term = {g : [ont.terms[t] for t in t_list]
                                   for g, t_list in ont.gene_2_term.items()}
                ont.terms = [t for t in ont.terms if t not in terms]
                ont.terms_index = make_index(ont.terms)
                ont.gene_2_term = {g : [ont.terms_index[t] for t in t_list if t not in terms]
                                   for g, t_list in tmp_gene_2_term.items()}
                ont.parent_2_child = {p : [c for c in c_list if c not in terms]
                                      for p, c_list in ont.parent_2_child.items()
                                      if p not in terms}
                ont._update_fields()

        # Update node/edge attributes
        to_keep = (set(ont.terms) | set(ont.genes)) - genes - terms                    
        ont.edge_attr = ont.edge_attr[ont.edge_attr.index.get_level_values(0).isin(to_keep) | \
                                      ont.edge_attr.index.get_level_values(1).isin(to_keep)]
        ont.node_attr = ont.node_attr[ont.node_attr.index.isin(to_keep)]
                
        return ont
        
    def rename(self,
               genes=lambda x: x,
               terms=lambda x: x,
               inplace=False):
        """Rename gene and/or term names. 

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
            new_genes = set()
            new_gene_2_term = {}
            for g in ont.genes:
                new_g = genes.get(g, g)
                if hasattr(new_g, '__iter__') and not isinstance(new_g, str):
                    for new_gg in new_g:
                        new_genes.add(new_gg)
                        new_gene_2_term[new_gg] = ont.gene_2_term[g]
                else:
                    new_genes.add(new_g)
                    new_gene_2_term[new_g] = ont.gene_2_term[g]

            ont.genes = sorted(new_genes)
            ont.gene_2_term = new_gene_2_term
            ont.genes_index = make_index(ont.genes)
            ont._update_fields()
        if terms:
            ont.parent_2_child = {terms.get(p, p) : [terms.get(c, c) for c in c_list]
                                  for p, c_list in ont.parent_2_child.items()}
            old_term_names = ont.terms
            ont.terms = [terms.get(t,t) for t in ont.terms]

            # Retain a unique set of term names
            ont.terms = sorted(set(ont.terms))
            ont.terms_index = make_index(ont.terms)

            ont.gene_2_term = {g : [ont.terms_index[terms.get(t,t)] for t in [old_term_names[t] for t in t_list]] for g, t_list in ont.gene_2_term.items()}

            ont._update_fields()

        conversions = genes.copy()
        conversions.update(terms)
        # Remove identities
        conversions = {k : v for k, v in conversions.items() if k!=v}
        
        f = lambda x: conversions.get(x,x)
        
        # Update node attributes
        index = ont.node_attr.index
        ont.node_attr.index = pd.Series(index).map(f)

        # Update edge attributes
        idx = ont.edge_attr.index
        idx.set_levels([idx.levels[0].map(f), idx.levels[1].map(f)], inplace=True)

        ont._check_valid()
        
        return ont

    def _check_valid(self):
        if not self.is_dag():
            print('Found cycle:', nx.find_cycle(self._to_networkx_no_layout()))
            raise Exception('Not a directed acyclic graph')

        assert len(self.genes) == len(set(self.genes))
        assert len(self.terms) == len(set(self.terms))
        assert set(self.genes) == set(self.gene_2_term.keys())
        assert set(self.terms) == set(self.child_2_parent.keys())
        assert set(self.terms) == set(self.parent_2_child.keys())
        assert set(self.terms) == set(self.term_2_gene.keys())
        assert self.edge_attr.index.duplicated().sum()==0
        assert self.node_attr.index.duplicated().sum()==0
        
    def to_table(self,
                 output=None,
                 term_2_term=True,
                 gene_2_term=True,
                 edge_attr=False,
                 header=True,
                 parent_child=True,
                 clixo_format=False):
        """Convert Ontology to a table representation. Return a
        pandas.DataFrame and, optionally, write it to a file as a
        tab-delimited file.

        Parameters
        ----------

        output : filepath or file-like

            File to write table. If None, then only return a
            pandas.DataFrame

        term_2_term : bool

            Include (child term, parent term) pairs

        gene_2_term : bool

            Include (gene, term) pairs

        edge_attr : array-like or bool

            List of extra edge attributes to include. If True, then
            include all attributes. If False, then don't include any
            attribute.

        header : bool

            If True (default), then write the column names as the
            first row of the table.

        parent_child : bool

            If True, then the first column is the parent term and the
            second column is the child term or gene. If False, then
            the columns are reversed.

        clixo_format : bool

            If True, the table is the same format used the CLIXO C++
            implementation. In particular, the table has three columns:

            Column 1) Parent Term
            Column 2) Child Term or Gene
            Column 3) The string "gene" if the row is a
                      gene-term mapping, otherwise the string "default".

        Returns
        -------

        : pandas.DataFrame

            Contains at least three columns: (1) "Parent", (2)
            "Child", and (3) "EdgeType".

        """

        if clixo_format:
            df = self.to_table(output=None,
                               term_2_term=True,
                               gene_2_term=True,
                               edge_attr=False,
                               header=False,
                               parent_child=True,
                               clixo_format=False)

            df.replace({self.EDGETYPE_ATTR : {self.GENE_TERM_EDGETYPE : 'gene', self.CHILD_PARENT_EDGETYPE : 'default'}}, inplace=True)
            
            if output is not None:
                df.to_csv(output, header=False, index=False, sep='\t')
            return df

        df = pd.DataFrame(columns=['Parent','Child',self.EDGETYPE_ATTR])
        if term_2_term:
            df = df.append(self._hierarchy_to_pandas(), ignore_index=True)
        if gene_2_term:
            df = df.append(self._mapping_to_pandas(), ignore_index=True)

        if edge_attr and self.edge_attr.shape[1] > 0:
            if edge_attr==True:
                edge_attr = df.columns
            df = df.merge(self.edge_attr,
                          how='left',
                          left_on=['Child', 'Parent'],
                          right_index=True)

        first_two = ['Parent', 'Child'] if parent_child else ['Child', 'Parent']
        df = df[first_two + [x for x in df.columns if x not in first_two]]

        if output is not None:
            df.to_csv(output, header=header, index=False, sep='\t')        
        
        return df

    def _hierarchy_to_pandas(self):
        triples = [(p,c) for p, c_list in self.parent_2_child.items() for c in c_list]
        df = pd.DataFrame(triples, columns=['Parent', 'Child'])
        df[self.EDGETYPE_ATTR] = self.CHILD_PARENT_EDGETYPE                  
        return df

    def _mapping_to_pandas(self):
        pairs = [(self.terms[t], g) for g, t_list in self.gene_2_term.items() for t in t_list]
        df = pd.DataFrame(pairs, columns=['Parent', 'Child'])
        df[self.EDGETYPE_ATTR] = self.GENE_TERM_EDGETYPE
        return df

    def copy(self):
        """Create a deep copy of the Ontology object"""

        ont = Ontology(None, None, **{'empty' : True})
        
        for x in ['node_attr', 'edge_attr']:
            setattr(ont, x, getattr(self, x).copy())
        for x in ['genes', 'terms']:
            setattr(ont, x, getattr(self, x)[:])
        if self._term_sizes is None:
            ont._term_sizes = None
        else:
            ont._term_sizes = self._term_sizes[:]
        for x in ['genes_index', 'terms_index']:
            setattr(ont, x, getattr(self, x).copy())            
        for x in ['gene_2_term', 'term_2_gene', 'child_2_parent', 'parent_2_child']:
            copy_val = {k : v[:] for k, v in getattr(self, x).items()}
            setattr(ont, x, copy_val)
    
        return ont

    def flatten(self,
                include_genes=True,
                include_terms=False,
                similarity='Resnik'):
        """Flatten the hierarchy into a node-node similarity matrix by
        calculating a similarity between pair of genes in
        `genes_subset`. Currently, only the Resnik semantic similarity
        measure is implemented.
        
        Parameters
        -----------        

        include_genes : bool

            If True, then calculate pairwise similarities between
            genes. If `include_terms` is also True, then also
            calculate similarities between genes and terms.

        include_terms : bool

            If True, then calculate pairwise similarities between
            terms.  If `include_genes` is also True, then also
            calculate similarities between genes and terms.                        
        
        similarity : str

            Type of semantic similarity. (default: "Resnik")

            The Resnik similarity s(g1,g2) is defined as
            :math:`-log_2(|T_{sca}| / |T_{root}|)` where :math:`|T|` is
            the number of genes in `genes_subset` that are under term
            T. :math:`T_{sca}` is the "smallest common ancestor", the
            common ancestral term with the smallest term
            size. :math:`T_{root}` is the root term of the ontology.

            Resnik, P. (1999). Semantic similarity in a taxonomy: An
            information-based measured and its application to problems
            of ambiguity in natural
            language. J. Artif. Intell. Res. 11,95-130.

        Returns
        -------

        : (sim, nodes)

            A 2-tuple consisting of `sim`, a node-by-node NumPy array,
            and `nodes`, a NumPy array of the node names in `sim`.

        """

        assert include_genes

        assert not include_terms, 'include_terms is not yet implemented'
            
        if similarity=='Resnik':
            sca, nodes = self.get_best_ancestors(include_genes=include_genes)    
            nodes_subset = self.genes if include_genes else []
            nodes_subset += self.terms if include_terms else []
            nodes_idx = ddot.utils.make_index(nodes)
            idx = [nodes_idx[v] for v in nodes_subset]
            sca = sca[idx, :][:, idx]

            ss = -1 * np.log2(np.array(self.term_sizes)[sca] / float(len(self.genes)))
            ss = ss.astype(np.float32)
            return ss, np.array(nodes_subset)
        else:
            raise Exception('Unsupported similarity type')
    
    def common_ancestors(self, nodes, min_nodes='all', minimal=True):
        """Return the common ancestors of a set of genes

        Parameters
        ----------
        
        nodes : list
        
            List of nodes (genes and/or terms) to find the common ancestors

        min_nodes : str or int

            If 'all', then return only terms that contain all of the
            input genes. If an integer, then return only terms that
            contain at least <nodes> of the input genes

        minimal : bool

            If True, then do NOT return the terms that are themselves
            ancestors of the other common ancestors. This filter
            leaves only the 'minimal' set of common ancestors.

        Returns
        -------

        : list
        
            List of common ancestors

        """
        
        if min_nodes=='all':
            min_nodes = len(nodes)

        conn = self.connected(nodes, self.terms)
        anc_bool = conn.sum(0) >= min_nodes
        anc = np.array(self.terms)[anc_bool]

        if minimal:
            anc_conn = self.connected(anc, anc, sparse=False)
            np.fill_diagonal(anc_conn, 0)
            anc = anc[anc_conn.sum(0) == 0]

        return anc

    def _get_term_2_gene(self, verbose=False): 
        if verbose: print('Calculating term_2_gene')

        term_2_gene = invert_dict(
            self.gene_2_term,
            keymap=make_index(self.genes),
            valmap=dict(enumerate(self.terms)))

        for t in self.terms:
            if not t in term_2_gene:
                term_2_gene[t] = []
        return term_2_gene

    @property
    def term_sizes(self):
        if self._term_sizes is None:
            self._term_sizes = self._get_term_sizes(propagate=True)
        return self._term_sizes
    
    def _get_term_sizes(self, propagate=True):
        """Returns an array of term sizes in the same order as self.terms"""

        if propagate:
            ont = self.propagate(inplace=False)
            gene_2_term = ont.gene_2_term
#            gene_2_term = self._propagate_forward()
        else:
            gene_2_term = self.gene_2_term

        tmp = Counter([x for y in gene_2_term.values() for x in y])
        term_sizes = [tmp[x] for x in range(len(self.terms))]
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

    def get_tree(self, ret='edges', verbose=False):
        """Identify a spanning tree of the DAG (including genes as part of the
        DAG).

        Parameters
        ------------
        ret : str

            If 'edges', then return a list of (u, v) edges in the
            tree. If 'ontology', return an Ontology object consisting
            of only the tree edges.

        Returns
        -------
        : array-like or Ontology

        """

        tree = self.to_igraph(include_genes=True, spanning_tree=True)

        if ret=='edges':
            tree_edges = set([(tree.vs[e.source]['name'],
                               tree.vs[e.target]['name']) 
                              for e in tree.es if e['Is_Tree_Edge']=='Tree'])
            return tree_edges
        elif ret=='ontology':
            tree.delete_edges([e.index for e in tree.es if e['Is_Tree_Edge']=='Not_Tree'])
            return Ontology.from_igraph(tree, verbose=verbose)

    def is_dag(self):
        """Return True if the Ontology is a valid directed acyclic graph,
        False otherwise.

        """

        return self.to_igraph(include_genes=True, spanning_tree=False).is_dag()

    def topological_sorting(self, top_down=True, include_genes=False):
        """Perform a topological sorting.

        top_down : 
        
            If True, then ancestral nodes (e.g. the root nodes) come
            before descendants in the sorting. If False, then reverse the sorting
        """
        
        graph = self.to_igraph(include_genes=include_genes, spanning_tree=False)

        topo = list(graph.vs[graph.topological_sorting(mode='out')]['name'])

        if not top_down:
            topo = topo[::-1]

        return topo
                
    def to_igraph(self, include_genes=True, spanning_tree=False):

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
            gene_term_edges = [(self.genes_index[g], terms_index_offset[self.terms[t]])
                               for g in self.genes
                               for t in self.gene_2_term[g]]
            child_parent_edges = [(terms_index_offset[c], terms_index_offset[p]) 
                                 for p, children in self.parent_2_child.items()
                                 for c in children]

            vertex_attrs = self.node_attr.reindex(index=self.genes + self.terms).loc[self.genes + self.terms].to_dict(orient='list')
            vertex_attrs.update({
                'name':self.genes + self.terms,
                self.NODETYPE_ATTR:[self.GENE_NODETYPE for x in self.genes] + [self.TERM_NODETYPE for x in self.terms]
            })
            
            graph = igraph.Graph(n=len(self.genes) + len(self.terms),
                                 edges=gene_term_edges + child_parent_edges,
                                 directed=True,
                                 vertex_attrs=vertex_attrs,
                                 edge_attrs={self.EDGETYPE_ATTR : [self.GENE_TERM_EDGETYPE for x in gene_term_edges] + \
                                                                   [self.CHILD_PARENT_EDGETYPE for x in child_parent_edges]})                                             
        else:
            edges = [(self.terms_index[c], self.terms_index[p]) for p, children in self.parent_2_child.items() for c in children]
            graph = igraph.Graph(n=len(self.terms),
                                 edges=edges,
                                 directed=True,
                                 vertex_attrs={'name':self.terms},
                                 edge_attrs={self.EDGETYPE_ATTR : [self.CHILD_PARENT_EDGETYPE for x in edges]})
        if spanning_tree:
            parent_priority = [self.term_sizes[self.terms_index[v['name']]] if (v['name'] in self.terms_index) else 1 for v in graph.vs]

            # Identify spanning tree
            graph = self._make_tree_igraph(
                graph,
                parent_priority=parent_priority,
                optim=min,
                edge_name='Is_Tree_Edge')
            graph.es['Is_Tree_Edge'] = ['Tree' if x else 'Not_Tree' for x in graph.es['Is_Tree_Edge']]

        return graph

    def shortest_paths(self,
                       descendants=None,
                       ancestors=None,
                       sparse=False,
                       weights=None,
                       chunk_size=500):
        """Calculate the length of the shortest paths from descendant nodes to
        ancestor nodes.

        Parameters
        ----------
        sparse : bool

            If True, return a scipy.sparse matrix. If False, return a
            NumPy array

        weights : dict

            Dictionary mapping (child term, parent term) or (gene,
            term) edges to weights. Any edge with no given weight is
            assigned a weight of 0 by default.

            (default) If weights is None, then a uniform weight is
            assumed.

        chunk_size : int (optional)

            Computational optimization: shortest paths are calculated in batches.

        Returns
        -------

        d : np.ndarray or scipy.sparse.spmatrix

            d[x,y] is the length of the shortest directed path from a
            descendant node x to ancestor node y. d[x,y]==numpy.inf if
            no directed path exists. The rows are in the same order as
            <descendants>, and the columns are in the same order as
            <ancestors>.

        """
        
        graph = self.to_igraph(include_genes=True, spanning_tree=False)

        import numbers
        if weights is None:
            weights = 1
            
        if weights is not None and not isinstance(weights, numbers.Number):
            # Assume dictionary
            weights = [weights.get((graph.vs[e.source]['name'],
                                    graph.vs[e.target]['name']), 0) for e in graph.es]
        graph.es['weight'] = weights
        
        if descendants is None:
            descendants = graph.vs
        if ancestors is None:
            ancestors = descendants

        tmp = [graph.shortest_paths(
            descendants[x[0]:x[1]],
            ancestors,
            weights='weight',
            mode='out')
               for x in split_indices_chunk(len(descendants), chunk_size)]

        if sparse:
            return scipy.sparse.vstack([scipy.sparse.csr_matrix(x) for x in tmp])
        else:
            return np.vstack(tmp)
        
    def longest_paths(self,
                      descendants=None,
                      ancestors=None,
                      sparse=False,
                      weights=None,
                      chunk_size=500):
        """Computes the lengths of the longest directed paths between all pairs
        of terms.

        Returns
        -------
        d : np.ndarray or scipy.sparse.spmatrix

            d[x,y] is the length of the longest directed path from a
            descendant term with index x to an ancestral term with
            index y, where indices are defined by
            self.terms. d[x,y]==numpy.inf if no directed path exists.

        """
        
        d = self.shortest_paths(descendants=descendants,
                                ancestors=ancestors,
                                sparse=sparse,
                                weights=-1,
                                chunk_size=chunk_size)
        
        if sparse:
            d.data = -1 * d.data
        else:
            d = -1 * d

        return d

    def connected(self,
                  descendants=None,
                  ancestors=None,
                  sparse=False):
        """Calculate which genes or terms are descendants of other genes or
        terms.

        Parameters
        -----------

        descendants: list

            A list of genes and/or terms. Default: A list of all genes
            followed by a list of all terms, in the same order as
            `self.genes` and `self.terms`.

        ancestors: list

            A list of genes and/or terms. Default: Same as the
            ``descendants`` parameter.

        sparse : bool

            If True, return a scipy.sparse matrix. If False (default),
            return a NumPy array.

        Returns
        -------
        d : np.ndarray or scipy.sparse.matrix

            A descendants-by-ancestors matrix. ``d[i,j]`` is 1 if term
            i is a descendant of term j, and 0 otherwise. Note that
            ``d[i,i]==1`` and ``d[root,i]==0``, for every i.

        """
            
        d = self.shortest_paths(descendants=descendants,
                                ancestors=ancestors,
                                sparse=sparse)
        if sparse:            
            d.data = np.isfinite(d.data)
        else:
            d = np.isfinite(d)
        return d

    # def get_leaves(self, terms_list, children_list=None):
    #     """Returns terms in ``terms_list`` that are not ancestors of any term in
    #     ``children_list``.
        
    #     Parameters
    #     ----------
    #     terms_list : list

    #     children_list : list

    #         If ``children_list`` is None, then select the terms in
    #         <terms_list> that are not ancestors of any of the other
    #         terms in <terms_list>.

    #     """
        
    #     connectivity_matrix_nodiag = self.get_connectivity_matrix_nodiag()
        
    #     terms_list = np.array(terms_list)
    #     if children_list is None:
    #         children_list = terms_list
    #     else:
    #         children_list = np.array(children_list)

    #     return terms_list[~ np.any(connectivity_matrix_nodiag[children_list, :][:, terms_list], axis=0)]

    def propagate(self,
                  direction='forward',
                  gene_term=True,
                  term_term=False,
                  verbose=False,
                  inplace=False):
        """Propagate gene-term annotations through the ontology.

        As an example, consider an ontology with one gene ``g``, three terms
        ``t1, t2, t3`` and the following connections:

        ::

            t1-->t2
            t2-->t3
            g-->t1
            g-->t2

        In "forward" propagation, a new relation ``g-->t3`` is added. In
        "reverse" propagation, the relation "g-->t2" is deleted
        because it is an indirect relation inferred from "g-->t1" and
        "t1-->t2".
        
        Parameters
        ----------
        direction : str

            The direction of propagation. Either 'forward' or 'reverse'

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

        assert direction in ['forward', 'reverse'], "Propagation direction must be forward or backward"

        forward = direction=='forward'
        
        if not forward:
            # This is needed to ensure that the pruning to a parent's
            # gene set can be based on the gene sets of its direct
            # children
            ont = ont.propagate(gene_term=gene_term, term_term=term_term, direction='forward', inplace=True)

        if gene_term:
            term_2_gene_set = {t : set(g) for t, g in ont.term_2_gene.items()}
        if term_term:
            parent_2_child_set = {p : set(c) for p, c in ont.parent_2_child.items()}
            
        # # TODO: have this topological sorting be a part of the code below
        # graph = ont.to_igraph(include_genes=False, spanning_tree=False)
        # for c_idx in graph.topological_sorting(mode='in'):
        #     child = graph.vs[c_idx]['name']

        for child in ont.topological_sorting(top_down=forward, include_genes=False):
            for parent in ont.child_2_parent[child]:
                if gene_term:
                    if forward:
                        term_2_gene_set[parent] |= term_2_gene_set[child]
                    else:
                        term_2_gene_set[parent] -= term_2_gene_set[child]
                    
                if term_term:
                    if forward:
                        parent_2_child_set[parent] |= parent_2_child_set[child]
                    else:
                        parent_2_child_set[parent] -= parent_2_child_set[child]

        if gene_term:
            ont.gene_2_term = invert_dict(term_2_gene_set,
                                          keymap=make_index(ont.terms),
                                          valmap=dict(enumerate(ont.genes)))
            ont.term_2_gene = {a : list(b) for a, b in term_2_gene_set.items()}

        if term_term:
            ont.parent_2_child = {a : list(b) for a, b in parent_2_child_set.items()}
            ont.child_2_parent = ont._get_child_2_parent()

        ont._check_valid()
        
        return ont

    def get_ontotype(self,
                     genotypes,
                     input_format='gene_list',
                     output_format='dataframe',
                     matrix_columns=None):
        """Transform genotypes to ontotypes.

        .. [1] Yu, M.K., Kramer, M., Dutkowski, J., Srivas, R., Licon,
               K., Kreisberg, J.F., Ng, C.T., Krogan, N., Sharan,
               R. and Ideker, T., 2016. "Translation of genotype to
               phenotype by a hierarchy of cell subsystems". *Cell
               Systems*, 2(2), pp.77-88.

        Parameters
        ----------
        genotypes : list, np.ndarray, scipy.sparse.spmatrix, pd.DataFrame

        input_format : str

            If "gene_list", then ``genotypes`` is a list of genotypes,
            where genotype is itself a list of genes mutated. Each
            gene is assumed to have a mutation value of 1.

            If 'matrix', then ``genotypes`` is a genotype-by-gene
            matrix, where the value at position (i,j) represents the
            mutation value of gene j in genotype i. ``genotypes`` can
            be a NumPy array, SciPy sparse matrix, or Pandas
            dataframe.

        output_format : str

            If 'sparse', then return a sparse matrix as a
            scipy.sparse.csr_matrix object. (default)

            If 'dataframe', then return a pandas.DataFrame object.

            If 'array', then return a numpy.ndarray object.

        matrix_columns : list
            
            represents a list of the genes that are represented by the
            columns of ``genotypes``. Only used when input_format is
            "matrix" and ``genotypes`` is a NumPy array or SciPy sparse
            matrix.

        Returns
        -------
        : scipy.sparse.csr_matrix, pandas.DataFrame, numpy.ndarray

            genotype-by-term matrix, where the ordering of rows and
            terms is the same as ``genotypes`` and ``self.terms``

        """

        genotypes_names = None
        
        if input_format=='gene_list':
            gene_2_term = {k: np.array(v) for k, v in self.gene_2_term.items()}
            genotypes_x = [np.concatenate([gene_2_term[g] for g in gset]) if len(gset)>0 else np.array([]) for gset in genotypes]
            indices = np.concatenate(genotypes_x)
            indptr = np.append(0, np.cumsum([gset.size for gset in genotypes_x]))
            data = np.ones((indices.size, ), dtype=np.int64)
            ontotypes = scipy.sparse.csr_matrix(
                (data, indices, indptr),
                (len(genotypes), len(self.terms)))
            ontotypes.sum_duplicates()        
        elif input_format=='matrix':            
            if isinstance(genotypes, pd.DataFrame):
                matrix_columns = genotypes.columns
                genotypes_names = genotypes.index
                genotypes = genotypes.values
            elif isinstance(genotypes, np.ndarray) or scipy.sparse.issparse(genotypes):
                assert matrix_columns is not None
            else:
                raise Exception("Parameter <genotypes> must be a genotype-by-gene matrix "
                                "represented as a Pandas dataframe, NumPy array, or SciPy sparse matrix. "
                                "Consider changing the <input_format> parameter")
            contained = np.array([g in self.genes_index for g in matrix_columns])
            genotypes = scipy.sparse.csc_matrix(genotypes)[:,contained]            
            gene_2_term_matrix = scipy.sparse.csr_matrix(self.get_gene_2_term_matrix())
            gene_2_term_matrix = scipy.sparse.csr_matrix(gene_2_term_matrix)[contained,:]
            ontotypes = genotypes.dot(gene_2_term_matrix)
        else:
            raise Exception('Invalid input format')

        if output_format=='dataframe':
            ontotypes = pd.DataFrame(ontotypes.toarray(), columns=self.terms)
            if genotypes_names is not None:
                ontotypes.index = genotypes_names
        elif output_format=='sparse':
            pass
        elif output_format=='array':
            ontotypes = ontotypes.toarray()
        else:
            raise Exception('Invalid output format')
            
        return ontotypes

    def get_gene_2_term_matrix(self):
        """Returns a gene-by-term matrix stored as a scipy.sparse.coo_matrix
        
        Returns
        -------
        : scipy.sparse.coo_matrix

        """

        # Convert gene names to indices
        gene_2_term = [(self.genes_index[g], t_list) 
                       for g, t_list in self.gene_2_term.items()]

        gene_2_term_matrix = scipy.sparse.coo_matrix(
            ([1 for g, t_list in gene_2_term for t in t_list],
             ([g for g, t_list in gene_2_term for t in t_list],
              [t for g, t_list in gene_2_term for t in t_list])),
            shape=(len(self.genes), len(self.terms)))
        
        return gene_2_term_matrix

    def summary(self):
        """Summarize the Ontology's contents with respect to number of genes,
        terms, and connections.

        Returns
        --------
        : str

        """

        if self.node_attr is None:
            node_attr_names = []
        else:
            node_attr_names = self.node_attr.columns.tolist()
            # node_attr_names = ', '.join(map(str, self.node_attr.columns))

        if self.edge_attr is None:
            edge_attr_names = []
        else:
            edge_attr_names = self.edge_attr.columns.tolist()
            # edge_attr_names = ', '.join(map(str, self.edge_attr.columns))
            
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
    def infer_ontology(cls,
                       graph,
                       method,
                       **kwargs):        
        if method.lower()=='clixo_0.3':
            kwargs["clixo_version"] = 0.3
            return cls.run_clixo(graph, **kwargs)
        elif method.lower()=='clixo_1.0':
            kwargs["clixo_version"] = 1.0
            return cls.run_clixo(graph, **kwargs)
        elif method.lower()=='scipy.linkage':
            return cls.run_scipy_linkage(graph, **kwargs)
        else:
            raise Exception("Unsupported method for inferring ontology")

    @classmethod
    def run_scipy_linkage(cls,
                          graph,
                          linkage_method='single',                          
                          leaves=None,
                          **kwargs):
        """Runs scipy's agglomerative clustering methods for inferring binary
        trees. Arguments are passed directly to the
        `scipy.cluster.hierarchy.linkage` function.

        Parameters
        ----------
        graph

            Passed as the `y` argument to the `linkage` function.

        linkage_method

            Passed as the `method` argument to the `linkage` function.

        leaves : list

            A list of the gene names represented in `graph` (listed in
            the same order).

        kwargs

            Passed as extra arguments with the same argument names to the `linkage` function.
        
        Returns
        --------
        : ddot.Ontology.Ontology
        
        See `https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage` for more details.

        """

        import scipy.cluster.hierarchy
        from scipy.cluster.hierarchy import linkage

        Z = linkage(graph, method=linkage_method, **kwargs)

        ont = cls.from_scipy_linkage(Z)

        if leaves is None:
            ont.rename(terms=lambda x: 'S:%s' % x,
                       inplace=True)
        else:
            leaves_dict = {str(i) : x for i, x in enumerate(leaves)}
            ont.rename(genes=leaves_dict,
                       terms=lambda x: 'S:%s' % x,
                       inplace=True)

        return ont
    
    @classmethod
    def run_clixo(cls,
                  graph,
                  alpha=0.0,
                  beta=None,
                  newman_modularity=None,
                  miyauchi_modularity=None,
                  stop_score=None,
                  min_dt=-10000000,
                  timeout=100000000,
                  square=False,
                  square_names=None,
                  output=None,
                  output_log=None,
                  clixo_cmd=None,
                  clixo_version=None,
                  verbose=False, 
                  debug=False):
        """Runs the CLIXO algorithm and returns the result as an Ontology object.

        Acts as a wrapper for the C++ packages for CLIXO v0.3 (`https://mhk7.github.io/clixo_0.3/`) and v1.0 (`https://github.com/fanzheng10/CliXO`).
        
        Parameters
        ----------

        graph : pandas.DataFrame

           Gene-gene similarities represented as a 3-column
           pandas.DataFrame (sparse) or square-shaped pandas.DataFrame
           (square format)

        alpha : float

           CLIXO alpha parameter

        beta : float

           CLIXO beta parameter

        min_dt : float

           Minimum similarity score

        timeout : int

           Maximum time (in seconds) allowed to run CLIXO

        square : bool
        
           If True, then <graph> is interpreted as a square-shaped
           DataFrame. Otherwise, it is interpreted as a 3-column
           DataFrame.

        square_names : array-like (optional)

           The names of the rows and columns in <graph>. Only used
           when square==True.

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
                print('temp output:', output)
            rerun, delete_output = True, True
        else:
            delete_output = False
                
        if not (isinstance(graph, str) and os.path.exists(graph)):

            if square:
                assert graph.shape[1] == graph.shape[0]
                graph_sq = pd.DataFrame(graph, index=square_names, columns=square_names)
                graph = melt_square(graph_sq)

            # Write graph into a temporary file.
            # Assumes that <graph> is a list of 3-tuples (parent, child, score)
            with tempfile.NamedTemporaryFile('w', delete=False) as graph_file:
                try:
                    graph.to_csv(graph_file, sep='\t', header=False, index=False)
                except:
                    graph_file.write('\n'.join(['\t'.join([str(x[0]),str(x[1]),str(x[2])]) for x in graph]) + '\n')
            graph = graph_file.name
            if verbose:
                print('temp graph:', graph)
            rerun, delete_graph = True, True
        else:
            delete_graph = False

        if debug:
            delete_graph = False
            
        if not (isinstance(output_log, str) and os.path.exists(output_log)):
            output_log_file = tempfile.NamedTemporaryFile('w', delete=False)
            output_log = output_log_file.name
            if verbose:
                print('temp output log:', output_log)
            rerun, delete_output_log = True, True
        else:
            delete_output_log = False

        if debug:
            delete_output_log = False

        if rerun:
            try:
                return cls.run_clixo(
                    graph,
                    alpha=alpha,
                    beta=beta,
                    newman_modularity=newman_modularity,
                    miyauchi_modularity=miyauchi_modularity,
                    stop_score=stop_score,
                    min_dt=min_dt,
                    timeout=timeout,
                    output=output,
                    output_log=output_log,
                    clixo_cmd=clixo_cmd,
                    clixo_version=clixo_version,
                    verbose=verbose,
                    debug=debug)
            finally:
                if delete_output:
                    os.remove(output)
                if delete_output_log:
                    os.remove(output_log)
                if delete_graph:
                    os.remove(graph)

        top_level = os.path.dirname(os.path.abspath(inspect.getfile(ddot)))

        if clixo_version is None:
            clixo_version = 0.3

        if clixo_version == 0.3:
            if clixo_cmd is None:
                clixo_cmd = os.path.join(top_level, 'clixo_0.3', 'clixo')

            if beta is None:
                beta = 1.0
                
            cmd = ("""{0} {1} {2} {3} """.format(clixo_cmd, graph, alpha, beta) +
                   """ | tee {}""".format(output_log))            
        elif clixo_version == 1.0:
            if clixo_cmd is None:
                raise Exception(
                    """Please install CliXO 1.0 from https://github.com/fanzheng10/CliXO,"""
                    """and specify the path on your machine to the CliXO executable file.""")

            if alpha is None:
                raise Exception(
                    """For CliXO v1.0, alpha must be set to a positive (non-zero) value.""")
                
            # Format optional parameters
            if beta is None:
                beta = ""
            else:
                beta = '-b %s' % beta                
            if newman_modularity is None:
                newman_modularity = ""            
            else:
                newman_modularity = '-m %s' % newman_modularity
            if miyauchi_modularity is None:
                miyauchi_modularity = ""                
            else:
                miyauchi_modularity = '-z %s' % miyauchi_modularity
            if stop_score is None:
                stop_score = ""
            else:
                stop_score = '-s %s' % stop_score

            cmd = ("""{} -i {} -a {} {} {} {} {}""".format(clixo_cmd, graph, alpha, beta, newman_modularity, miyauchi_modularity, stop_score) + 
                   """ | tee {}""".format(output_log))
            
        if verbose:
            print('CLIXO command:', cmd)

        p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, bufsize=1)
        
        curr_dt = None
        start = time.time()

        ## Asynchronous readout of the command's output
        while p.poll() is None:

            # Break if passed the maximum processing time
            if time.time() - start > timeout:
                if verbose:
                    time_print(
                        ('Killing process %s (OUT OF TIME). '
                         'Current dt: %s: Output: %s') % (p.pid, curr_dt, output_log))
                break

            line = p.stdout.readline().decode().rstrip()
            
            # Break if the min_dt has been met
            if '# dt: ' in line:
                curr_dt = float(line.split('# dt: ')[1])
                if curr_dt < min_dt:
                    if verbose:
                        time_print(
                            ('Killing process %s (BEYOND MIN THRESHOLD). '
                             'Current dt: %s, min_dt: %s') % (p.pid, curr_dt, min_dt))
                    break

            # datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            # If line was empty, then sleep a bit
            if line=='':  time.sleep(0.01)

        # if not os.path.isfile(table):
        #     raise Exception(
        #         "CliXO script did not produce a file. One of the common reasons is the"
        #         "choice of alpha and beta parmaeters. Consider setting"
        #         "beta=0.5 and alpha to be between 1% to 1% of the"
        #         "distance from the minimum or maximum gene similarity"
        #         "scores. For example, if the gene similarities range"
        #         "between 0 and 1, then consider alpha's between 0.01"
        #         "and 0.1")

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

        if verbose: time_print('Elapsed time (sec): %s' % (time.time() - start))
        
        ont = cls.from_table(output, clixo_format=True)
        ont.rename(terms=lambda x: 'S:%s' % x, inplace=True)
        ont.edge_attr.rename(columns={'3':'CLIXO_score'}, inplace=True)

        # Copy the "CliXO score" edge attribute into a "Parent weight" node attribute
        ont.node_attr = ont.edge_attr.copy()
        ont.node_attr.index = ont.node_attr.index.droplevel(0)
        ont.node_attr.rename(columns={'CLIXO_score':'Parent weight'}, inplace=True)
        ont.node_attr.index.name = 'Node'
        ont.node_attr = ont.node_attr.loc[~ ont.node_attr.index.duplicated()]

        if verbose: print('Ontology:', ont)

        return ont

    def to_ndex(self,
                ndex_user,
                ndex_pass,
                ndex_server=None,
                name=None,
                description=None,
                network=None,
                main_feature=None,
                subnet_max_term_size=None,
                visible_term_attr=None,
                layout='bubble',
                propagate='reverse',
                style=None,
                node_alias='Original_Name',
                term_2_uuid=None,
                visibility='PUBLIC',
                verbose=False):
        """Upload an Ontology object to NDEx. The Ontology can be preformatted in
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
            Ontology as a graph. Node positions are stored in the
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

        node_alias : str

        visibility : str


        Returns
        -------
        : ndex.networkn.NdexGraph

        """

        if propagate is not None:
            ont = self.propagate(direction=propagate, inplace=False)
        else:
            ont = self

        if ndex_server is None:
            ndex_server = ddot.config.ndex_server
            
        if (network is not None) and (term_2_uuid is None):
            if subnet_max_term_size is None:
                terms = ont.terms
            else:
                terms = [t for t,s in zip(ont.terms, ont.term_sizes) if s <= subnet_max_term_size]
            
            # Only upload subnets for the unique set of the original
            # terms
            if node_alias in ont.node_attr.columns:
                orig_2_new = {a : b.index.values for a, b in ont.node_attr.loc[terms, [node_alias]].groupby(node_alias)}
                terms = [b[0] for b in orig_2_new.values()]
                        
            term_2_uuid = ont.upload_subnets_ndex(
                network,
                main_feature,
                name,
                ndex_user,
                ndex_pass,
                ndex_server=ndex_server,
                terms=terms,
                visibility=visibility,
                verbose=verbose
            )

            if node_alias in ont.node_attr.columns:
                term_2_uuid = {s : term_2_uuid[orig_2_new[t][0]] for t in orig_2_new for s in orig_2_new[t] if orig_2_new[t][0] in term_2_uuid}
        elif term_2_uuid is None:
            term_2_uuid = {}

        if verbose: print('Creating NdexGraph')
        G = ont.to_NdexGraph(
                name=name,
                description=description,
                term_2_uuid=term_2_uuid,
                layout=layout,
                style=style)

        if visible_term_attr is not None:
            df = ddot.utils.nx_nodes_to_pandas(G, visible_term_attr)
            df.rename(columns=lambda x: 'Display:' + x, inplace=True)
            ddot.utils.set_node_attributes_from_pandas(G, df)
            G.set_network_attribute('Display', '|'.join(visible_term_attr))
            
        if verbose:  print('Uploading to NDEx')
        ont_url = G.upload_to(ndex_server, ndex_user, ndex_pass, visibility=visibility)

        return ont_url, G
            
    def to_NdexGraph(self,
                     name=None,
                     description=None,
                     term_2_uuid=None,
                     spanning_tree=True,
                     layout='bubble',
                     style=None,
                     verbose=False):
        """Formats an Ontology object into a NetworkX object with extra node
        attributes that are accessed by the hierarchical viewer.

        Parameters
        -----------
        name : str
        
            Name of Ontology, as would appear if uploaded to NDEx.

        description : str

            Description of Ontology, as would appear if uploaded to NDEx.

        term_2_uuid : dict

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
        G = self.to_networkx(layout=layout, spanning_tree=spanning_tree)
                
        if style is None:
            style = 'passthrough'
        
        # Set extra attributes for passthrough visual styling
        if style=='passthrough':
            for v, data in G.nodes(data=True):
                is_gene = data[self.NODETYPE_ATTR]==self.GENE_NODETYPE
                if 'Vis:Shape' not in data:
                    data['Vis:Shape'] = 'Rectangle' if is_gene else 'Circle'
                if 'Vis:Fill Color' not in data:
                    data['Vis:Fill Color'] = '#FFFFFF'
                if 'Vis:Border Paint' not in data:
                    data['Vis:Border Paint'] = '#000000'

                # # Deprecated: use of 'Vis:Size' attribute. HiView currently sets size according to 'Size' attribute
                # if 'Vis:Size' not in data:
                #     if 'Size' in data and data['Size'] > 0:
                #         data['Vis:Size'] = 20 * float(1. + np.log10(data['Size']))
                #     else:
                #         data['Vis:Size'] = 20

                # Set links to subnetworks supporting each term
                if term_2_uuid:
                    if v in term_2_uuid:
                        data['ndex:internalLink'] = '[%s](%s)' % (data['Label'], term_2_uuid[v])
                    elif ('Original_Name' in data) and (data['Original_Name'] in term_2_uuid):
                        data['ndex:internalLink'] = '[%s](%s)' % (data['Label'], term_2_uuid[data['Original_Name']])

#                 if ('Hidden' in data) and (data['Hidden']):
#                     print('hiding %s' % v)
# #                    data['Vis:Size'] = 100
#                     data['Size'] = 1

            for u, v, data in G.edges(data=True):
                if 'Vis:Visible' not in data and 'Is_Tree_Edge' in data:
                    data['Vis:Visible'] = data['Is_Tree_Edge']=='Tree'

            # # Deprecated: use of 'Vis:Size' attribute. HiView currently sets size according to 'Size' attribute
            # # Need to change this later
            # if hasattr(G, 'pos'):
            #     # Rescale node sizes so that they don't overlap
            #     base_size = 20.
            #     desired_ratio = 0.5                
            #     pos = G.pos
            #     max_ratio = 0
            #     for v, v_data in G.nodes(data=True):
            #         v_size = v_data['Vis:Size']
            #         for u in G.predecessors(v):                        
            #             u_size = G.node[u]['Vis:Size']
            #             dist = math.sqrt(sum([(a-b)**2 for a, b in zip(pos[v], pos[u])]))
            #             if dist > 0:
            #                 ratio = ((v_size + u_size) / 2.) / dist
            #                 max_ratio = max(max_ratio, ratio)
            #     if max_ratio > desired_ratio:
            #         for v, data in G.nodes(data=True):
            #             data['Vis:Size'] = (data['Vis:Size'] - base_size) * (desired_ratio / max_ratio) + base_size
                    
            style = ddot.config.get_passthrough_style()
        else:
            raise Exception('Unsupported style')
                
        # # Set links to subnetworks supporting each term
        # if term_2_uuid:
        #     for t in self.terms:        
        #         if t in term_2_uuid:
        #             G.node[t]['ndex:internalLink'] = '[%s](%s)' % (G.node[t]['Label'], term_2_uuid[t])

        # # Change Original_Name to node indices
        # name_2_idx = {data['name'] : v for v, data in G.nodes(data=True)}
        # for v, data in G.nodes(data=True):
        #     if 'Original_Name' in data and 'Hidden' in data and data['Hidden']==True:
        #         data['Original_Name'] = name_2_idx[data['Original_Name']]

        G = nx_to_NdexGraph(G)
        
        if name is not None:
            G.set_name(name)
        if description is not None:
            G.set_network_attribute('description', description)
            
        if style:
            import ndex.beta.toolbox as toolbox
            toolbox.apply_network_as_template(G, style)
            
        return G

    def to_cx(self,
              output=None,
              name=None,
              description=None,
              term_2_uuid=None,
              spanning_tree=True,              
              layout='bubble',
              style=None):
        """Formats an Ontology object into a CX file format

        Parameters
        -----------
        output : str

            Filename or file-like object to write CX file. If None,
            then CX is returned as a JSON object, but not written to a
            file.

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
        : CX representation as a JSON-like dictionary

        """

        # Convert to NdexGraph
        G = self.to_NdexGraph(name=name,
                              description=description,
                              term_2_uuid=term_2_uuid,
                              spanning_tree=spanning_tree,
                              layout=layout,
                              style=style)
        cx = G.to_cx()

        if output is not None:
            if hasattr(output, 'write'):
                json.dump(cx, output)
            else:
                with io.open(output, 'w') as f:
                    json.dump(cx, f)

        return cx        

    def to_graphml(self,
                   output,
                   layout='bubble',
                   spanning_tree=True):
        """Writes an Ontology object in graphml format.

        Parameters
        -----------
        output : str

            Filename or file-like object to write CX file. If None,
            then CX is returned as a JSON object, but not written to a
            file.

        layout : str

            Layout the genes and terms in this Ontology. Stored in the
            node attributes 'x_pos' and 'y_pos'. If None, then do not
            perform a layout.

        """

        # Convert to NetworkX
        G = self.to_NdexGraph(spanning_tree=spanning_tree,
                              layout=layout)
        
        if hasattr(output, 'write'):
            nx.write_graphml(G, output)
        else:
            with io.open(output, 'w') as f:
                nx.write_graphml(G, f)

    def _force_directed_layout(self, G):
        """Force-directed layout on only the terms"""

        sub_nx = G.copy()
        sub_nx.remove_edges_from([(u,v) for u,v,attr in sub_nx.edges(data=True) if attr['Is_Tree_Edge']=='Not_Tree'])
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
                            main_feature,
                            name,
                            ndex_user,
                            ndex_pass,
                            ndex_server=None,
                            terms=None,
                            gene_columns=['Gene1', 'Gene2'],
                            propagate='forward',
                            visibility='PUBLIC',
                            node_attr=None,
                            node_alias='Original_Name',
                            z_score=False,
                            spring_feature=None, spring_weight=1.0,
                            edge_groups=None,
                            max_num_edges = -1,                            
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
            Ontology.propagate(). If None, then don't
            propagate annotations.
        
        public : bool

            Whether to make networks public on NDEx
        
        node_attr : pandas.DataFrame

        """

        if propagate:
            ont = self.propagate(direction=propagate, inplace=False)
        else:
            ont = self

        if ndex_server is None:
            ndex_server = ddot.config.ndex_server

        ndex = nc.Ndex(ndex_server, ndex_user, ndex_pass)
        term_2_uuid = {}
       
        start = time.time()
        g1, g2 = gene_columns[0] + '_lex', gene_columns[1] + '_lex'

        features = [f for f in network.columns if f not in gene_columns]
        assert main_feature in features, 'A main feature of the network must be specified'

        network = network[features + gene_columns].copy()
        network[gene_columns[0]] = network[gene_columns[0]].astype(str)
        network[gene_columns[1]] = network[gene_columns[1]].astype(str)
                
        # Filter dataframe for gene pairs within the ontology
        genes_set = set(ont.genes)
        tmp = [x in genes_set and y in genes_set
               for x, y in zip(network[gene_columns[0]], network[gene_columns[1]])]
        network = network.loc[tmp, :]
        
        # Lexicographically sort gene1 and gene2 so that gene1 < gene2
        network[g1], network[g2] = zip(*[(x,y) if x<y else (y,x) for x, y in zip(network[gene_columns[0]], network[gene_columns[1]])])
        network_idx = {x : i for i, x in enumerate(zip(network[g1], network[g2]))}

        if z_score:
            for feat in features:
                network[feat] = network[feat].astype(np.float64)

            # Normalize features into z-scores
            tmp = network[features]
            network[features] = (tmp - tmp.mean()) / tmp.std()

#        network_sq = ddot.utils.pivot_square(network, g1, g2, main_feature)
            
        # Calculate the min/max range of features
        numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
        def f(x):
            if str(x) in numerics:
                return 'numeric'
            elif str(x) == 'bool':
                return 'boolean'
            else:
                raise Exception()
        feature_types = network[features].dtypes.map(f)
        feature_mins = network[features].min().astype(np.str)
        feature_maxs = network[features].max().astype(np.str)

        # set an upper limit to the maximum number of edges uploaded to NDEx
        # (contributed by Fan Zheng)
        if max_num_edges > 0:
            network.sort_values(by = main_feature, ascending=False, inplace=True)
            network = network.iloc[:max_num_edges, :]

            # Lexicographically sort gene1 and gene2 so that gene1 < gene2
            # actually this may be redundant
        network[g1], network[g2] = zip(
            *[(x, y) if x < y else (y, x) for x, y in zip(network[gene_columns[0]], network[gene_columns[1]])])
        network_idx = {x: i for i, x in enumerate(zip(network[g1], network[g2]))}

        if terms is None:
            terms = ont.terms

        if verbose: print('Uploading %s terms' % len(terms))

        for upload_idx, t in enumerate(terms):
            start = time.time()

            if node_alias in ont.node_attr.columns:
                genes = ont.node_attr.loc[genes, node_alias].values
            else:
                genes = [ont.genes[g] for g in ont.term_2_gene[t]]
            
            genes.sort()
            gene_pairs_idx = [network_idx[gp] for gp in itertools.combinations(genes, 2) \
                              if gp in network_idx]

            # New (Parent weight)
            # (contributed by Fan Zheng)
            children = ont.parent_2_child[t]
            min_children_term_weights = -1
            format_parent_weight = ('Parent weight' in ont.node_attr.columns.tolist()) and (len(children) >0)
            if format_parent_weight:                
                children_term_weights = []
                for c in children:
                    if ont.node_attr.loc[c, 'Parent weight'] >0:
                        children_term_weights.append(ont.node_attr.loc[c, 'Parent weight'])
                if len(children_term_weights):
                    children_term_weights = np.array(children_term_weights)
                    min_children_term_weights = np.min(children_term_weights)

            if len(gene_pairs_idx) > 0:
                network_sub = network.iloc[gene_pairs_idx, :]

                # New: apply some minimum string force so nodes will not fly away
                if spring_feature != None:
                    network_sub.loc[network_sub[spring_feature] < min_children_term_weights, spring_feature] = 0.5*min_children_term_weights
                    network_sub[spring_feature] = network_sub[spring_feature] ** spring_weight

                G_nx = nx.from_pandas_dataframe(network_sub, g1, g2,
                                                edge_attr=features)
                if node_attr is not None:
                    set_node_attributes_from_pandas(G_nx, node_attr)

                G_nx.add_nodes_from(list(set(genes) - set(G_nx.nodes())))
                
                # Annotate the membership of each gene to every child term
                children = ont.parent_2_child[t]
                df = pd.DataFrame({c : False for c in children}, index=genes, dtype=bool)                
                for c in children:
                    genes_in = [ont.genes[g] for g in ont.term_2_gene[c]]
                    df.loc[genes_in, c] = True
                    
                # # For each gene that is directly connected to this
                # # term, create a "Group" that contains only that gene.
                # direct_genes = set(ont.term_2_gene[t]) - set(g for c in children for g in ont.term_2_gene[c])
                # for g in direct_genes:
                #     g_name = ont.genes[g]
                #     df[g_name] = False
                #     df.loc[g_name, g_name] = True
                df.rename(columns=lambda x: 'Group:'+x, inplace=True)
                ddot.utils.set_node_attributes_from_pandas(G_nx, df)
                
                # # If a gene belongs to multiple children, then place it where it is most similar
                # for g_i in (df.sum(1) > 0).nonzero():
                #     g = genes[g_i]
                #     choices = df.loc[g, :].nonzero()
                #     network_sq.loc[g, :].argmax()                    

                G = nx_to_NdexGraph(G_nx)
                G.set_name('%s supporting network for %s' % (name, t))
                G.set_network_attribute('description', '%s supporting network for %s' % (name, t))
                G.set_network_attribute('Main Feature', main_feature)
                for f in features:

                    G.set_network_attribute('%s type' % f, feature_types[f])
                    if feature_types[f] == 'numeric':
                        G.set_network_attribute('%s min' % f, feature_mins[f])
                        G.set_network_attribute('%s max' % f, feature_maxs[f])
#                for c in children:
#                    G.set_network_attribute('Group:' + c, True)
                G.set_network_attribute('Group', '|'.join(children))

                if format_parent_weight:
                    # New: calculate the score threshold of this subnetwork
                    # (contributed by Fan Zheng)
                    G.set_network_attribute('Main Feature Default Cutoff', 0.4)
                    G.set_network_attribute('Parent weight', float(ont.node_attr.loc[t, 'Parent weight']))

                    if min_children_term_weights > 0:
                        G.set_network_attribute('Children weight', '|'.join(['{:.3f}'.format(w) for w in children_term_weights]))
                        G.set_network_attribute('Main Feature Default Cutoff', float(min_children_term_weights))

                # New: annotate edge groups
                # (contributed by Fan Zheng)
                if isinstance(edge_groups, dict) and (len(edge_groups.keys()) > 0):
                    edge_group_string = []
                    for k, vs in edge_groups.iteritems():
                        vs.sort()
                        edge_group_string.append(','.join([k] + vs))
                    edge_group_string = '|'.join(edge_group_string)
                    G.set_network_attribute('edge groups', edge_group_string)


                # New: compute a pre-layout to networks
                # (contributed by Fan Zheng)
                if spring_feature != None:
                    G_cx = G.to_cx()
                    G = NdexGraph(G_cx)
                    layouts.apply_directed_flow_layout(G, node_width=50, weight=spring_feature)


                start_upload = time.time()
                ndex_url = G.upload_to(ndex_server, ndex_user, ndex_pass, visibility=visibility)
                term_2_uuid[t] = parse_ndex_uuid(ndex_url)
                upload_time = time.time() - start_upload

                if verbose:
                    print(upload_idx,
                          'Term:', t,
                          'Gene pairs:', len(gene_pairs_idx),
                          'Genes:', len(genes),
                          'Time:', round(time.time() - start, 4),
                          'Upload time:', round(upload_time, 4),
                          'NDEx URL:', ndex_url)
            else:
                if verbose:
                    print(upload_idx, 'No data provided for gene pairs in Term: %s' % t)
                    
        return term_2_uuid

    def get_best_ancestors(self, node_order=None, verbose=False, include_genes=True):
        """Compute the 'best' ancestor for every pair of terms. 'Best' is
        specified by a ranking of terms. For example, if terms are
        ranked by size, from smallest to largest, then the smallest
        common ancestor is calculated.

        Parameters
        ----------
        node_order : list

           A list of terms, ordered by their rank with the 'best' term at the beginning.

        include_genes : bool

        Returns
        --------
        ancestors : np.ndarray

            ancestors[a,b] = the best common ancestor of terms a and
            b, represented as a 0-based index of self.terms

        nodes : list

            List of the row and column names. Rows and columns are the
            same.

        """

        ont = self.propagate(direction='reverse', inplace=False)
        graph = ont.to_igraph(include_genes=include_genes, spanning_tree=False)
        
        if node_order is None:
            # By default, sort from smallest to largest terms
            node_order = [self.terms[t] for t in np.argsort(ont.term_sizes)]
             
        d = np.int8(np.isfinite(np.array(graph.shortest_paths(graph.vs, graph.vs, mode='out'), order='C')))

        ancestor_matrix = np.zeros(d.shape, dtype=np.int32)
        ancestor_matrix.fill(-1)

        if verbose: time_print('Iterating:')
        for t in node_order:
            i = graph.vs.find(t).index
            t_i = self.terms_index[t]
            
            # Note: includes self as a child
            children = np.where(d[:,i] == 1)[0]

            # For those descendants without a computed LCA yet, set their LCA to this term
            lca_sub = ancestor_matrix[children.reshape(-1,1), children]
            lca_sub[lca_sub == -1] = t_i
            ancestor_matrix[children.reshape(-1,1), children] = lca_sub

        # Check symmetry
        assert (ancestor_matrix.T == ancestor_matrix).all()
        assert (-1 == ancestor_matrix).sum() == 0, 'The ontology may have more than one root'

        return ancestor_matrix, graph.vs['name']

    @classmethod
    def _make_tree_igraph(self,
                          graph=None,
                          method='priority',
                          edge_name='smallest_parent',
                          parent_priority=None, edge_priority=None, default_priority=None, optim='max'):
        """Returns copy of graph with new edge attribute marking spanning
        tree

        """

        if graph is None:
            graph = self.to_igraph(include_genes=False, spanning_tree=True)

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
        """Saves Ontology object with the Python pickle protocol.""" 
        pandas.io.pickle.to_pickle(self, file, compression=compression)

    @classmethod
    def read_pickle(cls, file, compression='infer'):
        """Loads an Ontology object from a pickled state.""" 
        return pandas.io.pickle.read_pickle(file, compression=compression)

    # def to_adjacency(self, include_genes=False):
    #     """Returns the adjacency matrix between genes and terms.

    #     Parameters
    #     ----------
    #     include_genes : bool
                        
    #         Include genes in the adjacency matrix such that the rows
    #         (and columns) represent both genes and terms in the order
    #         (self.genes + self.terms).

    #     Returns
    #     -------
    #     : scipy.sparse.csr_matrix

    #     """

    #     edges = [(self.terms_index[c], self.terms_index[p]) for c in self.terms for p in self.child_2_parent.get(c, [])]
    #     i, j = zip(*edges)
    #     child_2_parent_adj = coo_matrix((np.ones(len(edges), np.bool), (i,j)), shape=(len(self.terms), len(self.terms)))
    #     child_2_parent_adj = child_2_parent_adj.tocsr()

    #     adj = child_2_parent_adj

    #     if include_genes:
    #         tmp = [self.gene_2_term[g] for g in self.genes]
    #         indices = np.concatenate(tmp)
    #         indptr = np.append(0, np.cumsum([len(x) for x in tmp]))
    #         data = np.ones(indices.size, np.bool)
    #         gene_2_term_adj = csr_matrix((data,indices,indptr), shape=(len(self.genes), len(self.terms)))

    #         empty1 = csr_matrix((len(self.genes),len(self.genes))).astype(np.bool)
    #         empty2 = csr_matrix((len(self.terms),len(self.genes))).astype(np.bool)
                    
    #         hstack, vstack = scipy.sparse.hstack, scipy.sparse.vstack
            
    #         adj = vstack([hstack([empty1, gene_2_term_adj]).tocsr(),
    #                       hstack([empty2, adj]).tocsr()])
    #         adj = adj.tocsr()

    #         nodes = self.genes + self.terms
    #     else:
    #         nodes = self.terms
            
    #     return adj, nodes

    def __repr__(self):
        return self.summary()
    
    def __str__(self):
        return self.summary()
