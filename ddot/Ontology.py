from __future__ import absolute_import

import itertools, multiprocessing, logging, os, collections, random, math, sys, time, scipy, scipy.sparse, igraph
import numpy as np, pandas as pd
from itertools import groupby, combinations
from operator import *
from scipy.stats import hypergeom
from collections import Counter
import networkx as nx
import tempfile
from subprocess import Popen, PIPE, STDOUT
import inspect
import shlex
import shutil
from scipy.sparse import csr_matrix

import ddot
import ddot.config
from ddot.utils import time_print, set_node_attributes_from_pandas, set_edge_attributes_from_pandas, nx_to_NdexGraph, parse_ndex_uuid, make_index, update_nx_with_alignment, bubble_layout_nx, split_indices_chunk, invert_dict

GENE_TERM_ATTR = 'Gene_or_Term'

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
    """Align two hierarchies using alignOntology's calculateFDRs script

    Parameters
    ----------

    hier1 : ddot.Ontology.Ontology
       First ontology

    hier2 : ddot.Ontology.Ontology
       Second ontology

    iterations : int
       Number of randomized iterations

    threads : int
       Number of CPU threads

    calculateFDRs : str

       Filename of calculateFDRs script from alignOntology package. If
       None, then it is inferred based on ddot.config.
       
    output
       Filename to write alignment. If None, then don't write.
    
    Returns
    --------

    """

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
        assert os.path.isdir(ddot.config.alignOntology)
        calculateFDRs = os.path.join(ddot.config.alignOntology, 'calculateFDRs')
    assert os.path.isfile(calculateFDRs)

    output_dir = tempfile.mkdtemp(prefix='tmp')
    cmd = '{5} {0} {1} 0.05 criss_cross {2} {3} {4} gene'.format(
              hier1, hier2, output_dir, iterations, threads, calculateFDRs)
    print 'Alignment command:', cmd

    try:
        p = Popen(shlex.split(cmd), shell=False)
        p.wait()
        shutil.copy(os.path.join(output_dir, 'alignments_FDR_0.1_t_0.1'), output)
    finally:
        if os.path.isdir(output_dir):
            shutil.rmtree(output_dir)

        if p.poll() is None:
            if verbose: time_print('Killing alignment process %s. Output: %s' % (p.pid, output))
            p.kill()  # Kill the process

        alignment = read_alignment_file(output)

        print update_hier1, update_hier2

        if update_hier1:
            print 'Updating hier1'
            tmp = alignment.copy()[['Term', 'FDR', 'Similarity']]
            tmp.columns = ['Aligned_%s' % x for x in tmp.columns]            
            hier1_orig.update_node_attr(tmp)
            
        if update_hier2:
            print 'Updating hier2'
            tmp = alignment.copy()[['Term', 'FDR', 'Similarity']]
            tmp.index, tmp['Term'] = tmp['Term'], tmp.index
            tmp.columns = ['Aligned_%s' % x for x in tmp.columns]            
            hier2_orig.update_node_attr(tmp)

        return alignment

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

    return df.loc[:, ['DB Object ID', 'GO ID']].values.tolist()

class Ontology:
    """A Python representation for constructing, analyzing, and
    manipulating ontologies.

    This class's attributes and methods are focused on the ontology's
    hierarchical structure, i.e. a directed acyclic graph.

    """

    def __init__(self,
                 hierarchy,
                 mapping,
                 edge_attr=None,
                 node_attr=None,
                 parent_child=False,
                 add_root_name=None,
                 propagate=False,
                 verbose=True):
        """Construct an Ontology object.

        Parameters
        ----------
        hierarchy
            Iterable of (child term, parent term). E.g. list of 2-tuples

        mapping
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
        
        propagate : bool

            Propagate gene-term annotations up the hierarchy

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
            self.node_attr = pd.DataFrame()
        else:
            self.node_attr = node_attr

        if edge_attr is None:
            self.edge_attr = pd.DataFrame()
        else:
            self.edge_attr = edge_attr

        self._update_fields()

        if propagate:
            self.propagate_annotations(inplace=True)
            self._update_fields()

    def _update_fields(self):
        self.child_2_parent, self.child_2_parent_indices = self._get_child_2_parent()
        self.term_2_genes = self._get_term_2_genes()
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

    def update_node_attr(self, node_attr):
        ####
        # TODO : make sure that renaming/deleting/collapsing of genes and columns respect the node_attr and edge_attr

        # Filter for genes and terms in the ontology
        nodes = set(self.genes) | set(self.terms)
        node_attr = node_attr.loc[[x for x in node_attr.index if x in nodes], :]

        # Update index
        self.node_attr = self.node_attr.reindex(node_attr.index)

        # Update columns
        for col in node_attr.columns:
            self.node_attr.loc[node_attr.index, col] = node_attr[col].values

    def update_edge_attr(self, edge_attr):
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

        # Update index
        self.edge_attr = self.edge_attr.reindex(edge_attr.index)

        # Update values for overlapping columns
        for col in edge_attr.columns:
            self.edge_attr.loc[edge_attr.index, col] = edge_attr[col].values
        
    def get_roots(self):
        """Returns a list of the root term(s)"""

        return sorted(set(self.parent_2_child.keys()) - set([y for x in self.parent_2_child.values() for y in x]))

    def align(self,
              hier,
              iterations,
              threads,
              update_self=False,
              update_ref=False,
              calculateFDRs=None,
              output=None):

        return align_hierarchies(
            self,
            hier,
            iterations,
            threads,
            update_hier1=update_self,
            update_hier2=update_ref,
            calculateFDRs=calculateFDRs,
            output=output)

    def to_networkx(self,
                    node_attr=None,
                    edge_attr=None,
                    layout='bubble',
                    spanning_tree=True):
        """Converts Ontology into a NetworkX object

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

        """

        import networkx as nx
        G = nx.DiGraph()

        #################################
        ### Add nodes and node attributes

        G.add_nodes_from([g for g in self.genes])
        G.add_nodes_from([t for t in self.terms])

        term_sizes = dict(zip(self.terms, self.get_term_sizes(propagate=True)))

        for t in self.terms:
            G.node[t][GENE_TERM_ATTR] = 'Term'
            G.node[t]['Size'] = term_sizes[t]
            G.node[t]['isRoot'] = False
        for g in self.genes:
            G.node[g][GENE_TERM_ATTR] = 'Gene'
            G.node[g]['Size'] = 1
            G.node[g]['isRoot'] = False

        # Identify the root
        root = self.get_roots()[0]
        G.node[root]['isRoot'] = True

        if node_attr is not None:
            set_node_attributes_from_pandas(G, node_attr)

        if edge_attr is not None:
            set_edge_attributes_from_pandas(G, edge_attr)
            
        #################################
        ### Add edges and edge attributes

        G.add_edges_from([(g, self.terms[t],
                           dict(EdgeType='Gene-Term')) \
                          for g in self.genes for t in self.gene_2_term[g]])

        G.add_edges_from([(c, p,  
                         dict(EdgeType='Child-Parent')) \
                          for p in self.terms for c in self.parent_2_child.get(p, [])])

        if spanning_tree:
            # Identify a spanning tree
            tree_edges = self.get_tree_edges()
            nx.set_edge_attributes(G,
                                   'Is_Tree_Edge',
                                   {(s,t) : 'Tree' if ((s,t) in tree_edges) else 'Not_Tree'
                                    for s, t in G.edges_iter(data=False)})

            if layout=='bubble':
                ont_dummy = self.make_dummy()
                G_tree = ont_dummy.to_networkx(layout=None, spanning_tree=False)
                for u, v, data in G.edges(data=True):
                    if data['Is_Tree_Edge']=='Not_Tree':
                        if data['EdgeType']=='Gene-Term':
                            G_tree.remove_edge(u, 'dummy_' + v)
                        else:
                            G_tree.remove_edge(u, v)

                G.pos = bubble_layout_nx(G_tree)
                G.pos = {n : p for n, p in G.pos.items() if 'dummy_' not in n}
                nx.set_node_attributes(G, 'x_pos', {n : x for n, (x,y) in G.pos.items()})
                nx.set_node_attributes(G, 'y_pos', {n : y for n, (x,y) in G.pos.items()})
                
                # G_tree = G.copy()
                # for u, v, data in G_tree.edges(data=True):
                #     if data['Is_Tree_Edge']=='Not_Tree':
                #         G_tree.remove_edge(u, v)
                # G.pos = bubble_layout_nx(G_tree)
                # print 'nodes:', G.nodes()[:5]
                # print 'layout:', G.pos.items()[:5]
                # nx.set_node_attributes(G, 'x_pos', {n : x for n, (x,y) in G.pos.items()})
                # nx.set_node_attributes(G, 'y_pos', {n : y for n, (x,y) in G.pos.items()})

        return G

    def make_dummy(self):

        ont = self.propagate_annotations(direction='backward', inplace=False)

        new_gene_2_term = []
        new_child_2_parent = []
        for t in ont.terms:
            if len(ont.term_2_genes[t]) > 0:
                dummy_term = 'dummy_%s' % t
                for g in ont.term_2_genes[t]:
                    new_gene_2_term.append([ont.genes[g], dummy_term])
                new_child_2_parent.append([dummy_term, t])
            else:
                for g in ont.term_2_genes[t]:            
                    new_gene_2_term.append([ont.genes[g], term])
            for p in ont.child_2_parent.get(t, []):
                new_child_2_parent.append([t, p])

        ont_dummy = Ontology(new_child_2_parent, new_gene_2_term)
        return ont_dummy

    @classmethod
    def from_pandas(cls,
                    df,
                    parent_col='Parent',
                    child_col='Child',
                    is_mapping=lambda x: x['Relation']=='gene',
                    propagate=False,
                    verbose=False):

        return cls.from_table(df,
                              is_mapping=is_mapping,
                              parent_col=parent_col,
                              child_col=child_col,
                              mapping_file=None,
                              propagate=propagate,
                              verbose=verbose)

    @classmethod
    def from_table(cls,
                   table_file,
                   is_mapping=lambda x: x[2]=='gene',
                   parent_col = 0,
                   child_col = 1,
                   mapping_file=None,
                   propagate=False,
                   verbose=False):
        """Create Ontology from a tab-delimited table.

        Parameters
        ----------
        table_file : str
            Filename of table

        is_mapping : function
            Function applied on each row to determine if it is 

        parent_col : int or str
            Column for parent terms
        
        child_col : int or str
            Column for child terms and genes

        propagate : bool
            If True, then propagate gene-term mappings

        Returns
        -------
        : ddot.Ontology.Ontology

        """

        # Read table
        try:
            df = pd.read_table(table_file, comment='#', header=None)
        except:
            df = table_file

        assert isinstance(df, pd.DataFrame)
        df.loc[:, child_col] = df.loc[:, child_col].astype(str)
        df.loc[:, parent_col] = df.loc[:, parent_col].astype(str)

        # Get gene-term mappings
        if mapping_file is None:
            tmp = df.apply(is_mapping, axis=1)
            mapping = df.loc[tmp, :].loc[:,[child_col, parent_col]].values.tolist()
            tmp2 = df.loc[~ tmp, :]
        else:
            mapping = pd.read_table(mapping_file, comment='#', header=None)
            mapping.loc[:, child_col] = mapping.loc[:, child_col].astype(str)
            mapping.loc[:, parent_col] = mapping.loc[:, parent_col].astype(str)
            mapping = mapping.loc[:,[child_col,parent_col]].values.tolist()
            tmp2 = df

        # Get term-term hierarchy and attributes                                                                                                                               
        if tmp2.shape[1] > 2:
            hierarchy = tmp2.loc[:,[child_col, parent_col]].values.tolist()
            edge_attr = tmp2.loc[:,np.setdiff1d(tmp2.columns, [child_col, parent_col])]
        else:
            hierarchy = tmp2.values.tolist()
            edge_attr = None

        return cls(hierarchy,
                   mapping,
                   parent_child=False,
                   edge_attr=edge_attr,
                   propagate=propagate,
                   verbose=verbose)

    @classmethod
    def from_ndex(cls,
                  ndex_uuid,
                  ndex_server=ddot.config.ndex_server,
                  ndex_user=ddot.config.ndex_user,
                  ndex_pass=ddot.config.ndex_pass,                  
                  edgetype_attr='EdgeType',
                  edgetype_value='Gene-Term'):
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

        if '/' in ndex_uuid:
            ndex_server = parse_ndex_server(ndex_uuid)
            ndex_uuid = parse_ndex_uuid(ndex_uuid)

        from ndex.networkn import NdexGraph
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
                       edgetype_attr='EdgeType',
                       edgetype_value='Gene-Term'):
        """Converts a NdexGraph object to an Ontology object. Gene and terms
        are distinguished according by an edge attribute.

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

        hierarchy = []
        mapping = []
        for u, v, attr in G.edges_iter(data=True):
            if attr[edgetype_attr] == edgetype_value:
                mapping.append((G.node[u]['name'], G.node[v]['name']))
            else:
                hierarchy.append((G.node[u]['name'], G.node[v]['name']))

        return cls(hierarchy, mapping)

    def collapse_ontology(self,                          
                          method='mhkramer',
                          collapseRedundantNodes=None,
                          verbose=True,
                          default_relation='default',
                          min_term_size=2):
        """Remove redundant and empty terms. When a term T is removed,
        hierarchical relations are preserved by connecting every child
        of T with every parent of T. This removal operation has the
        nice property of being commutative, i.e. the order of removal
        does not matter.

        Parameters
        -----------
        method

            If "mhkramer", then use the collapseRedundantNodes script
            in the alignOntology package. If "python", then use an internal Python script

        propagate
        
            If True, propagate annotations beforehand

        min_term_size : int
        
            Remove terms that are below this size

        """

        if method=='mhkramer':

            # Propagate forward and then backwards
            ont = self.propagate_annotations(direction='backward', inplace=False)
            ont.propagate_annotations(direction='forward', inplace=True)

            # if propagate:
            #     if verbose: print 'Propagating annotations'
            #     ont = self.propagate_annotations(inplace=False)
            # else:
            #     ont = self

            # ont = self

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

            from StringIO import StringIO
            return Ontology.from_table(StringIO(collapsed))

        elif method=='python':
            g = self.to_igraph()
            if verbose: print len(g.vs), 'total nodes'

            if verbose: print 'Propagating annotations'
            self.propagate_annotations()

            parity = True
            while True:
                
                # Calculate a unique hash for every term based on its set of gense
                names_2_idx = make_index(g.vs['name'])
                term_hash = {names_2_idx[t] : (len(g_list), hash(tuple(g_list)))
                             for t, g_list in self.term_2_genes.items()
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
                        g = self.collapse_node(g,
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

    def delete_terms(self, terms_to_delete, inplace=False):
        """Delete terms from the ontology. Note that if a gene is directly
        connected to a term T but not to an ancestor of T, then
        deleting T will lose the information about the indirect
        connection to the ancestor. For this reason, it is suggested
        that propagate_annotations() is called first.

        Parameters
        ----------
        terms_to_delete : iterable of str
           
            Names of terms to delete

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
        
        terms_to_delete = set(terms_to_delete)
        tmp_gene_2_term = {g : [ont.terms[t] for t in t_list]
                           for g, t_list in ont.gene_2_term.items()}
        ont.terms = [t for t in ont.terms if t not in terms_to_delete]
        ont.terms_index = make_index(ont.terms)
        ont.gene_2_term = {g : [ont.terms_index[t] for t in t_list if t not in terms_to_delete]
                           for g, t_list in tmp_gene_2_term.items()}
        ont.parent_2_child = {p : [c for c in c_list if c not in terms_to_delete]
                              for p, c_list in ont.parent_2_child.items()
                              if p not in terms_to_delete}
        ont._update_fields()

        return ont

    def delete_genes(self, genes_to_delete, inplace=False):
        """Delete genes from ontology

        """

        if inplace:
            ont = self
        else:
            ont = self.copy()

        genes_to_delete = set(genes_to_delete)
        ont.genes = [g for g in ont.genes if g not in genes_to_delete]
        ont.genes_index = make_index(ont.genes)
        ont.gene_2_term = {g : t for g, t in ont.gene_2_term.items() 
                           if g not in genes_to_delete}
        ont._update_fields()

        return ont
        
    def rename(self,
               genes=None,
               terms=None,
               gene_prefix=None,
               term_prefix=None,
               inplace=False):
        """Rename gene and/or term names. 

        Parameters
        ----------
        genes : dict

           Dictionary mapping current gene names to new names. If
           specified, genes not in dictionary are deleted.

        terms : dict

           Dictionary mapping current term names to new names. If
           specified, terms not in dictionary are deleted.

        gene_prefix : str
        
           Prefix to append to every gene's name

        term_prefix : str
           
           Prefix to append to every term's name

        inplace : bool

           If True, then modify the ontology. If False, then create
           and modify a copy.

        Returns
        -------
        : ddot.Ontology.Ontology

        """

        if term_prefix is not None:
            assert terms is None
            terms = {t : '%s%s' % (term_prefix, t) for t in self.terms}
        if gene_prefix is not None:
            assert genes is None
            genes = {g : '%s%s' % (gene_prefix, g) for g in self.genes}

        if inplace:
            ont = self
        else:
            ont = self.copy()

        if genes is not None:
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

            # ont.genes = [genes[g] for g in ont.genes]
            # ont.genes_index = make_index(ont.genes)
            # ont.gene_2_term = {genes[g] : t for g, t in ont.gene_2_term.items()}
        if terms is not None:
            ont.delete_terms([t for t in ont.terms if not terms.has_key(t)], inplace=True)
            ont.terms = [terms[t] for t in ont.terms]
            ont.terms_index = make_index(ont.terms)
            ont.parent_2_child = {terms[p] : [terms[c] for c in c_list]
                                  for p, c_list in ont.parent_2_child.items()}
            ont._update_fields()
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
            term, (3) Relation type.  The relation type of all (gene,
            term) pairs is "gene".

        """

        df = pd.DataFrame(columns=['Parent','Child','Relation'])
        if term_2_term:
            df = df.append(self._hierarchy_to_pandas(default_relation),
                           ignore_index=True)
        if gene_2_term:
            tmp = self._mapping_to_pandas()
            tmp.rename(columns={'Gene':'Child', 'Term':'Parent'},
                       inplace=True)
            tmp['Relation'] = 'gene'
            df = df.append(tmp, ignore_index=True)
        return df

    def _hierarchy_to_pandas(self, default_relation=u'default'):
        if hasattr(self, 'relation_dict'):
            relation_dict = self.relation_dict
        else:
            relation_dict = {}

        triples = [(p,c, relation_dict.get((c, p), default_relation)) \
                   for p, c_list in self.parent_2_child.items() for c in c_list]
        df = pd.DataFrame(triples, columns=['Parent', 'Child', 'Relation'])
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
            
            Relation type for all (child term, parent term)
            pairs. Note that the relation type for (gene, term) pairs
            is "gene".

        """

        df = self.to_pandas(default_relation=default_relation)
        if parent_child:
            df = df[['Parent','Child','Relation']]
        else:
            df = df[['Child','Parent','Relation']]
        return df.to_csv(output, header=header, index=False, sep='\t')

    def copy(self):
        """Create a deep copy of the Ontology object"""

        hierarchy = [(c, p) for p, c_list in self.parent_2_child.items() for c in c_list]
        mapping = [(g, self.terms[t]) for g, t_list in self.gene_2_term.items() for t in t_list]
        return Ontology(hierarchy,
                        mapping,
                        edge_attr=None if (self.edge_attr is None) else self.edge_attr.copy(),
                        parent_child=False,
                        verbose=False)

    def transitive_closure(self):
        """Computes the transitive closure on (child term, parent term)
        relations. Transitivity rules are defined at

        http://www.geneontology.org/page/ontology-relations

        """
        
        relations = ['is_a', 'regulates', 'positively_regulates', 'negatively_regulates', 'has_part', 'part_of', 'gene']
        relations_index = make_index(relations)
        go_reasoning = dict()
        for r in relations:
            go_reasoning[('is_a', r)] = r
            go_reasoning[(r, 'is_a')] = r
            if r != 'has_part':
                go_reasoning[(r, 'part_of')] = 'part_of'
        for r in ['regulates', 'positively_regulates', 'negatively_regulates']:
            go_reasoning[(r, 'part_of')] = 'regulates'
        go_reasoning[('has_part', 'has_part')] = 'has_part'

        g = self.to_igraph()

        # Sort from leaves up
        for i in g.topological_sorting(mode='out'):
        
            # Get parents
            for j in g.neighbors(i, mode='out'):
                
                # Infer GO relations of new edges
                base_relation = g.es[g.get_eid(i, j)]['relation']

                # Iterate over grandparents
                for p in g.neighbors(j, mode='out'):
                    r = go_reasoning.get((base_relation, g.es[g.get_eid(j, p)]['relation']), None)

                    # If a relation can't be inferred, then don't add a new edge
                    if r is not None:
                        if -1 != g.get_eid(i, p, error=False):
                            # Edge already exists, so take the higher-ranked relation
                            e = g.es[g.get_eid(i, p)]
                            e['relation'] = r if relations_index[r] > relations_index[e['relation']] else e['relation']
                        else:
                            # Add new edge with relation
                            g.add_edge(i, p, relation=r)                                        

        return g
        
        ## Update parent_2_child , child_2_parent, child_2_parent_indices
        
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
            term_2_genes = self.term_2_genes
            genes_subset_set = set(genes_subset)
            term_sizes = np.array([len(set([self.genes[x] for x in term_2_genes[t]]) & genes_subset_set) for t in self.terms])
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
        
    def _get_term_2_genes(self, verbose=False): 
        if verbose: print 'Calculating term_2_genes'

        term_2_genes = invert_dict(
            self.gene_2_term,
            keymap=make_index(self.genes),
            valmap=dict(enumerate(self.terms)))

        # term_2_genes = {t : [] for t in self.terms}
        # for g, t_list in self.gene_2_term.items():
        #     for t in t_list:
        #         term_2_genes[t].append(g)
        # for t in self.terms:
        #     term_2_genes[t].sort()

        # term_2_genes = {self.terms[c]: [self.genes_index[x[0]] for x in d] \
        #                 for c, d in itertools.groupby(sorted([(a,t) for a, terms in self.gene_2_term.items() for t in terms],
        #                                                      key=lambda x:x[1]),
        #                                               key=lambda x:x[1])}

        for t in self.terms:
            if not term_2_genes.has_key(t):
                term_2_genes[t] = []
        return term_2_genes

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
                                     GENE_TERM_ATTR:['Gene' for x in self.genes] + ['Term' for x in self.terms]
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
        'backward' propagation, the relation g --> t2 is deleted
        because it is an indirect relation inferred from g --> t1 and
        t1 --> t2.

        TODO: consider 'reverse' instead of 'backward'
        
        Parameters
        ----------
        direction : str

            The direction of propgation. Either 'forward' or 'backward'

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

        elif direction=='backward':
            ont.propagate_annotations(direction='forward', inplace=True)

            term_2_genes_set = {}
            for t, g in ont.term_2_genes.items():
                term_2_genes_set[t] = set(g)

            graph = ont.to_igraph(spanning_tree=False)

            for c_idx in graph.topological_sorting(mode='in'):
                child = graph.vs[c_idx]['name']
                for parent in ont.child_2_parent.get(child, []):
                    term_2_genes_set[parent] -= term_2_genes_set[child]

            ont.gene_2_term = invert_dict(term_2_genes_set,
                                          keymap=make_index(ont.terms),
                                          valmap=dict(enumerate(ont.genes)))
            ont.term_2_genes = {a : list(b) for a, b in term_2_genes_set.items()}
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

    def get_ontotype(self, gset_sample, prop='genes', format='dict', dtype=np.int8, gene_ordering=None):
        """Calculates the ontotypes of genotypes. 

        Parameters
        ----------
        gset_sample : list

            List of gene sets (i.e. gsets) A "gset" is a set of genes,
            represented as a tuple

        prop : str

            If 'genes', then the value of a term is the number of
            genes in the term which are deleted.

        format : str

            The data format of the ontotype. Supported formats are
            'dict', 'scipy.coo'

        Returns
        -------

        If format=='dict', then return a list, where each element is a
        dictionary representation of an ontotype.

        If format=='scipy.coo', then return a genotype-by-term matrix
        as a scipy.sparse.coo_matrix

        """
        
        if prop=='genes':

            time_print('Creating features')
            print format, format=='scipy.csr'
            
            if format=='dict':
                from collections import Counter                
                onto_sample = [reduce(add, [self.gene_2_term[g] for g in gset]) for gset in gset_sample]                
                onto_sample = [dict(Counter(idxs)) for idxs in onto_sample]
                
            elif format=='scipy.coo':
                time_print('Making i, j, data')

                gene_2_term = {k: np.array(v) for k, v in self.gene_2_term.items()}
                i = np.repeat(np.arange(len(gset_sample)), [sum(gene_2_term[g].size for g in gset) for gset in gset_sample])
                j = np.concatenate([gene_2_term[g] for gset in gset_sample for g in gset])

                #ij = np.array([(a, t) for a, gset in enumerate(gset_sample) for g in gset for t in self.gene_2_term[g]])
                #data = np.ones((ij.shape[0], ))
                data = np.ones((i.size, ))
                time_print('Making sparse COO matrix')
                #onto_sample = scipy.sparse.coo_matrix((data, (ij[:,0], ij[:,1])), (len(gset_sample), len(self.terms)), dtype=dtype)
                onto_sample = scipy.sparse.coo_matrix((data, (i, j)), (len(gset_sample), len(self.terms)), dtype=dtype)

            elif format=='scipy.csr':
                time_print('Making indices, indptr, data')
                gene_2_term = {k: np.array(v) for k, v in self.gene_2_term.items()}
                gset_sample_x = [np.concatenate([gene_2_term[g] for g in gset]) if len(gset)>0 else np.array([]) for gset in gset_sample]
                indices = np.concatenate(gset_sample_x)
                indptr = np.concatenate((np.array([0]),
                                         np.cumsum([gset.size for gset in gset_sample_x])))
                data = np.ones((len(indices), ), dtype=dtype)

                time_print('Making sparse CSR matrix')
                onto_sample = scipy.sparse.csr_matrix((data, indices, indptr), (len(gset_sample), len(self.terms)), dtype=dtype)

            elif format=='dense':
                dim = len(gset_sample[0])
                assert dim <= 127, "Gset dimension must be <= 127 to create an np.int8 array"

                # Check that all genes in gset_sample are in gset_2_terms
                assert set.union(*[set(x) for x in gset_sample]).issubset(set(self.gene_2_term.keys()))

                # Convert gene_2_term to a numpy genes-by-terms boolean matrix
                time_print('Creating gene_2_term matrix')
                gene_2_term_arr = np.zeros((len(self.genes), len(self.terms)), dtype=np.int8)
                for g, terms in self.gene_2_term.items():
                    gene_2_term_arr[self.genes_index[g], terms] = 1

                time_print('Converting gset_sample into matrix')
                gset_sample_arr = np.array([[self.genes_index[g] for g in gset] for gset in gset_sample])

                time_print('Creating gset_2_terms matrix')
                onto_sample = np.zeros((len(gset_sample), len(self.terms)), dtype=np.int8)
                for d in range(dim):
                    onto_sample += np.take(gene_2_term_arr,
                                           gset_sample_arr[:,d],
                                           axis=0)

            time_print('Done creating features')

            return onto_sample

        elif prop=='or':

            onto_sample = self.get_features(gset_sample, prop='genes', format=format)

            if format=='dict':
                # Set all values to 1
                onto_sample = [{k: 1 for k, v in d.items()} for d in onto_sample]
            elif format in ['scipy.coo', 'scipy.csr', 'dense']:
                # Do this to preserve data type of onto_sample,
                # instead of onto_sample = (onto_sample >= 1), which
                # would create a boolean
                onto_sample[ onto_sample >= 1 ] = 1
            else:
                raise Exception('Format %s not supported for prop' % (format, prop))

            return onto_sample

        elif prop=='children':

            import igraph, numpy, collections

            graph = igraph.Graph.Read_Ncol(self.ontology_prefix)
            g_vs_index = {b:a for a,b in enumerate(graph.vs['name'])}                
            gene_2_term = self.gene_2_term
            
            # For each gset, create a dictionary mapping a perturbed term
            # (i.e. its index in self.terms) to the number of genes
            # perturbed
            genes_features = [dict(collections.Counter([term for g in gset for term in gene_2_term[g]])) for gset in gset_sample]      
            
            for term_count in genes_features:
                parent_perturbations = [self.terms_index[x] for \
                                        x in graph.vs[[p for term in term_count.keys() for \
                                                       p in graph.neighbors(g_vs_index[self.terms[term]], mode='out')]]['name']]
                parent_perturbations = dict(Counter(parent_perturbations))                    
                term_count.update(parent_perturbations)
            perturbations = genes_features

            return perturbations

        elif prop=='min_cut':
            # TODO
            raise Exception('Not supported')

        elif prop=='gene_identity':
            tmp = [(i, self.genes_index[g]) for i, gset in enumerate(gset_sample) for g in gset]
            i, j = zip(*tmp)            
            data = np.ones(len(i), dtype=dtype)
            onto_sample = scipy.sparse.coo_matrix((data, (np.array(i), np.array(j))), (len(gset_sample), len(self.genes)))

            if format=='scipy.csr':
                onto_sample = scipy.sparse.csr_matrix(onto_sample)

            assert format in ['scipy.coo', 'scipy.csr']

            return onto_sample

        elif prop=='matrix_mult':
            assert isinstance(gset_sample, (np.ndarray, pd.DataFrame)) or scipy.sparse.issparse(gset_sample), 'gset_sample must be a sample-by-gene matrix'

            if isinstance(gset_sample, pd.DataFrame):
                gene_ordering = gset_sample.columns
                gset_sample = np.array(gset_sample)

            gset_sample = scipy.sparse.coo_matrix(gset_sample)
            annotation_matrix = self.get_annotation_matrix()

            if gene_ordering is not None:
                contained = np.array([self.genes_index.has_key(g) for g in gene_ordering])
                gset_sample = scipy.sparse.coo_matrix(scipy.sparse.csc_matrix(gset_sample)[:,contained])
                subset = np.array([self.genes_index[g] for g in gene_ordering if self.genes_index.has_key(g)])
                annotation_matrix = scipy.sparse.csc_matrix(scipy.sparse.csr_matrix(annotation_matrix)[subset,:])

            onto_sample = gset_sample * annotation_matrix

            if format==pd.DataFrame:
                # Create dataframe
                unused_terms = (onto_sample!=0).toarray().sum(0) == 0
                print 'Removing %s terms with no mutations among samples' % unused_terms.sum()
                onto_sample = pd.DataFrame(onto_sample.toarray()[:,~unused_terms], columns=np.array(self.terms)[~unused_terms])

            return onto_sample

        else:
            raise Exception('Invalid perturbation type')

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

        summary = ('%s genes, '
                    '%s terms, '
                    '%s gene-term relations, '
                    '%s term-term relations' %
                    (len(self.genes),
                     len(self.terms),
                     sum([len(x) for x in self.gene_2_term.values()]),
                     sum([len(x) for x in self.parent_2_child.values()])))
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

        Parameters
        ----------

        graph
           3-column pandas.DataFrame (long format) or square pandas.DataFrame (wide format)

        alpha
           CLIXO alpha parameter

        beta
           CLIXO beta parameter

        min_dt
           Minimum similarity score

        timeout
           Maximum time (in seconds) allowed to run CLIXO

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
        
        ont = cls.from_table(output, verbose=verbose)
        ont.rename(term_prefix='CLIXO:', inplace=True)

        if verbose:
            print 'Ontology:', ont.summary()

        return ont

    def to_NdexGraph(self,
                     name=None,
                     description=None,
                     term_2_uuid=False,
                     layout='bubble',
                     alignment=None,
                     represents=False,
                     node_attr=None,
                     edge_attr=None,
                     gene_prefix='',
                     term_prefix=''):
        """Formats an Ontology object into a NetworkX object with extra node
        attributes that are accessed by the hierarchical viewer.

        """

        # Convert to NetworkX
        G = self.to_networkx(node_attr=node_attr,
                             edge_attr=edge_attr,
                             layout=layout,
                             spanning_tree=True)

        set_label = (node_attr is None) or ('Label' not in node_attr.columns)
        
        for t in self.terms:        
            if set_label:
                G.node[t]['Label'] = '%s%s' % (term_prefix, t)
            if represents:
                G.node[t]['represents'] = t
            if term_2_uuid:
                G.node[t]['ndex:internalLink'] = '[%s](%s)' % (G.node[t]['Label'], term_2_uuid[t])
        for g in self.genes:
            if set_label:
                G.node[g]['Label'] = '%s%s' % (gene_prefix, g)

        if alignment is not None:
            update_nx_with_alignment(G, alignment, use_node_name=False)

        G = nx_to_NdexGraph(G)
        if name is not None:
            G.set_name(name)
        if description is not None:
            G.set_network_attribute('Description', description)

        try:
            for x, n in ont_ndex.nodes(data=True):
                n['x_pos'] = float(n['x_pos'])
                n['y_pos'] = float(n['y_pos'])
        except:
            pass
    
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
                            sim,
                            features,
                            name,
                            ndex_server=ddot.config.ndex_server,
                            ndex_user=ddot.config.ndex_user,
                            ndex_pass=ddot.config.ndex_pass,
                            gene_columns=['Gene1', 'Gene2'],
                            propagate=True,
                            public=False,
                            node_attr=None,
                            verbose=False):
        """Push subnetworks to NDEx"""

        import ndex.client as nc
        from ndex.networkn import NdexGraph

        if propagate:
            ontology = self.copy()
            ontology.propagate_annotations()
        else:
            ontology = self

        ndex = nc.Ndex(ndex_server, ndex_user, ndex_pass)
        term_2_uuid = {}

        start = time.time()
        g1, g2 = gene_columns[0] + '_lex', gene_columns[1] + '_lex'

        if verbose:
            print 'features:', features
            print 'gene_columns:', gene_columns

        sim = sim[features + gene_columns].copy()

        # Filter dataframe for gene pairs within the ontology
        genes_set = set(ontology.genes)
        tmp = [x in genes_set and y in genes_set
               for x, y in zip(sim[gene_columns[0]], sim[gene_columns[1]])]
        sim = sim.loc[tmp, :]

        for feat in features:
            sim[feat] = sim[feat].astype(np.float64)

        # Lexicographically sort gene1 and gene2 so that gene1 < gene2
        sim[g1] = sim[gene_columns].min(axis=1)
        sim[g2] = sim[gene_columns].max(axis=1)
        sim_idx = {x : i for i, x in enumerate(zip(sim[g1], sim[g2]))}

        if verbose:
            print 'Setup time:', time.time() - start

        # Normalize features into z-scores
        tmp = sim[features]
        sim[features] = (tmp - tmp.mean()) / tmp.std()

        # Calculate the min/max range of features
        feature_mins = sim[features].min().astype(np.str)
        feature_maxs = sim[features].max().astype(np.str)

        for t in ontology.terms:
            start = time.time()

            genes = [ontology.genes[g] for g in ontology.term_2_genes[t]]
            genes.sort()
            gene_pairs_idx = [sim_idx[gp] for gp in itertools.combinations(genes, 2) \
                              if sim_idx.has_key(gp)]

            G_nx = nx.from_pandas_dataframe(sim.iloc[gene_pairs_idx, :], g1, g2,
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
        """
        Compute common ancestor matrix.

        lca[a,b] = index of least common ancestor of terms a and b
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

    def collapse_node(self,
                      g,
                      v,
                      edge_filter=None,
                      use_v_name=False,
                      combine_attrs=None,
                      default_attr=None,
                      verbose=True,
                      fast_collapse=False,
                      delete=True):

        if use_v_name:
            assert isinstance(v, (unicode, str))
            v = g.vs.find(name_eq=v).index

        if fast_collapse:
            parents = g.neighbors(v, mode='out')
            children = g.neighbors(v, mode='in')

            if len(parents) > 0 and len(children) > 0:
                # A faster collapse that adds all new edges simultaneously. Ignores edge attributes        
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
