import itertools, multiprocessing, logging, os, collections, random, math, sys, time, scipy, scipy.sparse, igraph
import numpy as np, pandas as pd
from itertools import groupby, combinations
from operator import *
from scipy.stats import hypergeom
from collections import Counter
import networkx as nx
import tempfile
from subprocess import Popen, PIPE, STDOUT

from utils import time_print, set_node_attributes_from_pandas, nx_to_NdexGraph

def read_alignment_file(f, source='Term_1'):
    df = pd.read_table(f,
                       names=['Term_1', 'Term_2', 'Similarity', 'FDR', 'Term_1_Size'],
                       dtype={'Term_1':str, 'Term_2':str, 'Similarity':np.float64, 'FDR':np.float64, 'Term_1_Size':np.int64},
                       header=None)
    target = 'Term_2' if source=='Term_1' else 'Term_1'
    df.rename(columns={target : 'Term'}, inplace=True)
    df.set_index(source, inplace=True)
    return df

def align_hierarchies(hier1, hier2,
                      iterations, threads,
                      calculateFDRs_cmd,
                      output=None):
    if output is None:
        with tempfile.NamedTemporaryFile('w', delete=True) as output_file:
            return align_hierarchies(hier1, hier2, iterations, threads,
                                     output=output_file.name,
                                     calculateFDRs_cmd=calculateFDRs_cmd)

    def to_file(hier):
        if isinstance(hier, cls):
            with tempfile.NamedTemporaryFile('w', delete=False) as f:
                hier.to_3col_table(f, parent_child=True)
            hier = f.name
        else:
            assert isinstance(hier, file) or \
                (isinstance(hier, (str, unicode)) and os.path.exists(hier))
        return hier

    hier1 = to_file(hier1)
    hier2 = to_file(hier2)

    output_dir = tempfile.mkdtemp(prefix='tmp')
    cmd = '{5} {0} {1} 0.05 criss_cross {2} {3} {4} gene'.format(
           hier1, hier2, output_dir, iterations, threads, calculateFDRs_cmd)
    time_print(cmd)

    try:
        p = Popen(cmd, shell=True)
        p.wait()

        # Five columns in the output file
        # 1) Term from first ontology (Kramer calls it the "computed" ontology)
        # 2) Term from second ontology (Kramer calls it the "reference" ontology)
        # 3) Similarity value
        # 4) FDR
        # 5) Size of the term in the first ontology
        shutil.copy(os.path.join(output_dir, 'alignments_FDR_0.1_t_0.1'), output)
    finally:
        if os.path.isdir(output_dir):
            shutil.rmtree(output_dir)

        if p.poll() is None:
            if verbose: time_print('Killing alignment process %s. Output: %s' % (p.pid, output))
            p.kill()  # Kill the process

        return read_alignment_file(output)

def read_term_descriptions(description_table, alt_id=None):

    # '/cellar/users/mikeyu/DREAM_2015/prior networks/GO/goID_2_name.tab'
    # '/cellar/users/mikeyu/DREAM_2015/prior networks/GO/goID_2_alt_id.tab'

    ## Get table mapping GO terms to term descriptions
    go_descriptions = pd.read_table(description_table, header=None, names=['Term', 'Term_Description'], index_col=0)
    if alt_id is not None:
        alt_term = pd.read_table(alt_id, header=None, names=['Term', 'Alt_Term'], index_col=0)
        tmp = go_descriptions.loc[alt_term.index, 'Term_Description'].to_frame()
        tmp.index = alt_term['Alt_Term']
        go_descriptions = go_descriptions.append(tmp)
    go_descriptions = go_descriptions['Term_Description']

    return go_descriptions

def parse_obo(obo,
              output_file=None,
              id2name_file=None,
              id2namespace_file=None,
              alt_id_file=None):

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

    import pandas as pd

    pd.DataFrame(edges).to_csv(output_file, header=False, index=False, sep='\t')
    pd.Series(id2name).to_csv(id2name_file, sep='\t')
    pd.Series(id2namespace).to_csv(id2namespace_file, sep='\t')
    pd.Series(dict([(a, c) for a, b in alt_id.items() for c in b])).to_csv(alt_id_file, sep='\t')

def parse_gaf(gaf):
    """
    Read gene-term annotations from GAF file format:

    http://geneontology.org/page/go-annotation-file-gaf-format-21

    Returns a list of 2-tuples (gene, GO term)
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
    assert df['DB'].unique().size == 1 and df['DB'].unique()[0]=='UniProtKB'    

    return df.loc[:, ['DB Object ID', 'GO ID']].values.tolist()

class Ontology:

    def __init__(self,
                 hierarchy,
                 mapping,
                 hierarchy_attr=None,
                 parent_child=False,
                 add_root_name=None,
                 propagate=False,
                 verbose=True):
        """
        hierarchy: Iterable of (child term, parent) pairs
        mapping:  Iterable of (gene, term) pairs
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

        ## Check that the set of terms is the same according to parent_2_child and self.gene_2_term
        terms_A = set.union(set(self.parent_2_child.keys()), *[set(x) for x in self.parent_2_child.values()])
        terms_B = set.union(*[set(x) for x in self.gene_2_term.values()])

        if verbose and len(terms_B - terms_A)>0:
            print 'WARNING: There are {} terms that are annotated to genes but not connected to the rest of the ontology'.format(len(terms_B - terms_A))
        if verbose and len(terms_A - terms_B)>0:
            print 'WARNING: There are {} terms that have no direct gene annotations'.format(len(terms_A - terms_B))
            
        self.terms = sorted(list(terms_A | terms_B))
        self.genes = sorted(self.gene_2_term.keys())
        
        ## terms_index[<term_name>] --> index in self.terms
        self.terms_index = {b:a for a,b in enumerate(self.terms)}

        ## self.genes_index[<gene_name>] --> index in self.genes        
        self.genes_index = {b:a for a,b in enumerate(self.genes)}

        ## Convert self.gene_2_term to list term indices rather than term names
        for k, v in self.gene_2_term.items():
            self.gene_2_term[k] = [self.terms_index[x] for x in self.gene_2_term[k]]

        self.hierarchy_attr = hierarchy_attr
        self._update_fields()

        if propagate:
            self.propagate_annotations(inplace=True)
            self._update_fields()

    def _update_fields(self):
        self.child_2_parent, self.child_2_parent_indices = self.get_child_2_parent()
        self.term_2_genes = self.get_term_2_genes()
        self.term_sizes = self.get_term_sizes()

    def get_child_2_parent(self):
        ## Read term-to-term edges        
        # child_2_parent[<term_name>] --> list of <term_name>'s parent term names
        # child_2_parent_indices[<term_name>] --> list of indices of <term_name>'s parent terms
        child_2_parent = {r: tuple(p[1] for p in q) for r, q in \
                           itertools.groupby(sorted([(c,p) for p, c_list in self.parent_2_child.items() for c in c_list],
                                                    key=lambda a:a[0]),
                                             key=lambda a:a[0])}
        child_2_parent_indices = {r : [self.terms_index[x] for x in v] for r, v in child_2_parent.items()}
        return child_2_parent, child_2_parent_indices

    def get_roots(self):
        return sorted(set(self.parent_2_child.keys()) - set([y for x in self.parent_2_child.values() for y in x]))

    def to_networkx(self, gene_attr=None):
        import networkx as nx
        G = nx.DiGraph()

        G.add_nodes_from([(g, dict(name=g)) for g in self.genes])
        G.add_nodes_from([(t, dict(name=t)) for t in self.terms])
        G.add_edges_from([(g, self.terms[t],
                           dict(EdgeType='Gene-Term')) \
                          for g in self.genes for t in self.gene_2_term[g]])

        G.add_edges_from([(c, p,
                           dict(EdgeType='Child-Parent')) \
                          for p in self.terms for c in self.parent_2_child.get(p, [])])

        term_sizes = dict(zip(self.terms, self.get_term_sizes(propagate=True)))

        for t in self.terms:
            G.node[t]['Gene_or_Term'] = 'Term'
            G.node[t]['Size'] = term_sizes[t]
            G.node[t]['isRoot'] = False
        for g in self.genes:
            G.node[g]['Size'] = 1
            G.node[g]['Gene_or_Term'] = 'Gene'
            G.node[g]['isRoot'] = False

        if gene_attr is not None:
            set_node_attributes_from_pandas(G, gene_attr)

        # Identify the root
        root = self.get_roots()[0]
        G.node[root]['isRoot'] = True

        return G

    @classmethod
    def from_3col_table(cls, table_file, is_mapping=lambda x: x[2]=='gene', parent_child=True):
        
        df = pd.read_table(table_file, comment='#', header=None)
        tmp = df.apply(is_mapping, axis=1)
        mapping = df.loc[tmp, :].loc[:,:2].values.tolist()
        tmp2 = df.loc[~ tmp, :]
        if tmp2.shape[1] > 2:
            hierarchy, hierarchy_attr = tmp2.loc[:,:2].values.tolist(), tmp2.loc[:,2:]
        else:
            hierarchy, hierarchy_attr = tmp2.values.tolist(), None
        return cls(hierarchy, mapping, parent_child=parent_child, hierarchy_attr=hierarchy_attr)        

    @classmethod
    def from_NdexGraph(cls, G, gene_term='Gene-Term Annotation'):
        """Converts a NdexGraph object to an Ontology object"""

        tr = lambda x : 'gene' if x==gene_term else 'default'
        G_edges = [(G.node[u]['name'], G.node[v]['name'], tr(attr['Relation'])) for u, v, attr in G.edges_iter(data=True)]
        return Ontology(G_edges, combined_file=True, parent_child=False)

    def collapse_ontology(self,
                          method='mhkramer',
                          collapseRedundantNodes_path=None,
                          propagate=True,
                          verbose=True,
                          default_relation='default',
                          min_term_size=2):
    
        if method=='mhkramer':
            if propagate:
                if verbose: print 'Propagating annotations'
                ont = self.propagate_annotations(inplace=False)
            else:
                ont = self

            if collapseRedundantNodes_path is None:
                collapseRedundantNodes_path = os.path.join(os.getenv('ALIGN_ONTOLOGY'), 'collapseRedundantNodes')

            with tempfile.NamedTemporaryFile('w', delete=False) as f:
                ont.to_3col_table(f, parent_child=True, encoding=None, default_relation=u'default')
            try:                
                cmd = '%s %s' % (collapseRedundantNodes_path, f.name)
                print 'collapse command:', cmd
                p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
                collapsed, err = p.communicate()
            finally:
                os.remove(f.name)

            from StringIO import StringIO
            return Ontology.from_3col_table(StringIO(collapsed))

            #ontology = [x.split('\t') for x in collapsed.splitlines()]
            #return Ontology(ontology, parent_child=True, combined_file=True)

        elif method=='python':
            ## Returns a new ontology where redundant and empty terms have been collapsed
            g = self.to_igraph()
            if verbose: print len(g.vs), 'total nodes'

            if verbose: print 'Propagating annotations'
            self.propagate_annotations()

            parity = True
            while True:
                names_2_idx = {b : a for a, b in enumerate(g.vs['name'])}
                term_hash = {names_2_idx[t] : (len(g_list), hash(tuple(g_list))) for t, g_list in self.term_2_genes.items() if names_2_idx.has_key(t)}

                if verbose: time_print('Identify nodes to collapse')
                node_order = g.topological_sorting(mode='out')

                small_terms = [v for v in node_order if term_hash[v][0]<min_term_size]

                same_as_all_parents = [v for v in node_order \
                                          if (len(g.neighbors(v, mode='out'))>0 and all(term_hash[v]==term_hash[y] for y in g.neighbors(v, mode='out')))]
                same_as_all_children = [v for v in node_order \
                                          if (len(g.neighbors(v, mode='in'))>0 and all(term_hash[v]==term_hash[y] for y in g.neighbors(v, mode='in')))]
                if verbose: time_print('%s empty terms, %s (%s) terms that are redundant with all their parents (children)' % \
                           (len(small_terms), len(same_as_all_parents), len(same_as_all_children)))

                to_delete = list(set(small_terms) | set(same_as_all_children if parity else same_as_all_parents))

                if verbose: time_print('Collapsing %s empty terms and %s terms that redundant with its %s' % \
                                       (len(small_terms),
                                        len(same_as_all_children) if parity else len(same_as_all_parents),
                                        'children' if parity else 'parents'))
                parity = not parity

                if len(to_delete)==0:
                    break
                else:
                    for v in to_delete:
                        g = self.collapse_node(g, v, use_v_name=False, verbose=False, fast_collapse=True, delete=False)
                    g.delete_vertices(to_delete)

            g.es(relation_eq=None)['relation'] = default_relation

            remaining_terms = set([self.terms_index[x] for x in g.vs['name']])
            return Ontology(g,
                            mapping_file=[(gene, self.terms[t]) for gene, t_list in self.gene_2_term.items() for t in t_list if t in remaining_terms],
                            combined_file=False, parent_child=False)
    
    @classmethod
    def mutual_collapse(cls, ont1, ont2, verbose=False):
        common_genes = set(ont1_ont.genes) & set(ont2_ont.genes)

        ont1_ont.delete_genes(set(ont1_ont.genes) - common_genes)
        ont1_collapsed = ont1_ont.collapse_ontology(method='mhkramer')        
        ont2_ont.delete_genes(set(ont2_ont.genes) - common_genes)
        ont2_collapsed = ont2_ont.collapse_ontology(method='mhkramer')

        if verbose:
            print 'ont1_collapsed:', ont1_collapsed.summary()
            print 'ont2_collapsed:', ont2_collapsed.summary()

        return ont1_collapsed, ont2_collapsed

    def delete_terms(self, terms_to_delete):
        terms_to_delete = set(terms_to_delete)
        print 'Returning new ontology with %s terms removed' % len(terms_to_delete)
        from tempfile import NamedTemporaryFile
        ontology_file = NamedTemporaryFile()
        ontology_file.write(
            '\n'.join([t1+'\t'+t2+'\t'+self.relation_dict[(t2,t1)] for t1,t2_list in self.parent_2_child.items() for t2 in t2_list if (t1 not in terms_to_delete and t2 not in terms_to_delete)]) + '\n' +
            '\n'.join([g + '\t' + self.terms[t] + '\tgene' for g, t_list in self.gene_2_term.items() for t in t_list if self.terms[t] not in terms_to_delete]) + '\n')
        ontology_file.flush()
        return Ontology(ontology_file.name, combined_file=True, parent_child=False)

    def delete_genes(self, genes_to_delete):
        genes_to_delete = set(genes_to_delete)
        self.genes = [g for g in self.genes if g not in genes_to_delete]
        self.genes_index = {b : a for a, b in enumerate(self.genes)}
        self.gene_2_term = {g : t for g, t in self.gene_2_term.items() if g not in genes_to_delete}
        self._update_fields()
        
    def convert_gene_names(self, new_names):
        """
        new_names : Dictionary mapping current gene names to new names
        
        Delete genes not in dictionary"""        

        ont = self.copy()
        ont.delete_genes([g for g in ont.genes if not new_names.has_key(g)])
        ont.genes = [new_names[g] for g in ont.genes]
        ont.gene_2_term = {new_names[g] : t for g, t in ont.gene_2_term.items()}
        return ont

    def to_pandas(self, parent_2_child=True, gene_2_term=True, default_relation=u'default'):
        df = pd.DataFrame(columns=['Parent','Child','Relation'])
        if parent_2_child:
            df = df.append(self.term_hierarchy_to_pandas(default_relation),
                           ignore_index=True)
        if gene_2_term:
            tmp = self.annotations_to_pandas()
            tmp.rename(columns={'Gene':'Child', 'Term':'Parent'},
                       inplace=True)
            tmp['Relation'] = 'gene'
            df = df.append(tmp, ignore_index=True)
        return df

    def term_hierarchy_to_pandas(self, default_relation=u'default'):
        if hasattr(self, 'relation_dict'):
            relation_dict = self.relation_dict
        else:
            relation_dict = {}

        triples = [(p,c, relation_dict.get((c, p), default_relation)) \
                   for p, c_list in self.parent_2_child.items() for c in c_list]
        df = pd.DataFrame(triples, columns=['Parent', 'Child', 'Relation'])
        return df

    def annotations_to_pandas(self):
        pairs = [(g, self.terms[t]) for g, t_list in self.gene_2_term.items() for t in t_list]
        df = pd.DataFrame(pairs, columns=['Gene', 'Term'])
        return df

    def to_3col_table(self, output, header=False,
                      parent_child=True, encoding=None, default_relation=u'default'):        
        df = self.to_pandas(default_relation=default_relation)
        if parent_child:
            df = df[['Parent','Child','Relation']]
        else:
            df = df[['Child','Parent','Relation']]
        return df.to_csv(output, header=header, index=False, sep='\t')

    def copy(self):
        hierarchy = [(c, p) for p, c_list in self.parent_2_child.items() for c in c_list]
        mapping = [(g, self.terms[t]) for g, t_list in self.gene_2_term.items() for t in t_list]
        return Ontology(hierarchy, mapping,
                        hierarchy_attr=None if (self.hierarchy_attr is None) else self.hierarchy_attr.copy(),
                        parent_child=False)

    def transitive_closure(self, g=None):

        relations = ['is_a', 'regulates', 'positively_regulates', 'negatively_regulates', 'has_part', 'part_of', 'gene']
        relations_index = {b : a for a, b in enumerate(relations)}
        go_reasoning = dict()
        for r in relations:
            go_reasoning[('is_a', r)] = r
            go_reasoning[(r, 'is_a')] = r
            if r != 'has_part':
                go_reasoning[(r, 'part_of')] = 'part_of'
        for r in ['regulates', 'positively_regulates', 'negatively_regulates']:
            go_reasoning[(r, 'part_of')] = 'regulates'
        go_reasoning[('has_part', 'has_part')] = 'has_part'

        if g is None:
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

    def semantic_similarity(self, genes_subset=None, term_sizes='subset', between_terms=False, output='Resnik'):
        """Computes the semantic similarity between pair of genes in
        <genes_subset>. Similarity s(g1,g2) is defined as
        -log_2(|T_sca| / |T_root|) where |T| is the number of genes in
        <genes_subset> that are under term T. T_sca is the "smallest
        common ancestor", the common ancestral term with the smallest
        term size. T_root is the root term of the ontology.

        genes_subset : The set of genes, over which pairs the
        similarity will be calculated. If <term_sizes>='subset', then
        term sizes will be recalculated according to only these genes,
        rather than all genes in the ontology

        between_terms : if True, then output similarity between all
        terms

        output: type of semantic similarity
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

    def get_term_2_genes(self, verbose=False): 
        if verbose: print 'Calculating term_2_genes'
        term_2_genes = {self.terms[c]: [self.genes_index[x[0]] for x in d] \
                        for c, d in itertools.groupby(sorted([(a,t) for a, terms in self.gene_2_term.items() for t in terms],
                                                             key=lambda x:x[1]),
                                                      key=lambda x:x[1])}
        for t in self.terms:
            if not term_2_genes.has_key(t):
                term_2_genes[t] = []
        return term_2_genes

    def get_term_sizes(self, propagate=False):
        "Returns an array of term sizes in the same order as self.terms"
        if propagate:
            ont = self.propagate_annotations(verbose=False, inplace=False)
        else:
            ont = self

        tmp = Counter([x for y in ont.gene_2_term.values() for x in y])
        term_sizes = [tmp[x] for x in range(len(ont.terms))]
        return term_sizes

    def get_term_namespace(self,f=os.path.join(os.getenv('HOME'),'GI/data/GO/filtered-GO/goID_2_namespace.tab'), quote=':'):
        if not hasattr(self, 'term_namespace'):
            self.term_namespace = {x.split('\t')[0].replace(quote, '_'):x.split('\t')[1] for x in open(f).read().splitlines()}
            for t in self.terms:
                if not self.term_namespace.has_key(t): self.term_namespace[t] = 'NA'

        return self.term_namespace

    def shuffle_genes(self, shuffle_genes=None, verbose=False):

        if shuffle_genes is None:
            shuffle_genes = self.genes.copy()
        else:
            # Only shuffle a subset of genes
            assert isinstance(shuffle_genes, collections.Iterable)
            shuffle_genes = list(set(self.genes) & set(shuffle_genes))

        if verbose: time_print('Shuffling %s genes' % len(shuffle_genes))
        gene_permutation = range(len(shuffle_genes))
        random.shuffle(gene_permutation)
        gene_2_term = {shuffle_genes[gene_permutation[i]] : gene_2_term[shuffle_genes[i]] for i in range(len(shuffle_genes))}

        return gene_2_term

    def get_tree_edges(self):
        """Identify a spanning tree of the DAG (including genes as part of the
        DAG), and return a list of (u, v) edges in the tree.
        """
        ont = self.propagate_annotations(inplace=False)
        graph = ont.to_igraph(include_genes=True)        
        tmp = ont.get_term_sizes()
        parent_priority = [tmp[self.terms_index[v['name']]] if self.terms_index.has_key(v['name']) else 1 for v in graph.vs]

        # Identify spanning tree
        tree = self._make_tree_igraph(graph, parent_priority=parent_priority, optim=min)
        tree_edges = set([(tree.vs[e.source]['name'], tree.vs[e.target]['name']) for e in tree.es if e['smallest_parent']])
        return tree_edges

    def is_dag(self):
        return self.to_igraph().is_dag()

    def to_igraph(self, include_genes=False):
        """The vertices in the resulting graph will be in the same order as
        self.terms
        """
        if include_genes:
            terms_index_offset = {t : v + len(self.genes) for t, v in self.terms_index.items()}
            edges = [(self.genes_index[g], t) for g in self.genes for t in self.gene_2_term[g]] + \
                    [(terms_index_offset[c], terms_index_offset[p]) for p, children in self.parent_2_child.items() for c in children]
            graph = igraph.Graph(n=len(self.genes) + len(self.terms),
                                 edges=edges,
                                 directed=True,
                                 vertex_attrs={'name':self.genes + self.terms})
        else:
            edges = [(self.terms_index[c], self.terms_index[p]) for p, children in self.parent_2_child.items() for c in children]
            graph = igraph.Graph(n=len(self.terms),
                                 edges=edges,
                                 directed=True,
                                 vertex_attrs={'name':self.terms})
        
        return graph
        
    def get_common_ancestors(self, x, y):
        """
        Returns list of indices of the common ancestors between terms with index a and b.

        Ancestors are sorted by term size.
        """
        
        d = self.get_connectivity_matrix()
        term_sizes = self.get_term_sizes()

        return sorted(np.intersect1d(d[x,:].nonzero()[0], d[y,:].nonzero()[0], assume_unique=True),
                      key=lambda a: term_sizes[a])

    def get_shortest_paths(self, sparse=False, chunk_size=500):
        """
        d[a,b] = length of the shortest path from a to b, bottom-up the ontology
        """
        
        graph = self.to_igraph()

        if sparse:
            print len(graph.vs), chunk_size
            return scipy.sparse.vstack([scipy.sparse.csr_matrix(graph.shortest_paths(graph.vs[x[0]:x[1]], graph.vs, mode='out')) \
                                        for x in split_indices_chunk(len(graph.vs), chunk_size)])
        else:
            return np.array(graph.shortest_paths(graph.vs, graph.vs, mode='out'), order='C')
        
    def get_longest_paths(self):
        """
        d[a,b] = length of the longest path from a to b, bottom-up the ontology
        """

        graph = self.to_igraph()

        return -1 * np.array(graph.shortest_paths(graph.vs, graph.vs, weights=[-1 for x in graph.es], mode='out'), order='C')

    def get_connectivity_matrix(self, sparse=False):
        """
        Creates a term-by-term matrix d where d[a,b] is 1 if term b is
        an ancestor of term a, and 0 otherwise.

        Note:
           d[a,a] == 1
           d[root,j] == 0, for every j

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
        """
        Returns a similar matrix as in Ontology.get_connectivity_matrix(),
        but the diagonal of the matrix is 0.

        Note: !!!!!!!!!!!!!!!!!!!!!!!!
            d[a, a] == 0 instead of 1

        """
        
        d = self.get_connectivity_matrix()
        
        d[np.diag_indices(d.shape[0])] = 0
        assert not np.isfortran(d)
        return d

    def get_leaves(self, terms_list, children_list=None):
        """
        Select the terms in <terms_list> that are not ancestors of any term in <children_list>.
        
        If <children_list> is None, then select the terms in <terms_list> that are not ancestors of any of the other terms in <terms_list>.

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
                              method='iterative_union',
                              verbose=True, inplace=True):
        """Propagates gene-term annotations through the ontology"""
            
        if inplace:
            ont = self
        else:
            ont = self.copy()

        if direction=='forward':
            if method=='iterative_union':

                child_2_parent_idx = {ont.terms_index[c] : [ont.terms_index[p] for p in p_list] for c, p_list in ont.child_2_parent.items()}
                gene_2_term_set = {g : set(t_list) for g, t_list in ont.gene_2_term.items()}

                genes_to_update = set(ont.gene_2_term.keys())
                count = 0
                while len(genes_to_update) > 0:
                    # Iterate over a copy of genes_to_update
                    for g in genes_to_update.copy():
                        curr_terms = gene_2_term_set[g]
                        num_old = len(curr_terms)
                        curr_terms.update(set([p for t in curr_terms for p in child_2_parent_idx.get(t, [])]))
                        if len(curr_terms) == num_old:
                            genes_to_update.remove(g)                        

                    if verbose: print count,
                    count +=1
                    if count == 1000:
                        raise Exception('ERROR: Ontology depth >1000. Stopping in case of bug in code')
                if verbose: print
                ont.gene_2_term = {g : sorted(t_set) for g, t_set in gene_2_term_set.items()}
            else:
                ancestor_matrix = np.array(ont.get_connectivity_matrix(), dtype=np.int32)
                ont.gene_2_term = {g : ancestor_matrix[t, :].sum(0).nonzero()[0].tolist() for g, t in ont.gene_2_term.items()}
            
            self._update_fields()

        elif direction=='backward':
            print 'WARNING: assumes that annotations are already forward propagated'
            # parent_2_children_idx = {ont.terms_index[p] : [ont.terms_index[c] for c in c_list] for p, c_list in ont.parent_2_child.items()}
            # gene_2_term_set = {g : set(t_list) for g, t_list in ont.gene_2_term.items()}

            # graph = ont.to_igraph()
            # for parent in graph.vs[graph.topological_sorting(mode='in')]['name']:
            #     for c in parent_2_children_idx[parent]

            # ont.gene_2_term = {g : sorted(t_set) for g, t_set in gene_2_term_set.items()}
        else:
            raise Exception('Unsupported direction')

        return ont

    def propagate_ontotypes(self, ontotypes, prop, ontotype_size, max_ontotype, method='fixed_size'):
        """
        Propagates a list of base ontotypes.
        
        @params
        ontotypes: An array of ontotypes in scipy.sparse.csr_matrix format.  Each row is a separate ontotype.
        method: If 'fixed_size', then the sum of values in every ontotype is exactly the same, namely <ontotype_size>
        ontotype_size: If method=='fixed_size', then this is the sum of values in every ontotype
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

    def get_ontotype(self, gset_sample, prop, format='dict', dtype=np.int8, gene_ordering=None):
        """ Computes the ontotype for each gset in <gset_sample>
        according to the propagation scheme <prop>.  Returns a list,
        where each element is a dictionary representation of an
        ontotype.
        
        gset_sample = list of gene sets (i.e. gsets)
        A "gset" is a set of genes, represented as a tuple

        prop == 'genes' means the value of a term is the number of
        genes in the term which are deleted.

        scipy.coo : The Scipy Sparse COO format http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html#scipy.sparse.coo_matrix

        When format=='scipy.coo', this is how you convert it to a
        dense matrix in a text file.  Let the return matrix be x.  Run
        np.savetxt(<filename.txt>, x.toarray(), sep='\t').  Look up
        toarray() function in Scipy Sparse doc page and the np.savetxt
        function.
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

    def format_features(self, features, fmt='table'):

        if fmt=='table':
            return '\n'.join(['\t'.join([str(x)+'\t'+str(y) for x,y in sorted(z.items(), key=lambda a:a[0])]) for z in features])

        elif fmt in ['scipy.coo', 'scipy.csr', 'numpy.dense']:
            time_print('Getting element indices for sparse representation')
            i = [ind for ind, feat in enumerate(features) for x in range(len(feat))]
            tmp = [y for feat in features for x in feat.items() for y in x]
            j = tmp[0:len(tmp):2]
            data = tmp[1:len(tmp):2]
            
            time_print('Creating sparse COO matrix')
            features = scipy.sparse.coo_matrix((data, (i, j)), shape=(len(features), len(self.terms)), dtype=np.int8)

            time_print('Converting to specified format')
            if fmt=='scipy.csr':
                features = scipy.sparse.csr_matrix(features)
            elif fmt=='numpy.dense':
                features = features.todense()
            
            return features        

    def get_annotation_matrix(self):
        # Convert gene names to indices
        gene_2_term = [(self.genes_index[g], t_list) for g, t_list in self.gene_2_term.items()]
        self.annotation_matrix = scipy.sparse.coo_matrix(([1 for g, t_list in gene_2_term for t in t_list],
                                                          ([g for g, t_list in gene_2_term for t in t_list],
                                                           [t for g, t_list in gene_2_term for t in t_list])),
                                                         shape=(len(self.genes), len(self.terms)))
        return self.annotation_matrix

    def summary(self):
        return '%s genes, %s terms, %s gene-term relations, %s term-term relations' % \
            (len(self.genes), len(self.terms), sum([len(x) for x in self.gene_2_term.values()]), sum([len(x) for x in self.parent_2_child.values()]))
    
    @classmethod
    def run_clixo(cls, graph, alpha, beta, dt_thresh, max_time,
                  clixo_folder, output=None, output_log=None, verbose=True):

        rerun = False

        if output is None:
            output_file = tempfile.NamedTemporaryFile('w', delete=False)
            output = output_file.name
            print 'temp output:', output
            rerun, delete_output = True, True

        if not (isinstance(graph, str) and os.path.exists(graph)):
            # Write graph into a temporary file.
            # Assumes that <graph> is a list of 3-tuples (parent, child, score)
            with tempfile.NamedTemporaryFile('w', delete=False) as graph_file:
                try:
                    graph.to_csv(graph_file, sep='\t', header=False, index=False)
                except:
                    graph_file.write('\n'.join(['\t'.join([str(x[0]),str(x[1]),str(x[2])]) for x in graph]) + '\n')
            graph = graph_file.name
            print 'temp graph:', graph
            rerun, delete_graph = True, True

        if not (isinstance(output_log, str) and os.path.exists(output_log)):
            output_log_file = tempfile.NamedTemporaryFile('w', delete=False)
            output_log = output_log_file.name
            print 'temp output log:', output_log
            rerun, delete_output_log = True, True

        if rerun:            
            try:
                return cls.run_clixo(graph, alpha, beta, dt_thresh, max_time,
                                     clixo_folder,
                                     output=output,
                                     output_log=output_log, verbose=verbose)
            finally:
                if delete_output:
                    os.remove(output)
                if delete_output_log:
                    os.remove(output_log)
                if delete_graph:
                    os.remove(graph)

        if verbose: time_print('\t'.join(map(str, [graph, alpha, beta, dt_thresh])))

        # '/cellar/users/mikeyu/mhk7-clixo_0.3-cec3674'
        clixo_cmd = os.path.join(clixo_folder, 'clixo')
        extract_cmd = os.path.join(clixo_folder, 'extractOnt')

        # For timestamping everyline: awk '{ print strftime("%Y-%m-%d %H:%M:%S"), $0; fflush(); }'
        cmd = """{0} {1} {2} {3} | awk""".format(clixo_cmd, graph, alpha, beta) + \
              """ '{if ( $1 ~ /^#/ ) {print "\#", strftime("%Y-%m-%d %H:%M:%S"), $0 ; fflush() } else {print $0}}'""" + \
              """ | tee {}""".format(output_log)
        print 'CLIXO command:', cmd

        p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, bufsize=1)

        curr_dt = None
        start = time.time()

        # Asynchronous readout of the command's output
        while p.poll() is None:

            # Break if passed the maximum processing time
            if time.time() - start > max_time:
                if verbose: time_print('Killing process %s (OUT OF TIME). Current dt: %s: Output: %s' % (p.pid, curr_dt, output_log))
                break

            line = p.stdout.readline()

            # Remove newline character
            line = line[:-1]

            # Break if the dt_threshold has been met
            if '# dt: ' in line:
                curr_dt = float(line.split('# dt: ')[1])
                if curr_dt < dt_thresh:
                    if verbose: time_print('Killing process %s (BEYOND MIN THRESHOLD). Current dt: %s, dt_thresh: %s' % (p.pid, curr_dt, dt_thresh))
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

        time_print('Elapsed time (sec): %s' % (time.time() - start))
        
        return cls.from_3col_table(output)

    def format_for_hierarchical_viewer(self, term_2_uuid, gene_attr=None):
        """Formats an Ontology object into a NetworkX object with extra node
        attributes that are accessed by the hierarchical viewer.

        """

        # Convert to NetworkX
        ontology_nx = self.to_networkx(gene_attr=gene_attr)

        for t in self.terms:
            ontology_nx.node[t]['name'] = 'CLIXO:%s' % t
            ontology_nx.node[t]['Label'] = 'CLIXO:%s' % t
            ontology_nx.node[t]['ndex:internalLink'] = term_2_uuid[t]
        for g in self.genes:
            ontology_nx.node[g]['name'] = g
            ontology_nx.node[g]['Label'] = g

        # Identify a spanning tree
        tree_edges = self.get_tree_edges()
        import networkx as nx
        nx.set_edge_attributes(ontology_nx,
                               'Is_Tree_Edge',
                               {(s,t) : 'Tree' if ((s,t) in tree_edges) else 'Not_Tree' \
                                for s, t in ontology_nx.edges_iter(data=False)})


        return ontology_nx

    def _force_directed_layout(self, ontology_nx):
        # Force-directed layout on only the terms
        sub_nx = ontology_nx.copy()
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

    def upload_subnetworks_2_ndex(self,
                                  sim,
                                  feature_columns,
                                  ndex_server, ndex_user, ndex_pass,
                                  name,
                                  gene_columns=['Gene1', 'Gene2'],
                                  propagate=True,
                                  public=False,
                                  gene_attr=None,
                                  verbose=False):
        """Push subnetworks"""

        import ndex.client as nc
        from ndex.networkn import NdexGraph

        if propagate:
            ontology = self.copy()
            ontology.propagate_annotations()
        else:
            ontology = self

        ndex = nc.Ndex(ndex_server, ndex_user, ndex_pass)
        term_2_url = {}

        start = time.time()
        g1, g2 = gene_columns[0] + '_lex', gene_columns[1] + '_lex'

        print sim.head()
        print feature_columns, gene_columns
        sim = sim[feature_columns + gene_columns].copy()

        # Filter dataframe for gene pairs within the ontology
        genes_set = set(ontology.genes)
        tmp = [x in genes_set and y in genes_set for x, y in zip(sim[gene_columns[0]], sim[gene_columns[1]])]
        sim = sim.loc[tmp, :]

        for feat in feature_columns:
            sim[feat] = sim[feat].astype(np.float64)

        # Lexicographically sort gene1 and gene2 so that gene1 < gene2
        sim[g1] = sim[gene_columns].min(axis=1)
        sim[g2] = sim[gene_columns].max(axis=1)
        sim_idx = {x : i for i, x in enumerate(zip(sim[g1], sim[g2]))}

        print 'Setup time:', time.time() - start

        # Normalize features into z-scores
        tmp = sim[feature_columns]
        sim[feature_columns] = (tmp - tmp.mean()) / tmp.std()
#        sim[feature_columns] = (tmp - tmp.mean()).values / tmp.std().values

        # Calculate the min/max range of features
        feature_mins = sim[feature_columns].min().astype(np.str)
        feature_maxs = sim[feature_columns].max().astype(np.str)

        for t in ontology.terms:
            start = time.time()

            genes = [ontology.genes[g] for g in ontology.term_2_genes[t]]
            genes.sort()
            gene_pairs_idx = [sim_idx[gp] for gp in itertools.combinations(genes, 2) \
                              if sim_idx.has_key(gp)]

            G_nx = nx.from_pandas_dataframe(sim.iloc[gene_pairs_idx, :], g1, g2,
                                            edge_attr=feature_columns)
            if gene_attr is not None:
                set_node_attributes_from_pandas(G_nx, gene_attr)

            G = nx_to_NdexGraph(G_nx)
            G.set_name('%s supporting network for CLIXO:%s' % (name, t))
            G.set_network_attribute('Description', '%s supporting network for CLIXO:%s' % (name, t))
            for f in feature_columns:
                G.set_network_attribute('%s min' % f, feature_mins[f])
                G.set_network_attribute('%s max' % f, feature_maxs[f])

            start_upload = time.time()
            ndex_url = G.upload_to(ndex_server, ndex_user, ndex_pass)
            term_2_url[t] = ndex_url
            upload_time = time.time() - start_upload

            if public:
                ndex_uuid = parse_ndex_uuid(ndex_url)
                ndex.make_network_public(ndex_uuid)

            if verbose:
                print(ontology.terms_index[t],
                      'Term:', t,
                      'Gene pairs:', len(gene_pairs_idx),
                      'Genes:', len(genes),
                      'Time:', round(time.time() - start, 4),
                      'Upload time:', round(upload_time, 4))

        return term_2_url

    def get_best_ancestor_matrix(self, graph=None, node_order=None, verbose=False):
        """
        Compute common ancestor matrix.

        lca[a,b] = index of least common ancestor of terms a and b
        """        
        if graph is None:
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

    def _make_tree_igraph(self,
                          graph=None,
                          method='priority',
                          edge_name='smallest_parent',
                          parent_priority=None, edge_priority=None, default_priority=None, optim='max'):
        """Returns copy of graph with new edge attribute marking spanning tree"""

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

    def collapse_node(self, g, v, edge_filter=None, use_v_name=False, combine_attrs=None, default_attr=None, verbose=True, fast_collapse=False, delete=True):

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
                                print 'Setting', key, g.vs[in_neigh]['name'], g.vs[out_neigh]['name'], 'to', combine_attrs[key](e_in, e_out, e), (e_in[key], e_out[key])

                    e['collapsed_length'] = e_in['collapsed_length'] + e_out['collapsed_length']
                    e['collapsed_terms'] = e_in['collapsed_terms'] + [g.vs[v]['name']] + e_out['collapsed_terms']

        if delete:
            g.delete_vertices(v)

        return g

    def map_gene_pairs_2_triangles(self, gsets):
        """Assign each gene pair to it's smallest common ancestor.  This is
        it's 'triangle'.  It also works out that this is both of their lowest
        terms, such that this is like a within-term interaction"""

        ontology = self

        graph_orig = ontology.to_igraph()
        graph = ontology.to_igraph().copy()
        graph.add_vertices(ontology.genes)
        assert graph.vs[-len(ontology.genes):]['name'] == ontology.genes
        assert graph.vs[:len(ontology.terms)]['name'] == ontology.terms
        graph.add_edges([(g, t) for g, terms in [(len(ontology.terms) + ontology.genes_index[g], terms) for g, terms in ontology.gene_2_term.items()] for t in terms])
        gene_sca_matrix = get_smallest_ancestor(graph, ontology.get_term_sizes().tolist() + [1 for i in ontology.genes], verbose=True)
        term_2_genes = {ontology.terms_index[t] : set(g) for t, g in ontology.get_term_2_genes().items()}

        # tmp = [(i, j, gene_sca_matrix[len(ontology.terms) + i, len(ontology.terms) + j], graph_orig.neighbors(gene_sca_matrix[len(ontology.terms) + i, len(ontology.terms) + j], mode='in')) for i, j in [(ontology.genes_index[g1], ontology.genes_index[g2]) for g1, g2 in gsets]]
        # triangles = [(sca, [t for t in terms if i in term_2_genes[t]], [t for t in terms if j in term_2_genes[t]]) for i, j, sca, terms in tmp]

        #--- Alternate triangle calculation.  Allow for multiple SCA in anti-chain paths
        # Make sure the order is 'C'.  Set the diagonal to True (even though it should technically be False) for easier use of np.all below
        connectivity = ontology.get_connectivity_matrix()
        not_connectivity = np.logical_not(np.array(connectivity, dtype=np.bool, order='C'))
        np.fill_diagonal(not_connectivity, True)

        children_dict = {v.index : np.array(graph_orig.neighbors(v.index, mode='in')) for v in graph_orig.vs}
        gene_2_term_set = {g : np.array(terms) for g, terms in ontology.gene_2_term.items()}
        term_2_genes = {ontology.terms_index[t] : set(g) for t, g in ontology.get_term_2_genes().items()}

        time_print('Calculating triangles')
        triangles = ((g1, g2, gene_2_term_set[g1], gene_2_term_set[g2]) for g1, g2 in gsets)
        triangles = ((g1, g2, g1_terms, g2_terms, np.intersect1d(g1_terms, g2_terms, assume_unique=True)) for g1, g2, g1_terms, g2_terms in triangles)
        # This could be potentially double-counting a term pair if for a single gene pair, they have more than one sca
        # triangles = [[(sca, np.intersect1d(g1_terms, children_dict[sca], assume_unique=True), np.intersect1d(g2_terms, children_dict[sca], assume_unique=True)) for sca in terms[not_connectivity[terms, :][:, terms].all(axis=0)]] for g1, g2, g1_terms, g2_terms, terms in triangles]
        triangles = [[(sca,
                       tuple(np.intersect1d(g1_terms, children_dict[sca], assume_unique=True)),
                       tuple(np.intersect1d(g2_terms, children_dict[sca], assume_unique=True))) \
                          for sca in terms[not_connectivity[terms, :][:, terms].all(axis=0)]] \
                     for g1, g2, g1_terms, g2_terms, terms in triangles]

        ## Each element corresponds to a gene pair.  Each element is a
        ## list of triplets (LCA, children terms spanned by gene 1,
        ## children terms spanned by gene 2)
        return triangles

    # def count_siblings(triangles, num_terms, term_pair_triangles=None):
    #     """Returns the number of times that a sibling occurs in a set of triangles"""

    #     triu_indices = np.triu_indices(num_terms, k=1)

    #     if term_pair_triangles is None:
    #         term_pair_triangles = [set([ (t1, t2) if t1 < t2 else (t2, t1) for sca, side1, side2 in path for t1 in side1 for t2 in side2]) \
    #                                for path in triangles]

    #     return scipy.sparse.coo_matrix(([1 for term_pairs in term_pair_triangles for t1, t2 in term_pairs],
    #                                     ([t1 for term_pairs in term_pair_triangles for t1, t2 in term_pairs],
    #                                      ([t2 for term_pairs in term_pair_triangles for t1, t2 in term_pairs]))),
    #                                    shape=(num_terms, num_terms)).toarray()[triu_indices]


    # def siblings_2_gene_pairs(triangles):
    #     """Returns dictionary mapping a sibling (t1, t2), where the term
    # indices are t1<t2, to the indices of gene pairs.  This is the same as
    # indices to elements in <triangles>"""

    #     return {a : [y for x, y in b] for a, b in groupby(sorted([((t1, t2) if t1<t2 else (t2, t1), i) for i, tri in enumerate(triangles) for \
    #                                                               sca, t1_list, t2_list in tri for t1 in t1_list for t2 in t2_list],
    #                                                              key=lambda a:a[0]), key=lambda a:a[0])}
    # def lcas_2_gene_pairs(triangles):
    #     """Returns dictionary mapping a term to the indices of gene pairs.
    # This is the same as indices to elements in <triangles> """

    #     return {a : [y for x, y in b] for a, b in groupby(sorted([(sca, i) for i, tri in enumerate(triangles) for \
    #                                                                 sca, t1_list, t2_list in tri],
    #                                                                key=lambda a:a[0]), key=lambda a:a[0])}

    # def explain_gene_pairs(triangles, lcas, siblings):
    #     """Given a list of enriched siblings and lcas, returns those siblings
    # and lcas that a gene pair falls under.

    #     For each element in triangles, returns a pair (x, y) where x are
    #     terms in lcas and y are term pairs in siblings.
    #     """

    #     lcas, siblings = set(lcas), set([(t1,t2) if t1<t2 else (t2, t1) for t1, t2 in siblings])
    #     return [(set([sca for sca, t1_list, t2_list in tri]) & lcas,
    #              set([(t1,t2) if t1<t2 else (t2,t1) for sca, t1_list, t2_list in tri for t1 in t1_list for t2 in t2_list])) \
    #             for tri in triangles]
