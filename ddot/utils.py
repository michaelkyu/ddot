import numpy as np
import sys
from datetime import datetime
import pandas as pd
import networkx as nx
import simplejson
import base64
import time

def make_index(it):
    """Create a dictionary mapping elements of an iterable to the index
    position of that element

    """
    return {b : a for a, b in enumerate(it)}

def time_print(*s):
    print ' '.join(map(str, s)), datetime.today()
    sys.stdout.flush()

def pivot_square(df, index, columns, values, fill_value=0):
    """Convert a dataframe into a square compact representation"""

    df = df.pivot(index=index, columns=columns, values=values)    
    index = df.index.union(df.columns)
    df = df.reindex(index=index, columns=index, copy=False)

    # Make square matrix symmetric
    tmp = df.values.copy()
    tmp[np.isnan(tmp)] = np.inf
    tmp = np.minimum(tmp, tmp.T)
    tmp[np.isinf(tmp)] = fill_value

    df.iloc[:,:] = tmp

    return df

def melt_square(df, columns=['Gene1', 'Gene2'], similarity='similarity', empty_value=0, upper_triangle=True):
    """Melts square dataframe into sparse representation"""

    assert df.shape[0] == df.shape[1]
    
    if upper_triangle:
        assert np.all(df.index.values == df.columns.values)
        df = df.copy()
        tmp = df.values
        tmp[np.tril_indices(tmp.shape[0], k=0)] = np.nan
        df.iloc[:,:] = tmp
    
    tmp = df.stack()
    tmp.dropna(inplace=True)
    tmp.index.rename(columns, inplace=True)
    tmp.rename(similarity, inplace=True)
    return tmp.reset_index()

def get_gene_name_converter(genes, scopes='symbol', fields='entrezgene', species='human'):
    """Query mygene.info to get a dictionary mapping gene names in the ID                                                                                                    
    namespace scopes to the ID namespace in fields                                                                                                                           

    Weird behavior with mygene.info: for Entrez genes, use fields 'entrezgene'. For ENSEMBL genes, use fields "ensembl"
    """

    if hasattr(genes, '__iter__') and not isinstance(genes, (str, unicode)):
        genes = ','.join(genes)
        
    import requests
    r = requests.post('http://mygene.info/v3/query',
                      data={'q': genes,
                            'scopes': scopes,
                            'fields': fields,
                            'species': species})
    
    def parse_field(x):
        if isinstance(x, dict):
            return [unicode(x['gene'])]
        elif isinstance(x, list):
            return [unicode(y['gene']) for y in x]
        else:
            return unicode(x)
    
    dic = {x['query'] : parse_field(x[fields]) for x in r.json() if x.has_key(fields)}
    return dic

def update_nx_with_alignment(G,
                             alignment,
                             term_descriptions=None,
                             use_node_name=True):
    """Add node attributes to a NetworkX graph. 

    Parameters 
    ----------
    G
      NetworkX object

    alignment
       pandas.DataFrame where the index is the name of terms,
       and where there are 3 columns: 'Term', 'Similarity', 'FDR'

    use_node_name : bool

    term_descriptions : dict

    """

    alignment = alignment.copy()
    alignment['FDR'] = alignment['FDR'].astype(str)
    alignment['Similarity'] = alignment['Similarity'].astype(str)

    for node_idx, node_attr in G.nodes(data=True):
        if use_node_name:
            node_name = node_attr['name']
        else:
            node_name = node_idx

        if node_attr['Gene_or_Term']=='Term' and node_name in alignment.index:
            row = alignment.loc[node_name, :]
            if term_descriptions is not None:
                descr = term_descriptions[row['Term']]
                node_attr['Label'] = '%s\n%s' % (node_name, descr)
                node_attr['Aligned_Term_Description'] = descr

            node_attr['represents'] = row['Term']
            node_attr['Aligned_Term'] = row['Term']
            node_attr['Aligned_Similarity'] = row['Similarity']
            node_attr['Aligned_FDR'] = row['FDR']

    return G
    
###################################################
# NetworkX, NdexGraph, and NDEx format converters #
###################################################

def set_node_attributes_from_pandas(G, node_attr):
    G_nodes = set(G.nodes())
    node_attr = node_attr.loc[[x for x in node_attr.index if x in G_nodes], :]
    if node_attr is not None:
        for feat in node_attr.columns:
            nx.set_node_attributes(G, feat, node_attr[feat].to_dict())

def nx_nodes_to_pandas(G, attr_list=None):
    """Create pandas.DataFrame of node attributes of a NetworkX graph"""

    if attr_list is None:
        attr_list = list(set([a for d in G.nodes(data=True)
                              for a in d[1].keys()]))
    return pd.concat([pd.Series(nx.get_node_attributes(G,a), name=a) for a in attr_list],
                     axis=1)

def nx_edges_to_pandas(G, attr_list=None):
    """Create pandas.DataFrame of edge attributes of a NetworkX graph"""

    if attr_list is None:
        attr_list = list(set([a for d in G.edges(data=True) 
                              for a in d[2].keys()]))
    return pd.concat([pd.Series(nx.get_edge_attributes(G,a), name=a) for a in attr_list],
                     axis=1)

def nx_to_NdexGraph(G_nx):
    """Converts a NetworkX into a NdexGraph object"""

    from ndex.networkn import NdexGraph
    G = NdexGraph()
    node_id = 0
    node_dict = {}
    G.max_edge_id = 0
    for node_name, node_attr in G_nx.nodes_iter(data=True):
        if node_attr.has_key('name'):
            G.add_node(node_id, {k : str(v) for k, v in node_attr.items()})
        else:
            G.add_node(node_id, {k : str(v) for k, v in node_attr.items()}, name=node_name)
        node_dict[node_name] = node_id
        node_id += 1
    for s, t, edge_attr in G_nx.edges_iter(data=True):
        G.add_edge(node_dict[s], node_dict[t], G.max_edge_id, edge_attr)
        G.max_edge_id += 1

    if hasattr(G_nx, 'pos'):
        G.pos = {node_dict[a] : b for a, b in G_nx.pos.items()}
        G.subnetwork_id = 1
        G.view_id = 1

    return G

def NdexGraph_to_nx(G):
    return nx.DiGraph(nx.relabel_nodes(G, nx.get_node_attributes(G, 'name'), copy=True))

def parse_ndex_uuid(ndex_url):
    # print 'ndex_url:', ndex_url
    # print 'ndex uuid:', ndex_url.split('v2/network/')[1]
    return ndex_url.split('v2/network/')[1]

def create_edgeMatrix(X, X_cols, X_rows, verbose=True, G=None):
    import ndex.client as nc
    from ndex.networkn import NdexGraph

    if not X.flags['C_CONTIGUOUS']:
        X = np.ascontiguousarray(X)
    
    # Use base64 encoding of binary to text. More efficient than
    # pickle(*, protocol=0)
    start = time.time()
    serialized = base64.b64encode(X)
    if verbose:
        print 'base64 encoding time (sec):', time.time() - start
        print 'Size of numpy array (MB):', X.nbytes / 1e6
        print 'Size of serialization (MB):', sys.getsizeof(serialized) / 1e6
        print 'Constant factor overhead:', float(sys.getsizeof(serialized)) / X.nbytes

    if G is None:
        G = NdexGraph()
    G.unclassified_cx.append(
        {'matrix': serialized,
         'matrix_cols' : X_cols,
         'matrix_rows' : X_rows,
         'matrix_dtype' : X.dtype.name})
    
    return G

def load_edgeMatrix(ndex_uuid,
                    ndex_server,
                    ndex_user,
                    ndex_pass,
                    ndex=None,
                    json=simplejson,
                    verbose=True):
    import ndex.client as nc
    from ndex.networkn import NdexGraph

    if ndex is None:
        ndex = nc.Ndex(ndex_server, ndex_user, ndex_pass)

    start = time.time()
    response = ndex.get_network_as_cx_stream(ndex_uuid)
    if verbose:
        print 'NDEx download time (sec):', time.time() - start

    start = time.time()
    # cx = json.loads(response.text)
    cx = json.loads(response.content)
    if verbose:
        print('Read HTTP response as JSON',
              json.__name__, 'time (sec):',
              time.time() - start)

    start_loop = time.time()

    for aspect in cx:
        if 'matrix' in aspect:
            assert 'matrix_dtype' in aspect
            assert 'matrix_cols' in aspect
            assert 'matrix_rows' in aspect

            # Convert text back into binary data
            start = time.time()
            binary_data = base64.decodestring(aspect.get('matrix'))
            if verbose:
                print 'base64 decoding time (sec):', time.time() - start

            dtype = np.dtype(aspect.get('matrix_dtype'))
            rows = aspect.get('matrix_rows')
            cols = aspect.get('matrix_cols')
            dim = (len(rows), len(cols))

            # Create a NumPy array, which is nothing but a glorified
            # pointer in C to the binary data in RAM
            X = np.frombuffer(binary_data, dtype=dtype).reshape(dim)

    if verbose:
        print 'loop time (sec):', time.time() - start_loop
    
    return X, rows, cols

def sim_matrix_to_NdexGraph(sim, names, similarity, output_fmt, node_attr=None):
    """Convert similarity matrix into NdexGraph object

    Parameters
    -----------
    sim
        Square-shaped NumPy array representing similarities
    names
        Genes names, in the same order as the rows and columns of sim
    similarity
        Edge attribute name for similarities in the resulting NdexGraph object
    output_fmt
        Either 'cx' (Standard CX format), or 'cx_matrix' (custom edgeMatrix aspect)
    node_attr
        Node attributes, as a pandas.DataFrame, to be set in NdexGraph object

    Returns
    ---------
    NdexGraph object

    """

    if output_fmt == 'cx_matrix':
        # G = nx.Digraph()
        # if node_attr is not None:
        #     set_node_attributes_from_pandas(G, gene_attr)

        names = list(names)
        return create_edgeMatrix(sim, names, names)

    elif output_fmt == 'cx':
        # Keep only upper-right triangle
        sim[np.tril_indices(sim.shape[0],k=0)] = 0

        # Create NetworkX graph
        nnz = sim.nonzero()
        ebunch = [(a,b,float(c)) for (a,b),c in zip(zip(*nnz), sim[nnz])]
        G = nx.DiGraph()
        G.add_weighted_edges_from(ebunch, weight=similarity)
        nx.relabel_nodes(G, dict(enumerate(names)), copy=False)

        if node_attr is not None:
            set_node_attributes_from_pandas(G, node_attr)

        return nx_to_NdexGraph(G)
    else:
        raise Exception('Unsupported output_fmt: %s' % output_fmt)

def ndex_to_sim_matrix(ndex_uuid,
                       ndex_server,
                       ndex_user,
                       ndex_pass,
                       similarity='similarity',
                       input_fmt='cx_matrix',
                       output_fmt='matrix',
                       subset=None):
    """Read a similarity network from NDEx and return it as either a
    square np.array (compact representation) or a pandas.DataFrame of
    the non-zero similarity values (sparse representation)

    """

    from ndex.networkn import NdexGraph
    
    if input_fmt=='cx':
        # Read graph using NDEx client
        G = NdexGraph_to_nx(
              NdexGraph(
                  server=ndex_server, 
                  username=ndex_user,
                  password=ndex_pass,
                  uuid=ndex_uuid))

        # Create a DataFrame of similarity scores
        G_df = nx_edges_to_pandas(G)
        G_df.index.rename(['Gene1', 'Gene2'], inplace=True)
        G_df.reset_index(inplace=True)
        G_df[similarity] = G_df[similarity].astype(np.float64)

        nodes_attr = nx_nodes_to_pandas(G)

        if output_fmt=='matrix':
            G_sq = pivot_square(G_df, 'Gene1', 'Gene2', similarity)
            return G_sq.values, G_sq.index.values
        elif output_fmt=='sparse':
            return G_df, nodes_attr
        else:
            raise Exception('Unsupported output_fmt: %s' % output_fmt)

    elif input_fmt=='cx_matrix':
        sim, sim_names, sim_names_col = load_edgeMatrix(
            ndex_uuid,
            ndex_server,
            ndex_user,
            ndex_pass)
        assert sim_names == sim_names_col

        if subset is not None:
            idx = make_index(sim_names)
            idx = np.array([idx[g] for g in subset if idx.has_key(g)])
            sim = sim[idx, :][:, idx]
            sim_names = sim_names[idx]

        if output_fmt=='matrix':
            return sim, sim_names
        elif output_fmt=='sparse':
            G_sq = pd.DataFrame(sim, index=sim_names, columns=sim_names)
            G_df = melt_square(G_sq)            
            return G_df, None
        else:
            raise Exception('Unsupported output_fmt: %s' % output_fmt)
    else:
        raise Exception('Unsupported input_fmt: %s' % input_fmt)



def expand_seed(seed,
                sim,
                sim_names,
                agg='mean',
                min_sim=-np.inf,
                filter_perc=None,
                seed_perc=None,
                agg_perc=0.5,
                expand_size=None,
                include_seed=True,
                figure=False):
    """Take genes that are most similar to the seed set of genes."""

    assert sim.shape[0] == sim.shape[1] and sim.shape[0] == len(sim_names)
    sim_names = np.array(sim_names)
    
    index = make_index(sim_names)
    seed_idx = np.array([index[g] for g in seed])
    non_seed_idx = np.setdiff1d(np.arange(len(sim_names)), seed_idx)

    # Calculate a similarity score between each gene and the seed
    # set of genes
    if agg=='mean':
        # Average similarity to the seed set
        sim_2_seed = sim[seed_idx, :].mean(0)
    elif agg=='min':
        # The minimum similarity to any gene in the seed set
        sim_2_seed = sim[seed_idx, :].min(0)
    elif agg=='max':
        # The maximum similarity to any gene in the seed set
        sim_2_seed = sim[seed_idx, :].max(0)
    elif agg=='perc':
        # The <agg_perc> percentile of similarities to the seed set.
        # For example, if a gene has similarities of (0, 0.2, 0.4,
        # 0.6, 0.8) to five seed genes, then the 10% similarity is 0.2
        sim_2_seed = np.percentile(100 * agg_perc, sim[seed_idx, :], axis=0)
    else:
        raise Exception('Unsupported aggregation method: %s' % agg)

    # Temporarily remove seed from consideration (they'll be reinserted later), and decrease expand_size
    if include_seed:
        expand_idx = np.setdiff1d(np.arange(len(sim_names)), seed_idx)
        expand_idx = expand_idx[np.argsort(-1 * sim_2_seed[expand_idx])]
        if expand_size is not None:
            expand_size = max(0, expand_size - seed_idx.size)
    else:
        expand_idx = np.argsort(-1 * sim_2_seed)

    if expand_idx.size > 0:
        # Maximum limit on the number of returned genes
        if expand_size is not None:
            expand_idx = expand_idx[ : expand_size]

        # Filter based on a percentile of similarities between all genes and the seed set
        if filter_perc is not None:
            min_sim = max(min_sim, np.percentile(sim_2_seed, 100 * filter_perc))

        # Filter based on a percentile of similarities between seed set to itself
        if seed_perc is not None:
            min_sim = max(min_sim, np.percentile(sim_2_seed[seed_idx], 100 * seed_perc))

#         print 'sim_2_seed[expand_idx]:', sim_2_seed[expand_idx[:40]]
#         print 'expand_idx:', expand_idx[:40]
#         print 'expand_idx size:', expand_idx.size
        print 'min_sim:', min_sim

        expand_idx = expand_idx[sim_2_seed[expand_idx] >= min_sim]
        
    expand = np.array(sim_names)[expand_idx]
    expand_sim = sim_2_seed[expand_idx]

    # Include seed genes that didn't meet the criteria
    if include_seed:
        attach_idx = np.setdiff1d(seed_idx, expand_idx)
        if attach_idx.size > 0:
            expand = np.append(expand, sim_names[attach_idx])
            expand_idx = np.append(expand_idx, attach_idx)
            expand_sim = np.append(expand_sim, sim_2_seed[attach_idx])

    if figure:
        import matplotlib.pyplot as plt
        import seaborn as sns
        fig = plt.figure()
        ax = plt.gca()
        sns.distplot(sim_2_seed[non_seed_idx],
                     ax=ax,
                     axlabel='Similarity to seed set', label='Probability Density')
        sns.distplot(sim_2_seed[seed_idx],
                     ax=ax,
                     axlabel='Similarity to seed set', label='Probability Density')
    else:
        fig = None

    return expand, expand_idx, sim_2_seed, fig

