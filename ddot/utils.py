import numpy as np
import sys
from datetime import datetime
import pandas as pd
import networkx as nx
import simplejson
import base64
import time

from ddot import config

def invert_dict(dic, sort=True, keymap={}, valmap={}):
    """Inverts a dictionary of the form
    key1 : [val1, val2]
    key2 : [val1]

    to a dictionary of the form

    val1 : [key1, key2]
    val2 : [key2]

    Parameters
    -----------
    dic : dict

    Returns
    -----------
    dict
    
    """

    dic_inv = {}
    for k, v_list in dic.items():
        k = keymap.get(k, k)
        for v in v_list:
            v = valmap.get(v, v)
            if dic_inv.has_key(v):
                dic_inv[v].append(k)
            else:
                dic_inv[v] = [k]

    if sort:
        for k in dic_inv.keys():
            dic_inv[k].sort()

    return dic_inv

def transform_pos(pos, xmin=-250, xmax=250, ymin=-250, ymax=250):
    """Transforms coordinates to fit a bounding box.

    Parameters
    -----------
    pos : dict
        Dictionary mapping node names to (x,y) coordinates

    xmin : float, optional
       Minimum x-coordinate of the bounding box

    xmax : float, optional
       Maximum x-coordinate of the bounding box

    ymin : float, optional
       Minimum y-coordinate of the bounding box

    ymax : float, optional
       Maximum y-coordinate of the bounding box

    Returns
    -------
    dict:
        New dictionary with transformed coordinates
    """

    def make_transform(x_list):
        xmin_old, xmax_old = min(x_list), max(x_list)
        midpoint_old = (xmin_old + xmax_old) / 2.
        midpoint_new = (xmin + xmax) / 2.
        scale = float(xmax - xmin) / (xmax_old - xmin_old)

        def transform_x(x):
            return (x - midpoint_old) * scale + midpoint_new

        return transform_x

    transform_x = make_transform([x for x, y in pos.values()])
    transform_y = make_transform([y for x, y in pos.values()])
    pos = {g : (transform_x(x), transform_y(y)) for g, (x, y) in pos.items()}

    return pos

def bubble_layout_nx(G, xmin=-750, xmax=750, ymin=-750, ymax=750):
    """Bubble-tree Layout using the Tulip library.

    The input tree must be a graph. The layout is scaled so that it is
    fit exactly within a bounding box.
    
    Grivet, S., Auber, D., Domenger, J. P., & Melancon,
    G. (2006). Bubble tree drawing algorithm. In Computer Vision and
    Graphics (pp. 633-641). Springer Netherlands.

    Parameters
    ----------
    G : networkx.Graph
       Tree

    xmin : float, optional
       Minimum x-coordinate of the bounding box

    xmax : float, optional
       Maximum x-coordinate of the bounding box

    ymin : float, optional
       Minimum y-coordinate of the bounding box

    ymax : float, optional
       Maximum y-coordinate of the bounding box

    Returns
    -------
    dict
       Dictionary mapping nodes to 2D coordinates. pos[node_name] -> (x,y)

    """

    from tulip import tlp
#     from tulip import *
#     from tulipgui import *

    graph = tlp.newGraph()

    nodes = graph.addNodes(len(G.nodes()))
    nodes_idx = make_index(G.nodes())

    for x, y in G.edges():
        graph.addEdge(nodes[nodes_idx[x]], nodes[nodes_idx[y]])

    # apply the 'Bubble Tree' graph layout plugin from Tulip
    graph.applyLayoutAlgorithm('Bubble Tree')

    viewLayout = graph.getLayoutProperty("viewLayout")
    pos = {g : (viewLayout[i].x(), viewLayout[i].y()) for i, g in zip(nodes, G.nodes())}

    pos = transform_pos(pos, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

    return pos

def split_indices(n, k):
    try:
        tmp = iter(n)
        indices = n
    except TypeError, te:
        assert type(n)==type(int(1))
        indices = range(n)

    from math import ceil
    chunk_size = int(ceil(float(len(indices)) / k))

    return [(chunk_size * a, min(chunk_size * (a+1), len(indices))) for a in range(int(ceil(float(len(indices)) / chunk_size)))]

def split_indices_chunk(n, k):
    from math import ceil

    try:
        iter(n)
        n = len(n)
    except TypeError, te:
        assert isinstance(n, int) or isinstance(n, long)

    return [(k*i, min(k*(i+1), n)) for i in range(int(ceil(float(n) / k)))]

def make_index(it):
    """Create a dictionary mapping elements of an iterable to the index
    position of that element

    """
    return {b : a for a, b in enumerate(it)}

def time_print(*s):
    print ' '.join(map(str, s)), datetime.today()
    sys.stdout.flush()

def pivot_square(df, index, columns, values, fill_value=0):
    """Convert a dataframe into a square compact representation.

    Parameters
    ----------
    df : pandas.DataFrame
       DataFrame in long-format where every row represents one gene pair

    Returns
    -------
    df : pandas.DataFrame
       DataFrame with gene-by-gene dimensions

    """

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
    """Melts square dataframe into sparse representation.

    Parameters
    ----------
    df : pandas.DataFrame
        Square-shaped dataframe where df[i,j] is the value of edge (i,j)

    columns : iterable    
        Column names for nodes in the output dataframe

    similarity : string
        Column for edge value in the output dataframe
    
    empty_value
        Not yet supported

    upper_triangle : bool
        Only use the values in the upper-right triangle (including the diagonal) of the input square dataframe

    Returns
    -------
    pandas.DataFrame
        3-column dataframe that provides a sparse representation of
        the edges. Two of the columns indicate the node name, and the
        third column indicates the edge value

    """

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

    Weird behavior with mygene.info: for Entrez genes, use fields
    'entrezgene'. For ENSEMBL genes, use fields "ensembl"

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

    Returns
    -------
    None

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
    """Modify node attributes according to a pandas.DataFrame.

    Parameters
    ----------

    G : networkx.Graph

    node_attr : pandas.DataFrame

    """
    
    G_nodes = set(G.nodes())
    node_attr = node_attr.loc[[x for x in node_attr.index if x in G_nodes], :]
    if node_attr is not None:
        for feature_name, feature in node_attr.iteritems():
            for n, v in feature.dropna().iteritems():
                try:
                    # If v is actually a NumPy scalar type,
                    # e.g. np.float or np.int, then convert it to a
                    # fundamental Python type, e.g. int or float
                    v = v.item()
                except:
                    pass
                G.node[n][feature_name] = v

def set_edge_attributes_from_pandas(G, edge_attr):
    """Modify edge attributes according to a pandas.DataFrame.

    Parameters
    ----------

    G : networkx.Graph

    edge_attr : pandas.DataFrame

    """

    G_edges = set(G.edges())
    edge_attr = edge_attr.loc[[x for x in edge_attr.index if x in G_edges], :]
    if edge_attr is not None:
        for feature_name, feature in edge_attr.iteritems():
            for (e1,e2), v in feature.dropna().iteritems():
                try:
                    v = v.item()
                except:
                    pass
                G.edge[(e1,e2)][feature_name] = v

def nx_nodes_to_pandas(G, attr_list=None):
    """Create pandas.DataFrame of node attributes of a NetworkX graph.

    Parameters
    ----------
    G : networkx.Graph
    
    attr_list : list, optional
       Names of node attributes. Default: all node attributes

    Returns
    -------
    pandas.DataFrame
       Dataframe where index is the names of nodes and the columns are node attributes.

    """

    if attr_list is None:
        attr_list = list(set([a for d in G.nodes(data=True)
                              for a in d[1].keys()]))
    return pd.concat([pd.Series(nx.get_node_attributes(G,a), name=a) for a in attr_list],
                     axis=1)

def nx_edges_to_pandas(G, attr_list=None):
    """Create pandas.DataFrame of edge attributes of a NetworkX graph.

    Parameters
    ----------
    G : networkx.Graph
    
    attr_list : list, optional
       Names of edge attributes. Default: all edge attributes

    Returns
    -------
    pandas.DataFrame

       Dataframe where index is a MultIndex with two levels (u,v)
       referring to edges and the columns refer to edge
       attributes. For multi(di)graphs, the MultiIndex have three
       levels of the form (u, v, key).

    """

    if attr_list is None:
        attr_list = list(set([a for d in G.edges(data=True) 
                              for a in d[2].keys()]))
    if len(attr_list) > 0:
        return pd.concat([pd.Series(nx.get_edge_attributes(G,a), name=a) for a in attr_list],
                         axis=1)
    else:
        return pd.DataFrame(index=pd.MultiIndex.from_tuples(G.edges()))

def nx_to_NdexGraph(G_nx, discard_null=True):
    """Converts a NetworkX into a NdexGraph object.

    Parameters
    ----------
    G_nx : networkx.Graph

    Returns
    -------
    ndex.networkn.NdexGraph

    """

    from ndex.networkn import NdexGraph
    G = NdexGraph()
    node_id = 0
    node_dict = {}
    G.max_edge_id = 0
    for node_name, node_attr in G_nx.nodes_iter(data=True):
        if discard_null:
            node_attr = {k:v for k,v in node_attr.items() if not pd.isnull(v)}

        if node_attr.has_key('name'):
            G.add_node(node_id, node_attr)
        else:
            G.add_node(node_id, node_attr, name=node_name)
        node_dict[node_name] = node_id
        node_id += 1
    for s, t, edge_attr in G_nx.edges_iter(data=True):
        if discard_null:
            edge_attr = {k:v for k,v in edge_attr.items() if not pd.isnull(v)}
        G.add_edge(node_dict[s], node_dict[t], G.max_edge_id, edge_attr)
        G.max_edge_id += 1

    if hasattr(G_nx, 'pos'):
        G.pos = {node_dict[a] : b for a, b in G_nx.pos.items()}
        G.subnetwork_id = 1
        G.view_id = 1

    return G

def NdexGraph_to_nx(G):
    """Converts a NetworkX into a NdexGraph object.

    Parameters
    ----------
    G : ndex.networkn.NdexGraph

    Returns
    -------
    networkx.classes.DiGraph

    """

    return nx.DiGraph(nx.relabel_nodes(G, nx.get_node_attributes(G, 'name'), copy=True))

def parse_ndex_uuid(ndex_url):
    """Extracts the NDEx UUID from a URL

    Parameters
    ----------
    ndex_url : str
        URL for a network stored on NDEx

    Returns
    -------
    str
        UUID of the network

    """
    return ndex_url.split('v2/network/')[1]

def parse_ndex_server(ndex_url):
    tmp = ndex_url.split('//')
    if len(tmp) == 2:
        # e.g. 'http://dev2.ndexbio.org/v2/network/8bfa8318-55ed-11e7-a2e2-0660b7976219'
        return tmp[0] + '//' + tmp[1].split('v2/network/')[0]
    elif len(tmp) == 1:
        # e.g. 'dev2.ndexbio.org/v2/network/8bfa8318-55ed-11e7-a2e2-0660b7976219'
        return tmp[0].split('v2/network/')[0]
    elif len(tmp) == 0 or len(tmp) > 2:
        raise Exception()        

def create_edgeMatrix(X, X_cols, X_rows, verbose=True, G=None):
    """Converts an NumPy array into a NdexGraph with a special CX aspect
    called "edge_matrix". The array is serialized using base64 encoding.
    
    Parameters
    ----------
    X : np.ndarray
    
    X_cols : list
        Column names

    X_rows : list
        Row names

    Returns
    -------
    ndex.networkn.NdexGraph        

    """

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
    """Loads a NumPy array from a NdexGraph with a special CX aspect
    called "edge_matrix".
    
    Parameters
    ----------
    ndex_uuid : str
        NDEx UUID of ontology
    
    ndex_server : str
        URL of NDEx server

    ndex_user : str
        NDEx username

    ndex_pass : str
        NDEx password

    json : module
        JSON module with "loads" function

    Returns
    -------
    X : np.ndarray
    
    X_cols : list
        Column names

    X_rows : list
        Row names

    """

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
    sim : np.ndarray
        Square-shaped NumPy array representing similarities

    names : list
        Genes names, in the same order as the rows and columns of sim

    similarity : str
        Edge attribute name for similarities in the resulting NdexGraph object

    output_fmt : str
        Either 'cx' (Standard CX format), or 'cx_matrix' (custom edgeMatrix aspect)

    node_attr : pandas.DataFrame, optional
        Node attributes, as a pandas.DataFrame, to be set in NdexGraph object

    Returns
    ---------
    ndex.networkn.NdexGraph

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
                       ndex_server=config.ndex_server,
                       ndex_user=config.ndex_user,
                       ndex_pass=config.ndex_pass,
                       similarity='similarity',
                       input_fmt='cx_matrix',
                       output_fmt='matrix',
                       subset=None):
    """Read a similarity network from NDEx and return it as either a
    square np.array (compact representation) or a pandas.DataFrame of
    the non-zero similarity values (sparse representation)

    Parameters
    ----------
    ndex_uuid : str
        NDEx UUID of ontology
    
    ndex_server : str
        URL of NDEx server

    ndex_user : str
        NDEx username

    ndex_pass : str
        NDEx password

    similarity : str

        Name of the edge attribute that represents the
        similarity/weight between two nodes. If None, then the name of
        the edge attribute in the output is named 'similarity' and all
        edges are assumed to have a similarity value of 1.

    input_fmt : str
    
    output_fmt : str
        If 'matrix', return a NumPy array. If 'sparse', return a pandas.DataFrame

    subset : optional

    Returns
    --------
    np.ndarray or pandas.DataFrame

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
        G_df.index.rename(['Node1', 'Node2'], inplace=True)
        G_df.reset_index(inplace=True)
        
        if similarity is None:
            G_df['similarity'] = 1.0
        else:
            G_df[similarity] = G_df[similarity].astype(np.float64)

        nodes_attr = nx_nodes_to_pandas(G)

        if output_fmt=='matrix':
            G_sq = pivot_square(G_df, 'Node1', 'Node2', similarity)
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
    """Identify genes that are most similar to a seed set of genes.

    A gene is included in the expanded set only if it meets all of the
    specified criteria. If include_seed is True, then genes that are
    seeds will be included regardless of the criteria.

    Parameters
    ----------
    seed : list

    sim : np.ndarray

    sim_names : list of str

    agg : str or function

    min_sim : float
       Minimum similarity to the seed set.

    filter_perc : float
       

    seed_perc : float

    agg_perc : float

    expand_size : int

    include_seed : bool

    figure : bool

    Returns
    -------
    expand

    expand_idx
    
    sim_2_seed

    fig

    """

    # Check at least one seed
    assert len(seed) > 0

    # Check similarity matrix has correct dimensions
    assert sim.shape[0] == sim.shape[1] and sim.shape[0] == len(sim_names)

    sim_names = np.array(sim_names)
    
    index = make_index(sim_names)
    seed_idx = np.array([index[g] for g in seed])
    non_seed_idx = np.setdiff1d(np.arange(len(sim_names)), seed_idx)
        
    # if seed_idx.size == 1:
    #     seed_idx = np.array([seed_idx])
#        sim_slice = np.array([sim_slice])

    # Calculate a similarity score between each gene and the seed
    # set of genes
    sim_slice = sim[seed_idx, :]        
    if agg=='mean':
        # Average similarity to the seed set
        sim_2_seed = sim_slice.mean(0)
    elif agg=='min':
        # The minimum similarity to any gene in the seed set
        sim_2_seed = sim_slice.min(0)
    elif agg=='max':
        # The maximum similarity to any gene in the seed set
        sim_2_seed = sim_slice.max(0)
    elif agg=='perc':
        # The <agg_perc> percentile of similarities to the seed set.
        # For example, if a gene has similarities of (0, 0.2, 0.4,
        # 0.6, 0.8) to five seed genes, then the 10% similarity is 0.2
        sim_2_seed = np.percentile(100 * agg_perc, sim_slice, axis=0)
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

    try:
        if figure:
            import matplotlib.pyplot as plt
            import seaborn as sns
            fig = plt.figure()
            ax = plt.gca()
            if non_seed_idx.size > 1:
                sns.distplot(sim_2_seed[non_seed_idx], kde=False, norm_hist=True,
                             ax=ax,
                             axlabel='Similarity to seed set', label='Probability Density')
            if seed_idx.size > 1:
                sns.distplot(sim_2_seed[seed_idx], kde=False, norm_hist=True,
                             ax=ax,
                             axlabel='Similarity to seed set', label='Probability Density')
        else:
            fig = None
    except:
        fig = None

    return expand, expand_idx, sim_2_seed, fig

def ddot_pipeline(alpha,
                  beta,
                  gene_similarities,
                  genes,
                  seed,
                  ref,
                  name='Data-driven Ontology',
                  expand_kwargs={},
                  align_kwargs={},
                  ndex_args={'ndex_server' : config.ndex_server,
                             'ndex_user' : config.ndex_user,
                             'ndex_pass' : config.ndex_pass},
                  node_attr=None,
                  public=True,
                  verbose=False):
    """
    Assembles and analyzes a data-driven ontology to study a process or disease

    Parameters
    ----------

    alpha : float
       alpha parameter to CLIXO

    beta : float


    """

    ndex_server = ndex_args['ndex_server']
    ndex_user = ndex_args['ndex_user']
    ndex_pass = ndex_args['ndex_pass']

    ################
    # Expand genes #
    
    kwargs = {
        'agg':'mean',
        'min_sim':4,
        'filter_perc':None,
        'seed_perc':0,
        'agg_perc':None,
        'expand_size':200,
        'figure':True
    }
    kwargs.update(expand_kwargs)
        
    expand, expand_idx, sim_2_seed, fig = expand_seed(
        seed,
        gene_similarities,
        genes,
        **kwargs
    )
    
    expand = list(expand)
    if verbose:
        print 'Expanded gene set:', len(expand)
        # import matplotlib.pyplot as plt
        # from IPython.display import display
        # display(fig)
        # plt.close()
    
    try:
        similarity_uuid = gene_similarities
        gene_similarities, genes = ndex_to_sim_matrix(
            similarity_uuid,
            ndex_server,
            ndex_user,
            ndex_pass,
            similarity='similarity',
            input_fmt='cx_matrix',
            output_fmt='matrix',
            subset=None)
        similarity_from_ndex = True
    except:
        similarity_from_ndex = False
        assert isinstance(gene_similarities, np.ndarray)
    
    #############
    # Run CLIXO #
    
    df_sq = pd.DataFrame(gene_similarities[expand_idx, :][:, expand_idx], index=expand, columns=expand)
    df = melt_square(df_sq)

    from ddot.Ontology import Ontology    
    ont = Ontology.run_clixo(df, alpha, beta, verbose=verbose)

    try:
        ref = Ontology.from_ndex(go_uuid, ndex_server, ndex_user, ndex_pass)
    except:
        assert isinstance(ref, Ontology)
    
    #################
    # Align with GO #
    
    kwargs = {
        'iterations' : 3,
        'threads' : 4,        
    }
    kwargs.update(align_kwargs)
        
    alignment = ont.align(ref, **kwargs)
    if verbose:
        print 'Alignment: %s matches' % alignment.shape[0]
    
    ##################
    # Upload to NDEx #

    description = (
        'Data-driven ontology created by CLIXO '
        '(parameters: alpha={alpha}, beta={beta}). ').format(
            alpha=alpha,
            beta=beta)
    if similarity_from_ndex:
        description += (
            'Created from the similarity network '
            'at {ndex_server}/{ndex_uuid}').format(
                ndex_server=ndex_server,
                ndex_uuid=similarity_uuid)
    else:
        description += (
            'Created from a similarity network '
            'on a local file')
    
#     term_2_uuid = {t : [] for t in ont.terms}
    term_2_uuid = ont.upload_subnets_ndex(
        df,
        ['similarity'],
        name,
        ndex_server=ndex_server,
        ndex_user=ndex_user,
        ndex_pass=ndex_pass,
        propagate=True,
        public=public,
        verbose=False
    )
    
    # Annotate which genes were part of the seed set
    seed_set = set(seed)
    seed_attr = pd.DataFrame({'Seed' : [g in seed_set for g in ont.genes]}, index=ont.genes)
    if node_attr is None:
        node_attr = seed_attr
    else:
        node_attr = pd.concat([node_attr, seed_attr], axis=1)
    
    ont_ndex = ont.to_NdexGraph(
        name=name,
        description=description,
        term_2_uuid=term_2_uuid,
        layout='bubble',
        alignment=alignment,
        represents=True,
        node_attr=node_attr
    )
        
    ont_url = ont_ndex.upload_to(ndex_server, ndex_user, ndex_pass)
    ont_uuid = parse_ndex_uuid(ont_url)

    if public:
        make_network_public(ont_uuid)

    
    return ont, alignment, ont_ndex

def make_network_public(uuid,
                        ndex_server=config.ndex_server,
                        ndex_user=config.ndex_user,
                        ndex_pass=config.ndex_pass,
                        timeout=10,
                        error=False):
    import ndex.client as nc
    ndex = nc.Ndex(ndex_server, ndex_user, ndex_pass)

    start = time.time()
    while True:
        if time.time() - start > timeout:
            import traceback
            print traceback.print_exc()
            if error:
                raise Exception('Could not make the NDEX network %s public' % uuid)
            else:
                break
        else:
            try:
                ndex.make_network_public(uuid)
                break
            except:
                time.sleep(0.25)
