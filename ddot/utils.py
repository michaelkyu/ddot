from __future__ import absolute_import, print_function

import sys
import base64
import time
import traceback
import os
import io
from math import ceil
from datetime import datetime

import pandas as pd
import networkx as nx
import numpy as np

import ndex.client as nc
from ndex.networkn import NdexGraph
    
import ddot
import ddot.config

def print_time(*s):
    print(' '.join(map(str, s)), datetime.today())
    sys.stdout.flush()

def to_hiview_url(ndex_url, hiview_server='http://hiview.ucsd.edu'):
    if hiview_server=='test':
        hiview_server = 'http://hiview-test.ucsd.edu'

    ndex_uuid = parse_ndex_uuid(ndex_url)
    ndex_server = parse_ndex_server(ndex_url)
    ndex_server = ndex_server.rstrip('/')
    if ndex_server=='http://test.ndexbio.org' or ndex_server=='http://dev2.ndexbio.org':
        server_type = 'test'
    elif ndex_server=='http://public.ndexbio.org' or ndex_server=='http://ndexbio.org':
        server_type = 'public'
    else:
        raise Exception('Invalid NDEx server URL: %s' % ndex_server)
    
    hiview_url = "%s/%s?type=%s&server=%s" % (hiview_server, ndex_uuid, server_type, ndex_server)
    return hiview_url

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
            if v in dic_inv:
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

def bubble_layout_nx(G, xmin=-750, xmax=750, ymin=-750, ymax=750, verbose=False):
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
    # Apply the 'Bubble Tree' graph layout plugin from Tulip
    graph.applyLayoutAlgorithm('Bubble Tree')
        
    viewLayout = graph.getLayoutProperty("viewLayout")
    pos = {g : (viewLayout[i].x(), viewLayout[i].y()) for i, g in zip(nodes, G.nodes())}
    pos = transform_pos(pos, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    return pos

def split_indices(n, k):
    try:
        tmp = iter(n)
        indices = n
    except TypeError:
        assert type(n)==type(int(1))
        indices = list(range(n))

    chunk_size = int(ceil(float(len(indices)) / k))

    return [(chunk_size * a, min(chunk_size * (a+1), len(indices))) for a in range(int(ceil(float(len(indices)) / chunk_size)))]

def split_indices_chunk(n, k):
    try:
        iter(n)
        n = len(n)
    except TypeError:
        assert isinstance(n, int)

    return [(k*i, min(k*(i+1), n)) for i in range(int(ceil(float(n) / k)))]

def make_index(it):
    """Create a dictionary mapping elements of an iterable to the index
    position of that element

    """
    return {b : a for a, b in enumerate(it)}

def time_print(*s):
    print(' '.join(map(str, s)), datetime.today())
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

# def get_gene_name_converter(genes, scopes='symbol', fields='entrezgene', species='human', target='gene'):
#     """Query mygene.info to get a dictionary mapping gene names in the ID
#     namespace scopes to the ID namespace in fields

#     Weird behavior with mygene.info: for Entrez genes, use fields
#     'entrezgene'. For ENSEMBL genes, use fields "ensembl"

#     """

#     if hasattr(genes, '__iter__') and not isinstance(genes, (str, unicode)):
#         genes = ','.join(genes)
        
#     import requests
#     r = requests.post('http://mygene.info/v3/query',
#                       data={'q': genes,
#                             'scopes': scopes,
#                             'fields': fields,
#                             'species': species})
    
    
#     def parse_field(x):
#         if isinstance(x, dict):
#             return [unicode(x[target])]
#         elif isinstance(x, list):
#             return [unicode(y[target]) for y in x]
#         else:
#             return unicode(x)
    
#     dic = {x['query'] : parse_field(x[fields]) for x in r.json() if x.has_key(fields)}
#     return dic

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
                G[e1][e2][feature_name] = v
#                G.edge[(e1,e2)][feature_name] = v

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
       DataFrame where index is the names of nodes and the columns are node attributes.

    """

    if attr_list is None:
        attr_list = list(set([a for d in G.nodes(data=True)
                              for a in d[1].keys()]))
    try:
        # Used for pandas version >= 0.23
        return pd.concat([pd.Series(nx.get_node_attributes(G,a), name=a) for a in attr_list],
                         axis=1, sort=True)
    except:
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

       DataFrame where index is a MultIndex with two levels (u,v)
       referring to edges and the columns refer to edge
       attributes. For multi(di)graphs, the MultiIndex have three
       levels of the form (u, v, key).

    """

    if attr_list is None:
        attr_list = list(set([a for d in G.edges(data=True) 
                              for a in d[2].keys()]))

    # print 'attr_list:', attr_list
    # for a in attr_list:
    #     print 'attr:', a
    #     nx.get_edge_attributes(G,a)
    
    if len(attr_list) > 0:
        df = pd.concat([pd.Series(nx.get_edge_attributes(G,a), name=a) for a in attr_list],
                         axis=1)
    else:
        df = pd.DataFrame(index=pd.MultiIndex.from_tuples(G.edges()))

    if df.index.nlevels==2:
        df.index.rename(['Node1', 'Node2'], inplace=True)
    elif df.index.nlevels==3:
        df.index.rename(['Node1', 'Node2', 'EdgeID'], inplace=True)
    else:
        raise Exception('Invalid number of levels: %s' % df.index.nlevels)

    return df    

def ig_nodes_to_pandas(G, attr_list=None):
    """Create pandas.DataFrame of node attributes of a igraph.Graph object.

    Parameters
    ----------
    G : igraph.Graph
    
    attr_list : list, optional
       Names of node attributes. Default: all node attributes

    Returns
    -------
    pandas.DataFrame
       DataFrame where index is the names of nodes and the columns are node attributes.

    """

    if attr_list is None:
        attr_list = G.vertex_attributes()

        # Remove the "name" attribute because it will be the DataFrame's index
        if 'name' in attr_list:
            attr_list.remove('name')

    df = pd.DataFrame(index=G.vs['name'])
    for attr in attr_list:
        df[attr] = G.vs[attr]
        
    df.dropna(axis=0, how='all', inplace=True)
    return df

def ig_edges_to_pandas(G, attr_list=None):
    """Create pandas.DataFrame of edge attributes of a igraph Graph object.

    Parameters
    ----------
    G : igraph.Graph
    
    attr_list : list, optional
       Names of edge attributes. Default: all edge attributes

    Returns
    -------
    pandas.DataFrame

       DataFrame where index is a MultIndex with two levels (u,v)
       referring to edges and the columns refer to edge
       attributes.

    """

    if attr_list is None:
        attr_list = G.edge_attributes()

    edge_list = [(G.vs[e.source]['name'], G.vs[e.target]['name']) for e in G.es]
    df = pd.DataFrame(index=pd.MultiIndex.from_tuples(edge_list))
    df.index.rename(['Node1', 'Node2'], inplace=True)
    for attr in attr_list:
        df[attr] = G.es[attr]

    df.dropna(axis=0, how='all', inplace=True)
    
    return df    

def nx_to_NdexGraph(G_nx, discard_null=True):
    """Converts a NetworkX into a NdexGraph object.

    Parameters
    ----------
    G_nx : networkx.Graph

    Returns
    -------
    ndex.networkn.NdexGraph

    """

    G = NdexGraph()
    node_id = 0
    node_dict = {}
    G.max_edge_id = 0
    for node_name, node_attr in G_nx.nodes(data=True):
        if discard_null:
            node_attr = {k:v for k,v in node_attr.items() if not pd.isnull(v)}

        if 'name' in node_attr:
            #G.add_node(node_id, node_attr)
            G.add_node(node_id, **node_attr)
        else:
            #G.add_node(node_id, node_attr, name=node_name)
            G.add_node(node_id, name=node_name, **node_attr)
        node_dict[node_name] = node_id
        node_id += 1
    for s, t, edge_attr in G_nx.edges(data=True):
        if discard_null:
            edge_attr = {k:v for k,v in edge_attr.items() if not pd.isnull(v)}

        G.add_edge(node_dict[s], node_dict[t], G.max_edge_id, edge_attr)
        G.max_edge_id += 1

    if hasattr(G_nx, 'pos'):
        G.pos = {node_dict[a] : b for a, b in G_nx.pos.items()}
        # G.subnetwork_id = 1
        # G.view_id = 1

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
    if 'v2/network/' in ndex_url:
        return ndex_url.split('v2/network/')[1]
    elif '#/network/' in ndex_url:
        return ndex_url.split('#/network/')[1]
    else:
        raise Exception("Not a valid NDEx URL: %s" % ndex_url)

def parse_ndex_server(ndex_url):
    if 'v2/network/' in ndex_url:
        return ndex_url.split('v2/network/')[0]
    elif '#/network/' in ndex_url:
        return ndex_url.split('#/network/')[0]
    else:
        raise Exception("Not a valid NDEx URL: %s" % ndex_url)

    # tmp = ndex_url.split('//')
    # if len(tmp) == 2:
    #     # e.g. 'http://dev2.ndexbio.org/v2/network/8bfa8318-55ed-11e7-a2e2-0660b7976219'
    #     return tmp[0] + '//' + tmp[1].split('v2/network/')[0]
    # elif len(tmp) == 1:
    #     # e.g. 'dev2.ndexbio.org/v2/network/8bfa8318-55ed-11e7-a2e2-0660b7976219'
    #     return tmp[0].split('v2/network/')[0]
    # elif len(tmp) == 0 or len(tmp) > 2:
    #     raise Exception()        

def create_edgeMatrix(X, X_cols, X_rows, verbose=True, G=None, ndex2=True):
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

    if ndex2:
        import ndex2
        import ndex2.client
        
        if not isinstance(X, np.ndarray):
            raise Exception('Provided matrix is not of type numpy.ndarray')
        if not isinstance(X_cols, list):
            raise Exception('Provided column header is not in the correct format.  Please provide a list of strings')
        if not isinstance(X_rows, list):
            raise Exception('Provided row header is not in the correct format.  Please provide a list of strings')

        if not X.flags['C_CONTIGUOUS']:
            X = np.ascontiguousarray(X)

        X_bytes = X.tobytes()
        chunk_size = int(1e8) # 100MB
        serialized_list = [{'v': base64.b64encode(X_bytes[s:e])} for i, (s, e) in enumerate(ddot.split_indices_chunk(len(X_bytes), chunk_size))]
        if verbose:
            print('Broke up serialization into %s chunks' % len(serialized_list))

        if verbose:
            print('Size of numpy array (MB):', X.nbytes / 1e6)
            serialize_size = sum([sys.getsizeof(x['v']) for x in serialized_list])
            print('Size of serialization (MB):', serialize_size / 1e6)
            print('Constant factor overhead:', float(serialize_size) / X.nbytes)
            
#         serialized = base64.b64encode(X.tobytes())

#         if verbose:
#             print('Size of numpy array (MB):', X.nbytes / 1e6)
#             print('Size of serialization (MB):', sys.getsizeof(serialized) / 1e6)
#             print('Constant factor overhead:', float(sys.getsizeof(serialized)) / X.nbytes)

#         chunk_size = int(1e8)
# #        chunk_size = int(1e1)
#         serialized_list = [{'v': serialized[s:e]} for i, (s, e) in enumerate(ddot.split_indices_chunk(len(serialized), chunk_size))]
#         if verbose:
#             print('Broke up serialization into %s chunks' % len(serialized_list))


        nice_cx_builder = ndex2.NiceCXBuilder()
        # nice_cx_builder.set_name(name)
        nice_cx_builder.add_node(name='Matrix', represents='Matrix')

        # nice_cx_builder.add_opaque_aspect('matrix', [{'v': serialized}])
        #nice_cx_builder.add_opaque_aspect('matrix', [serialized_list])
        nice_cx_builder.add_opaque_aspect('matrix', serialized_list)
        nice_cx_builder.add_opaque_aspect('matrix_cols', [{'v': X_cols}])
        nice_cx_builder.add_opaque_aspect('matrix_rows', [{'v': X_rows}])
        nice_cx_builder.add_opaque_aspect('matrix_dtype', [{'v': X.dtype.name}])

        nice_cx = nice_cx_builder.get_nice_cx()

        return nice_cx
    else:
        if not X.flags['C_CONTIGUOUS']:
            X = np.ascontiguousarray(X)

        # Use base64 encoding of binary to text. More efficient than
        # pickle(*, protocol=0)
        start = time.time()
        serialized = base64.b64encode(X)
        base64_time = time.time() - start

        assert isinstance(X_cols, list)
        assert isinstance(X_rows, list)

        if sys.version_info.major==3:
            start = time.time()
            serialized = serialized.decode('utf-8')
            serialize_decode_time = time.time() - start
            if verbose:
                print('serialize_decode_time (sec):', serialize_decode_time)

        if verbose:
            print('base64 encoding time (sec):', base64_time)
            print('Size of numpy array (MB):', X.nbytes / 1e6)
            print('Size of serialization (MB):', sys.getsizeof(serialized) / 1e6)
            print('Constant factor overhead:', float(sys.getsizeof(serialized)) / X.nbytes)

        if G is None:
            G = NdexGraph()
        # G.unclassified_cx.append(
        #     {'matrix': serialized,
        #      'matrix_cols' : X_cols,
        #      'matrix_rows' : X_rows,
        #      'matrix_dtype' : X.dtype.name})

        G.unclassified_cx.append({'matrix': [{'v': serialized}]})
        G.unclassified_cx.append({'matrix_cols': [{'v': X_cols}]})
        G.unclassified_cx.append({'matrix_rows': [{'v': X_rows}]})
        G.unclassified_cx.append({'matrix_dtype': [{'v': X.dtype.name}]})

        G.add_new_node('Matrix')
        
        return G

def load_edgeMatrix(ndex_uuid,
                    ndex_server,
                    ndex_user,
                    ndex_pass,
                    ndex=None,
                    json=None,
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

        JSON module with "loads" function. Default: the simplejson
        package (must be installed)

    Returns
    -------
    X : np.ndarray
    
    X_cols : list
        Column names

    X_rows : list
        Row names

    """

    if json is None:
        # import ijson
        # json = ijson
        
        import simplejson
        json = simplejson
        
    if ndex is None:
        ndex = nc.Ndex(ndex_server, ndex_user, ndex_pass)

    start = time.time()
    response = ndex.get_network_as_cx_stream(ndex_uuid)
    response.raw.decode_content = True
    try:
        cx = json.load(response.raw)
    finally:
        response.close()
    if verbose:
        print('NDEx download and CX parse time (sec):', time.time() - start)

    start_loop = time.time()

    for aspect in cx:
        if 'matrix_cols' in aspect:
            cols = aspect.get('matrix_cols')[0].get('v')
        if 'matrix_rows' in aspect:
            rows = aspect.get('matrix_rows')[0].get('v')
        if 'matrix_dtype' in aspect:
            dtype = np.dtype(aspect.get('matrix_dtype')[0].get('v'))

    dim = (len(rows), len(cols))
    #X_buf = bytearray(dim[0] * dim[1] * np.dtype(dtype).itemsize)
    X_buf = []
    pointer = 0

    if verbose:
        print('Dim:', dim)
        # print('Bytes:', len(X_buf))
        print('Bytes:', dim[0] * dim[1] * np.dtype(dtype).itemsize)
    
    for aspect in cx:
        if 'matrix' in aspect:
            for x in aspect.get('matrix'):
                if sys.version_info.major==3:
                    binary_data = base64.decodebytes(x.get('v').encode('utf-8'))
                else:
                    binary_data = base64.b64decode(x.get('v'))
                del x['v']
                X_buf.append(binary_data)
                
                # X_buf[pointer : pointer + len(binary_data)] = binary_data
                # pointer += len(binary_data)

    X_buf = (b"").join(X_buf)
        
    # Create a NumPy array
    X = np.frombuffer(X_buf, dtype=dtype).reshape(dim)

    if verbose:
        print('Iterate through CX and construct array time (sec):', time.time() - start_loop)
    
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
        return create_edgeMatrix(sim, names, names, ndex2=False)

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

def ndex_to_sim_matrix(ndex_url,
                       ndex_server=None,
                       ndex_user=None,
                       ndex_pass=None,
                       similarity=None,
                       input_fmt='cx_matrix',
                       output_fmt='matrix',
                       subset=None,
                       verbose=True):
    """Read a similarity network from NDEx and return it as either a
    square np.array (compact representation) or a pandas.DataFrame of
    the non-zero similarity values (sparse representation)

    Parameters
    ----------
    ndex_url : str
        NDEx URL (or UUID) of ontology
    
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

    if ndex_server is None:
        ndex_server = ddot.config.ndex_server
    if ndex_user is None:
        ndex_pass = ddot.config.ndex_user
    if ndex_pass is None:
        ndex_pass = ddot.config.ndex_pass

    if 'http' in ndex_url:
        ndex_server = parse_ndex_server(ndex_url)
        ndex_uuid = parse_ndex_uuid(ndex_url)
    else:
        ndex_uuid = ndex_url
    
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
            ndex_pass,
            verbose=verbose)
        assert sim_names == sim_names_col

        if subset is not None:
            idx = make_index(sim_names)
            idx = np.array([idx[g] for g in subset if g in idx])
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
                figure=False,
                verbose=False):
    """Identify genes that are most similar to a seed set of genes.

    A gene is included in the expanded set only if it meets all of the
    specified criteria. If include_seed is True, then genes that are
    seeds will be included regardless of the criteria. At the same time,
    the number of genes returned is still limited by expand_size. One way
    to get n novel genes returned is therefore to set expand_size = n + |seed| and
    include_seed = True, and then to remove the seed list from expand.

    Parameters
    ----------
    seed : list

    sim : np.ndarray

    sim_names : list of str

    agg : str or function
    
       Aggregation method. Possible values are mean, min, max, perc.

    min_sim : float
    
       Minimum similarity to the seed set.

    filter_perc : float
    
       Filter based on a percentile of similarities between all genes and the seed set.
    
    seed_perc : float
    
       Filter based on a percentile of similarities between seed set to itself.

    agg_perc : float
    
       The <agg_perc> percentile of similarities to the seed set.
       For example, if a gene has similarities of (0, 0.2, 0.4,
       0.6, 0.8) to five seed genes, then the 10% similarity is 0.2

    expand_size : int
    
       Maximum limit on the number of returned genes.

    include_seed : bool
	Include the seed genes even if they didn't meet the criteria.

    figure : bool
    
       Generate a figure showing the average distances within the seed an d the average distances between seed and the background.

    Returns
    -------
    expand
    
       The list of expanded genes passing all filters.

    expand_idx
    
       Indices of the ranking. I.e. `expand_idx[0]` is the index of the top
       gene, so you can get the name of the top gene
       with `sim_names[expand_idx[0]]` where `sim_names` is the input parameter.

    sim_2_seed
    
       The returned array `sim_2_seed` is the calculated similarities of the
       genes to the seed set. So `sim_2_seed[0]` is the similarity of the gene 
    
    fig
    
       The generated figure. Can be saved like this: plt.savefig('foo.pdf')

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
        # Don't include self in the average.
        sim_2_seed = sim_slice.sum(0)
        sim_2_seed[seed_idx] -= sim[seed_idx, seed_idx]
        sim_2_seed[seed_idx] /= float(seed_idx.size - 1)
        sim_2_seed[np.setdiff1d(np.arange(sim_2_seed.size), seed_idx)] /= float(seed_idx.size)

        ## Just take a simple average including self.
        ## This method has the problem is that it's not clear what to set the similarity to self
        # sim_2_seed = sim_slice.mean(0)
        
        # print sim_2_seed.size, np.diagonal(sim_slice).size, sim_slice.shape        
        # sim_2_seed -= np.diagonal(sim_slice)
        # sim_2_seed = sim_2_seed / float(sim_2_seed.shape[0] - 1)
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
            
        if verbose: print('min_sim:', min_sim)

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

        try:
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
            try:
                figure.savefig(figure)
            except:
                pass
        except:
            fig = None
    else:
        fig = None

    return expand, expand_idx, sim_2_seed, fig

def make_seed_ontology(sim,
                       sim_names,
                       expand_kwargs={},
                       build_kwargs={},
                       align_kwargs={},
                       ndex_kwargs={},
                       node_attr=None,
                       verbose=False,
                       ndex=True):
    """Assembles and analyzes a data-driven ontology to study a process or disease

    Parameters
    ----------

    sim : np.ndarray

       gene-by-gene similarity array

    sim_names : array-like

       Names of genes as they appear in the rows and columns of <sim>

    expand_kwargs : dict

       Parameters for ddot.expand_seed() to identify an expanded set of genes

    build_kwargs : dict

       Parameters for Ontology.build_from_network(...) to build a data-driven ontology.

    align_kwargs : dict

       Parameters for Ontology.align() to align against a reference ontology.
    
    ndex_kwargs : dict

       Parameters for Ontology.to_ndex() to upload ontology to NDEx.

    node_attr : pd.DataFrame

       A DataFrame of node attributes to assign to the ontology.

    ndex : bool

       If True, then upload ontology to NDEx using parameters <ndex_kwargs>
    """

    assert 'seed' in expand_kwargs
    seed = expand_kwargs['seed']
    
    ################
    # Expand genes #
    ################
    if verbose:
        print('----------------')
        print('Expanding genes')
        print('----------------')
    
    kwargs = {'sim': sim,
              'sim_names': sim_names}
    kwargs.update(expand_kwargs)
        
    expand, expand_idx, sim_2_seed, fig = expand_seed(**kwargs)
    expand_results = {'expand' : expand, 'sim_2_seed' : sim_2_seed, 'fig' : fig}    
    expand = list(expand)

    if verbose:
        print('Seed genes:', len(seed))
        print('Expand genes:', len(expand))
        
    ##################
    # Build Ontology #
    ##################
    if verbose:
        print('-----------------')
        print('Building ontology')
        print('-----------------')
    
    # Slice the similarity matrix over the expanded gene set and
    # convert to a square dataframe
    df_sq = pd.DataFrame(sim[expand_idx, :][:, expand_idx], index=expand, columns=expand)
    # Convert square to long format
    df = melt_square(df_sq)

    # Build data-driven ontology
    ont = ddot.Ontology.infer_ontology(df, verbose=verbose, **build_kwargs)

    ###############################
    # Align to Reference Ontology #
    ###############################

    if 'hier' in align_kwargs:
        if verbose:
            print('------------------')
            print('Aligning Ontology')
            print('------------------')
                
        alignment = ont.align(**align_kwargs)
        if verbose:
            print('Alignment: %s alignment matches' % alignment.shape[0])

    #############################
    # Set other node attributes #
    #############################
    
    # Annotate which genes were part of the seed set
    seed_set = set(seed)
    seed_attr = pd.DataFrame({'Seed' : [g in seed_set for g in ont.genes]}, index=ont.genes)
    ont.update_node_attr(seed_attr)

    # Annotate the data similarity to the seed set
    tmp = make_index(sim_names)
    sim_attr = pd.DataFrame({'Similarity_2_Seed' : [sim_2_seed[tmp[g]] for g in ont.genes]}, index=ont.genes)
    ont.update_node_attr(sim_attr)    

    # Annotate user-specified node attributes
    if node_attr is not None:
        ont.update_node_attr(node_attr)

    # Color seed genes as green (hex #6ACC65)
    fill_attr = pd.DataFrame({'Vis:Fill Color' : '#6ACC65'}, index=seed)
    ont.update_node_attr(fill_attr)        

    # Color terms according to the exactness of alignment
    if 'Aligned_Similarity' in ont.node_attr.columns:
        fill_attr = ont.node_attr['Aligned_Similarity'].dropna().map(color_gradient)
        fill_attr = fill_attr.to_frame().rename(columns={'Aligned_Similarity' : 'Vis:Fill Color'})
        ont.update_node_attr(fill_attr)

    ##################
    # Upload to NDEx #
    ##################
        
    if ndex:
        if verbose:
            print('--------------------------')
            print('Uploading Ontology to NDEx')
            print('--------------------------')           
        
        if 'network' not in ndex_kwargs:
            ndex_kwargs['network'] = df
            ndex_kwargs['features'] = ['similarity']
            ndex_kwargs['main_feature'] = 'similarity'

        description = (
            'Data-driven ontology created by the function ddot.make_seed_ontology()'
            'in the DDOT Python package (https://github.com/michaelkyu/ontology)'
            '(parameters: %s' % ', '.join(['%s=%s' % (k,v) for k,v in build_kwargs.items()])
        )

        ont_url, ont_ndexgraph = ont.to_ndex(
            description=description,
            verbose=verbose,
            **ndex_kwargs
        )
    else:
        ont_url, ont_ndexgraph = None, None

    return ont, ont_url, ont_ndexgraph, expand_results

def make_network_public(uuid,
                        ndex_server,
                        ndex_user,
                        ndex_pass,
                        timeout=60,
                        error=False):
    ndex = nc.Ndex(ndex_server, ndex_user, ndex_pass)
            
    sleep_time = 0.25
    
    start = time.time()
    while True:
        if time.time() - start > timeout:
            print('Failed to make network public: error message:')
            print(traceback.print_exc())
            if error:
                raise Exception('Could not make the NDEX network %s public' % uuid)
            else:
                break
        else:
            try:
                ndex.make_network_public(uuid)
                break
            except:
                time.sleep(sleep_time)
                sleep_time = min(5, 2 * sleep_time)


def ig_unfold_tree_with_attr(g, sources, mode):
    """Call igraph.Graph.unfold_tree while preserving vertex and edge
    attributes.

    """
    
    g_unfold, g_map = g.unfold_tree(sources, mode=mode)
    
    g_eids = g.get_eids([(g_map[e.source], g_map[e.target]) for e in g_unfold.es])
    for attr in g.edge_attributes():
        g_unfold.es[attr] = g.es[g_eids][attr]
        
    for attr in g.vertex_attributes():
        g_unfold.vs[attr] = g.vs[g_map][attr]
                                        
    return g_unfold


def gridify(parents, pos, G):
    """Relayout leaf nodes into a grid.

    Nodes must be connected and already laid out in "star"-like
    topologies. In each "star", a set of nodes are positioned to form
    the shape of a circle and connect to a common parent node that is
    positioned at the circle's center.

    This function repositions the nodes in each start into a square
    grid that inscribes the circle.

    Parameters
    ----------
    parents : list
        
        For each parent, its children will be arranged in a grid.

    pos : dict

        Dictionary that maps names of nodes to their (x,y) coordinates
    
    G : nx.Graph

        Network

    Returns
    -------
    : None

        Modifies <pos> inplace

    """
    
    for v in parents:
        x_center, y_center = pos[v]
        children = list(G.predecessors(v))
        if len(children) <= 1:
            continue
    
        # Estimate radius by averaging the distance to the children
        radius = np.mean([np.sqrt((pos[c][0] - x_center)**2 + (pos[c][1] - y_center)**2) for c in children])
        
        # Width of the the square inscribing the circle is sqrt(2) *
        # radius
        width = np.sqrt(2) * radius

        # Reposition nodes into a square grid
        num_cols = int(np.ceil(np.sqrt(len(children))))
        col_space = width / (num_cols - 1)
        row_space = width / (num_cols - 1)
        topleft = x_center - width / 2, y_center - width / 2
        for i, c in enumerate(children):
            row, col = i / num_cols, i % num_cols
            pos[c] = (topleft[0] + col * col_space, topleft[1] + row * row_space)

        # Reposition the center node to be 1/3rd of the distance from
        # the outside of the circle to the center's parent
        p = list(G.successors(v))[0]
        x_parent, y_parent = pos[p]
        distance = np.sqrt((x_parent - x_center)**2 + (y_parent - y_center)**2)
        alpha = (2/3.) * (distance - radius) / distance
        
        pos[v] = (x_center * (1 - alpha) + x_parent * alpha,
                  y_center * (1 - alpha) + y_parent * alpha)
        
def nx_set_tree_edges(G, tree_edges):
    nx.set_edge_attributes(
        G,
        values={(s,t) : 'Tree' if ((s,t) in tree_edges) else 'Not_Tree'
         for s, t in G.edges(data=False)},
        name='Is_Tree_Edge'
    )

def color_gradient(ratio, min_col='#FFFFFF', max_col='#D65F5F', output_hex=True):
    """Calculate a proportional mix between two colors.

    """

    min_col_hex = min_col.lstrip('#')
    min_col_rgb = tuple(int(min_col_hex[i:i+2], 16) for i in (0, 2 ,4))
    max_col_hex = max_col.lstrip('#')
    max_col_rgb = tuple(int(max_col_hex[i:i+2], 16) for i in (0, 2 ,4))

    mix_col_rgb = [int(ratio*x + (1-ratio)*y) for x, y in zip(max_col_rgb, min_col_rgb)]
    if output_hex:
        return('#%02x%02x%02x' % tuple(mix_col_rgb)).upper()
    else:
        return mix_col_rgb

