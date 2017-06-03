import numpy as np
import sys
from datetime import datetime
import pandas as pd
import networkx as nx

def time_print(*s):
#    import sys
#    from datetime import datetime
    print ' '.join(map(str, s)), datetime.today()
    sys.stdout.flush()

def load_table_as_list(f, delimiter='\t', encoding=None):    
    if encoding is None:
        f = open(f)
    else:
        import codecs
        f = codecs.open(f, encoding=encoding)
    tmp = [x.strip().split(delimiter) for x in f.read().splitlines()]
    f.close()
    return tmp    

def pivot_square(df, index, columns, values, fill_value=0):
    df = df.pivot(index=index, columns=columns, values=values)    
    index = df.index.union(df.columns)
    df = df.reindex(index=index, columns=index, fill_value=fill_value, copy=False)
    return df

def get_gene_name_converter(genes, scopes='symbol', fields='entrezgene', species='human'):
    import requests
    r = requests.post('http://mygene.info/v3/query',
                      data={'q': ','.join(genes),
                            'scopes': scopes,
                            'fields': fields,
                            'species': species})
    dic = {x['query'] : x[fields] for x in r.json() if x.has_key(fields)}
    return dic

def update_nx_with_alignment(G, alignment):
    for node_idx, node_attr in G.nodes(data=True):
        node_name = node_attr['name']
        if node_attr['Gene_or_Term']=='Term' and node_name in alignment.index:
            row = alignment.loc[node_name, :]
            node_attr['Label'] = row['Term_Description']
            node_attr['Aligned_Term'] = row['Term']
            node_attr['represents'] = row['Term']
            node_attr['Aligned_Term_Description'] = row['Term_Description']
            node_attr['Aligned_Similarity'] = row['Similarity']
            node_attr['Aligned_FDR'] = row['FDR']
        else:
            node_attr['Aligned_Term'] = 'NA'
            node_attr['Aligned_Term_Description'] = 'NA'
            node_attr['Aligned_Similarity'] = 0
            node_attr['Aligned_FDR'] = 1

def set_node_attributes_from_pandas(G, node_attr):
    G_nodes = set(G.nodes())
    node_attr = node_attr.loc[[x for x in node_attr.index if x in G_nodes], :]
    if node_attr is not None:
        for feat in node_attr.columns:
            nx.set_node_attributes(G, feat, node_attr[feat].to_dict())

def nx_nodes_to_pandas(G, attr_list=None):
    if attr_list is None:
        attr_list = list(set([a for d in G.nodes(data=True) for a in d[1].keys()]))
    return pd.concat([pd.Series(nx.get_node_attributes(G,a), name=a) for a in attr_list], axis=1)

def nx_edges_to_pandas(G, attr_list=None):
    if attr_list is None:
        attr_list = list(set([a for d in G.edges(data=True) for a in d[2].keys()]))
    return pd.concat([pd.Series(nx.get_edge_attributes(G,a), name=a) for a in attr_list], axis=1)

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
