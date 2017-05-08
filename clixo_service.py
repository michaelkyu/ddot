import ndex.client as nc
import io
import json
from IPython.display import HTML
from time import sleep
import os, time, tempfile
import sys
import time
import logging

import grpc
import networkx as nx
import cx_pb2
import cx_pb2_grpc

import numpy as np
import inspect
from concurrent import futures
from itertools import combinations
from subprocess import Popen, PIPE, STDOUT
import pandas as pd

from ndex.networkn import NdexGraph

from clixo import run_clixo
from Ontology import Ontology, networkx_to_NdexGraph, make_tree
from utilities import pivot_2_square

_ONE_DAY_IN_SECONDS = 60 * 60 * 24

class CyServiceServicer(cx_pb2_grpc.CyServiceServicer):

    def StreamElements(self, element_iterator, context):
        try:
            clixo_params = {'ndex_uuid' : None,
                            'ndex_server' : 'http://public.ndexbio.org',
                            'similarity_attr' : 'similarity',
                            'name' : 'Data-Driven Ontology',
                            'dt_thresh': -100000,
                            'max_time': 100000,
                            'clixo_folder': os.getenv('CLIXO'),
                            'output_fmt': 'cx',
                            'verbose': True}

            # print 'context:'
            # print context
            # print dir(context)
            # 0 / asdf

            input_G, clixo_params, errors = self.read_element_stream(element_iterator, clixo_params)


            print 'Parameters:'
            print clixo_params

            ###############################

            if isinstance(clixo_params['ndex_uuid'], (str, unicode)):
                # Read graph using NDEx client
                input_G = NdexGraph(server=clixo_params['ndex_server'],
                                    username=clixo_params['ndex_user'],
                                    password=clixo_params['ndex_pass'],
                                    uuid=clixo_params['ndex_uuid'])
                similarity_attr = clixo_params['similarity_attr']
                graph = [(input_G.node[u]['name'],
                          input_G.node[v]['name'],
                          float(attr[similarity_attr])) for u, v, attr in input_G.edges_iter(data=True)]
                clixo_params['graph'] = graph
            else:
                # Read graph from CXmate stream
                #graph = [(nodes_dict[v1], nodes_dict[v2], s) for v1, v2, s in edges_dict.values()]
                graph = [(u, v, float(attr['similarity'])) for u, v, attr in input_G.edges_iter(data=True)]
                clixo_params['graph'] = graph

#            if len(errors) == 0:

            if True:
                ## Run CLIXO
#                with tempfile.NamedTemporaryFile('w', delete=False) as output:
                with open('/tmp/tmp.txt', 'w') as output:
                    clixo_params['output'] = output.name

                    clixo_argnames = inspect.getargspec(run_clixo).args

                    print 'Input to CLIXO:'
                    print {k : v for k, v in clixo_params.items() if k in clixo_argnames}.keys()

                    run_clixo(**{k : v for k, v in clixo_params.items() if k in clixo_argnames})

                    with open(clixo_params['output'], 'r') as f:
                        lines = [x.split('\t') for x in f.read().splitlines() if len(x)>0 and x[0]!='#']

                ## Read using Ontology class (to make tree, and to get term sizes)
                ontology_table = [x for x in lines if x[2]!='gene']
                mapping_table = [x for x in lines if x[2]=='gene']
                ontology = Ontology(ontology_table, mapping_table, parent_child=True)

                ontology_propagated = Ontology(ontology_table, mapping_table, parent_child=True)
                ontology_propagated.propagate_annotations()
                term_sizes = dict(zip(ontology_propagated.terms, ontology_propagated.get_term_sizes()))

                # Convert gene-gene similarity network into a square matrix
                graph_pd = pd.DataFrame(graph, columns=['Gene1', 'Gene2', 'Sim'])
                graph_pd['Sim'] = graph_pd['Sim'].astype(np.float32)
                mat = graph_pd.pivot(index='Gene1', columns='Gene2', values='Sim')
                arr, arr_genes_index = pivot_2_square(mat)

                # Create an NDEX network for every term's subnetwork
                term_2_url = self.upload_subnetworks_2_ndex(ontology_propagated, arr, arr_genes_index,
                                                            clixo_params['ndex_server'], clixo_params['ndex_user'], clixo_params['ndex_pass'],
                                                            clixo_params['name'])
                term_2_uuid = {t : url.split('http://public.ndexbio.org/v2/network/')[1] for t, url in term_2_url.items()}

                ont_tmp = Ontology(ontology_table, mapping_table, parent_child=True, genes_as_terms=True)
                ont_tmp.propagate_annotations()
                tree = make_tree(ont_tmp.get_igraph(),
                                 parent_priority=ont_tmp.get_term_sizes(), optim=min)
                tree_edges = set([(tree.vs[e.source]['name'], tree.vs[e.target]['name']) for e in tree.es if e['smallest_parent']])

                if clixo_params['output_fmt']=='cx':
                    # Use CXmate to stream to NDEX
                    for elt in self.stream_ontology(ontology, term_sizes, term_2_uuid, tree_edges):
                        yield elt

                elif clixo_params['output_fmt']=='ndex':
                    # Directly output to NDEX
                    ontology_ndex = ontology.to_networkx()

                    # import cPickle
                    # with open('/cellar/users/mikeyu/DeepTranslate/tmp.pickle', 'wb') as f:
                    #     cPickle.dump({'ontology_ndex':ontology_ndex, 'ontology':ontology}, f, protocol=2)

                    for t in ontology.terms:
                        ontology_ndex.node[t]['name'] = 'CLIXO:%s' % t
                        ontology_ndex.node[t]['Gene_or_Term'] = 'Term'
                        ontology_ndex.node[t]['Size'] = term_sizes[t]
                        ontology_ndex.node[t]['ndex:internalLink'] = term_2_uuid[t]
                    for g in ontology.genes:
                        ontology_ndex.node[g]['name'] = g
                        ontology_ndex.node[g]['Size'] = 1
                        ontology_ndex.node[g]['Gene_or_Term'] = 'Gene'
                    nx.set_edge_attributes(ontology_ndex,
                                           'Is_Tree_Edge',
                                           {(s,t) : 'Tree' if ((s,t) in tree_edges) else 'Not_Tree' \
                                            for s, t in ontology_ndex.edges_iter(data=False)})

                    ontology_ndex = networkx_to_NdexGraph(ontology_ndex)
                    ontology_ndex.set_name(clixo_params['name'])

                    ontology_url = ontology_ndex.upload_to(clixo_params['ndex_server'], clixo_params['ndex_user'], clixo_params['ndex_pass'])
                    ontology_uuid = ontology_url.split('http://public.ndexbio.org/v2/network/')[1]
                    print 'ontology_url:', ontology_url 

                    # description = 'Data-driven ontology created by CLIXO (parameters: alpha=%s, beta=%s).' % (clixo_params['alpha'], clixo_params['beta'])
                    # if isinstance(clixo_params['ndex_uuid'], (str, unicode)):
                    #     description += ' Created from similarity network at %s/%s' % (clixo_params['ndex_server'], clixo_params['ndex_uuid'])
                    # my_ndex=nc.Ndex(clixo_params['ndex_server'], clixo_params['ndex_user'], clixo_params['ndex_pass'])
                    # ontology_profile = {'description': description}
                    # my_ndex.update_network_profile(ontology_uuid, ontology_profile)

                    element = cx_pb2.Element()
                    netAttr = element.networkAttribute
                    netAttr.name = 'ndex_uuid'
                    netAttr.value = ontology_uuid
                    yield element

                    # param = element.parameter
                    # param.name = 'ndex_uuid'
                    # param.value = ontology_uuid
                    # yield element

            else:
                for caught_error in errors:
                    error = self.create_internal_crash_error(caught_error.message, 500)
                    log_error(error)
                    yield error
        except Exception as e:
            message = "Unexpected error: " + str(e)
            error = self.create_internal_crash_error(message, 500)
            log_error(error)
        
            import traceback
            print traceback.print_exc()

            yield error

    def stream_ontology(self, ontology, term_sizes, term_2_uuid, tree_edges):
        node_2_id = {}
        node_id = 0
        for node_name in ontology.genes:
            yield self.create_node(node_id, node_name)
            yield self.create_nodeAttribute(node_id, 'Gene_or_Term', 'Gene')
            yield self.create_nodeAttribute(node_id, 'Size', '1')
            node_2_id[node_name] = node_id
            node_id += 1
        for node_name in ontology.terms:
            yield self.create_node(node_id, node_name)
            yield self.create_nodeAttribute(node_id, 'Gene_or_Term', 'Term')
            yield self.create_nodeAttribute(node_id, 'ndex:internalLink', term_2_uuid[node_name])
            yield self.create_nodeAttribute(node_id, 'Size', str(term_sizes[node_name]))
            node_2_id[node_name] = node_id
            node_id += 1

        edge_id = 0
        for g in ontology.genes:
            for t_i in ontology.gene_2_terms[g]:
                t = ontology.terms[t_i]
                yield self.create_edge(edge_id, node_2_id[g], node_2_id[t])
                yield self.create_edgeAttribute(edge_id, 'EdgeType', 'Gene-Term Annotation')
                if (g, t) in tree_edges:
                    yield self.create_edgeAttribute(edge_id, 'Is_Tree_Edge', 'Tree')
                else:
                    yield self.create_edgeAttribute(edge_id, 'Is_Tree_Edge', 'Not_Tree')
                edge_id += 1
        for p, c_list in ontology.term_2_terms.iteritems():        
            for c in c_list:
                yield self.create_edge(edge_id, node_2_id[c], node_2_id[p])
                yield self.create_edgeAttribute(edge_id, 'EdgeType', 'Child-Parent Hierarchical Relation')
                if (c,p) in tree_edges:
                    yield self.create_edgeAttribute(edge_id, 'Is_Tree_Edge', 'Tree')
                else:
                    yield self.create_edgeAttribute(edge_id, 'Is_Tree_Edge', 'Not_Tree')
                edge_id += 1

    def upload_subnetworks_2_ndex(self, ontology, arr, arr_genes_index,
                                  ndex_server, ndex_user, ndex_pass, name):
        """Push subnetworks"""

        #print ontology.get_term_2_genes()

        term_2_url = {}
        for t in ontology.terms:
            #print 't:', t
            genes = np.array([ontology.genes[g] for g in ontology.get_term_2_genes()[t]])
            #print 'genes:', genes
            idx = np.array([arr_genes_index[g] for g in genes])
            #print 'idx:', idx
            subarr = arr[idx,:][:,idx]

            # Set nan to 0
            subarr[np.isnan(subarr)] = 0

            row, col = subarr.nonzero()
            row, col = row[row < col], col[row < col]
                        
            G = NdexGraph()    
            G.create_from_edge_list(zip(genes[row], genes[col]))
            for i in np.arange(row.size):
                G.set_edge_attribute(i+1, "similarity", str(subarr[row[i], col[i]]))
            G.set_name('%s supporting network for CLIXO:%s' % (name, t))
            G.set_network_attribute('Description', '%s supporting network for CLIXO:%s' % (name, t))

            ndex_url = G.upload_to(ndex_server, ndex_user, ndex_pass)
            term_2_url[t] = ndex_url

        return term_2_url

    def create_node(self, node_id, node_name):
        element = cx_pb2.Element()
        node = element.node
        node.id = node_id
        node.name = node_name
        return element

    def create_nodeAttribute(self, node_id, key, val):
        element = cx_pb2.Element()
        attr = element.nodeAttribute
        attr.nodeId = node_id
        attr.name = key
        attr.value = val
        return element

    def create_edge(self, edge_id, node1, node2):
        element = cx_pb2.Element()
        edge = element.edge
        edge.id = edge_id
        edge.sourceId = node1
        edge.targetId = node2
        return element

    def create_edgeAttribute(self, edge_id, key, val):
        element = cx_pb2.Element()
        attr = element.edgeAttribute
        attr.edgeId = edge_id
        attr.name = key
        attr.value = val
        return element

    # def create_output_attribute(self, node_id, value, attribute_name, suffix):
    #     element = cx_pb2.Element()
    #     attr = element.nodeAttribute
    #     attr.nodeId = node_id
    #     attr.name = attribute_name + suffix
    #     attr.value = value
    #     return element

    def create_internal_crash_error(self, message, status):
        element = cx_pb2.Element()
        error = element.error
        error.status = status
        error.code = 'clixo_service/' + str(status)
        error.message = message
        error.link = 'https://github.com/michaelkyu/ontology_cx'
        return element

    def read_element_stream(self, element_iter, parameters):
        errors = []
        edgesAttr_dict = {}
        nodesAttr_dict = {}
        edges_dict = {}
        nodes_dict = {}

        for element in element_iter:
            ele_type = element.WhichOneof('value')
            if ele_type == 'error':
                errors.append(element.error)
            elif ele_type == 'parameter':
                param = element.parameter
                parameters[param.name] = param.value
            elif ele_type == 'node':
                node = element.node
                nodes_dict[node.id] = node.name
            elif ele_type == 'edge':
                edge = element.edge
                edges_dict[edge.id] = (edge.sourceId, edge.targetId)
            elif ele_type == 'nodeAttribute':
                pass
            elif ele_type == 'edgeAttribute':
                edgeAttr = element.edgeAttribute
                if edgesAttr_dict.has_key(edgeAttr.name):
                    edgesAttr_dict[edgeAttr.name][edgeAttr.edgeId] = edgeAttr.value
                else:
                    edgesAttr_dict[edgeAttr.name] = {edgeAttr.edgeId : edgeAttr.value}

#        print 'nodes_dict:', nodes_dict


        G = nx.Graph()
        for n_id, u in nodes_dict.iteritems():
            G.add_node(u, node_id=n_id)
        edge_attributes_list = edgesAttr_dict.keys()
        for e_id, (u, v) in edges_dict.iteritems():
            G.add_edge(nodes_dict[u], nodes_dict[v],
                       attr_dict={k : edgesAttr_dict[k][e_id] for k in edge_attributes_list if edgesAttr_dict[k].has_key(e_id)},
                       edge_id=e_id)

#        print 'G nodes:', G.nodes()
#        0 / asdf

        return G, parameters, errors

    # def read_element_stream(self, element_iter, parameters):
    #     errors = []
    #     edges_dict = {}
    #     nodes_dict = {}

    #     for element in element_iter:
    #         ele_type = element.WhichOneof('value')
    #         if ele_type == 'error':
    #             errors.append(element.error)
    #         elif ele_type == 'parameter':
    #             param = element.parameter
    #             parameters[param.name] = param.value
    #         elif ele_type == 'node':
    #             node = element.node
    #             nodes_dict[node.id] = node.name
    #         elif ele_type == 'edge':
    #             edge = element.edge
    #             if edges_dict.has_key(edge.id):
    #                 edges_dict[edge.id][:2] = [edge.sourceId, edge.targetId]
    #             else:
    #                 edges_dict[edge.id] = [edge.sourceId, edge.targetId, None]
    #         elif ele_type == 'nodeAttribute':
    #             pass
    #         elif ele_type == 'edgeAttribute':
    #             edgeAttr = element.edgeAttribute
    #             if edgeAttr.name == 'similarity':
    #                 if edges_dict.has_key(edgeAttr.edgeId):
    #                     edges_dict[edgeAttr.edgeId][2] = float(edgeAttr.value)
    #                 else:
    #                     edges_dict[edgeAttr.edgeId] = [None, None, float(edgeAttr.value)]
            
    #     return (nodes_dict, edges_dict), parameters, errors

def log_info(message):
    logging.info(message)

def log_error(message):
    logging.error(message)

def serve():
    server = grpc.server(futures.ThreadPoolExecutor(max_workers=10))
    cx_pb2_grpc.add_CyServiceServicer_to_server(
            CyServiceServicer(), server)
    server.add_insecure_port('0.0.0.0:8080')
    server.start()
    try:
        while True:
            time.sleep(_ONE_DAY_IN_SECONDS)
    except KeyboardInterrupt:
        server.stop(0)

if __name__ == '__main__':
    logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
    log_info("Listening for requests on '0.0.0.0:8080'")
    serve()
