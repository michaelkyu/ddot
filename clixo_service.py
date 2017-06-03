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

from Ontology import Ontology
from utils import pivot_square, nx_to_NdexGraph, NdexGraph_to_nx, nx_nodes_to_pandas, nx_edges_to_pandas
from cx_utils import yield_ndex, parse_ndex_uuid
import config

_ONE_DAY_IN_SECONDS = 60 * 60 * 24

class CyServiceServicer(cx_pb2_grpc.CyServiceServicer):

    def StreamElements(self, element_iterator, context):
        try:
#            assert 'CLIXO' in os.environ

            params = {'similarity' : 'Similarity',
                      'name' : 'Data-Driven Ontology',
                      'dt_thresh': -100000,
                      'max_time': 100000}
            params.update(config.params)
            input_G, params, errors = self.read_element_stream(element_iterator, params)
            
            # Required parameters
            for x in ['ndex_uuid', 'alpha', 'beta', 'similarity', 'ndex_user', 'ndex_pass', 'ndex_server', 'output_fmt']:
                assert params.has_key(x)
            SIMILARITY = params['similarity']

            print 'Parameters:'
            print params

            ###############################

            if isinstance(params['ndex_uuid'], (str, unicode)):
                # Read graph using NDEx client
                input_G = NdexGraph(server=params['ndex_server'],
                                    username=params['ndex_user'],
                                    password=params['ndex_pass'],
                                    uuid=params['ndex_uuid'])
                assert len(input_G.nodes()) > 0
                input_G = NdexGraph_to_nx(input_G)

            gene_attr = nx_nodes_to_pandas(input_G)

            graph_df = nx_edges_to_pandas(input_G)
            #feature_columns = graph_df.columns.values.tolist()
            feature_columns = [SIMILARITY]

            graph_df.index.rename(['Gene1', 'Gene2'], inplace=True)
            graph_df.reset_index(inplace=True)
            graph_df[SIMILARITY] = graph_df[SIMILARITY].astype(np.float64)
            mat = pivot_square(graph_df, 'Gene1', 'Gene2', SIMILARITY)
            arr, arr_genes, arr_genes_index = mat.values, mat.index.values, pd.Series(mat.index).to_dict()
            assert arr.flags['C_CONTIGUOUS']

#            'SIMILARITY'

            # print 'errors:', [type(x) for x in errors]
            # for x in errors:
            #     print x

            errors = [e for e in errors if e.message != "Error decoding token from stream, EOF"]
            if len(errors) == 0:
#            if True:
                ## Run CLIXO
                ontology = Ontology.run_clixo(graph_df[['Gene1', 'Gene2', SIMILARITY]],
                                              params['alpha'],
                                              params['beta'],
                                              params['dt_thresh'],
                                              params['max_time'],
                                              params['clixo_folder'])

                print graph_df.head()
                print feature_columns

                # Create an NDEX network for every term's subnetwork
                term_2_url = ontology.upload_subnetworks_2_ndex(graph_df,
                                                                feature_columns,
                                                                params['ndex_server'], params['ndex_user'], params['ndex_pass'],
                                                                params['name'],
                                                                propagate=True)
                term_2_uuid = {t : parse_ndex_uuid(url) for t, url in term_2_url.items()}

                if params['output_fmt']=='cx':
                    # Use CXmate to stream to NDEX
                    for elt in self.stream_ontology(ontology, term_sizes, term_2_uuid, tree_edges):
                        yield elt

                elif params['output_fmt']=='ndex':
                    ontology_ndex = ontology.format_for_hierarchical_viewer(term_2_uuid)
                    ontology_ndex = nx_to_NdexGraph(ontology_ndex)
                    ontology_ndex.set_name(params['name'])
                    ontology_url = ontology_ndex.upload_to(params['ndex_server'], params['ndex_user'], params['ndex_pass'])
                    print 'ontology_url:', ontology_url

                    description = 'Data-driven ontology created by CLIXO (parameters: alpha=%s, beta=%s).' % (params['alpha'], params['beta'])
                    description += ' Created from similarity network at %s/%s' % (params['ndex_server'], params['ndex_uuid'])
                    ontology_ndex.set_network_attribute('Description', description)

                    for elt in yield_ndex(ontology_url):
                        yield elt
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
                yield self.create_edgeAttribute(edge_id, 'Relation', 'Gene-Term Annotation')
                if (g, t) in tree_edges:
                    yield self.create_edgeAttribute(edge_id, 'Is_Tree_Edge', 'Tree')
                else:
                    yield self.create_edgeAttribute(edge_id, 'Is_Tree_Edge', 'Not_Tree')
                edge_id += 1
        for p, c_list in ontology.term_2_terms.iteritems():        
            for c in c_list:
                yield self.create_edge(edge_id, node_2_id[c], node_2_id[p])
                yield self.create_edgeAttribute(edge_id, 'Relation', 'Child-Parent Hierarchical Relation')
                if (c,p) in tree_edges:
                    yield self.create_edgeAttribute(edge_id, 'Is_Tree_Edge', 'Tree')
                else:
                    yield self.create_edgeAttribute(edge_id, 'Is_Tree_Edge', 'Not_Tree')
                edge_id += 1

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
