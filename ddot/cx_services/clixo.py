import ndex.client as nc
from ndex.networkn import NdexGraph

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

from ddot import Ontology
from ddot.utils import ndex_to_sim_matrix, nx_nodes_to_pandas
from ddot.cx_services.cx_utils import yield_ndex, required_params, cast_params
from ddot.config import default_params

_ONE_DAY_IN_SECONDS = 60 * 60 * 24

verbose = True

class CyServiceServicer(cx_pb2_grpc.CyServiceServicer):

    def format_params(self, params):
        required = [
            'ndex_uuid',
            'ndex_user',
            'ndex_pass',
            'ndex_server',
            'alpha',
            'beta',
            'similarity',
            'input_fmt',
            'output_fmt'
        ]
        required_params(params, required)

        cast = [
            ('alpha', float),
            ('beta', float),
            ('min_dt', float),
            ('timeout', int),
            ('all_features', bool)
        ]
        cast_params(params, cast)

        if params['features']:
            params['features'] = params['features'].split(',')

        assert params['input_fmt'] in ['cx', 'cx_matrix', 'ndex', 'ndex_matrix']
        assert params['output_fmt'] in ['cx', 'cx_matrix', 'ndex', 'ndex_matrix']

    def StreamElements(self, element_iterator, context):
        try:
            params = {
                'similarity' : 'Similarity',
                'name' : 'Data-Driven Ontology',
                'min_dt' : -100000,
                'features' : '',
                'all_features' : False,
                'input_fmt' : 'cx',
                'timeout' : 100,
                'ndex_uuid' : None
            }
            params.update(default_params)
            G, params, errors = self.read_element_stream(element_iterator, params)
            self.format_params(params)

            if verbose:
                print('Parameters', params)

            input_fmt = params['input_fmt'].replace('ndex', 'cx')

            G_df, nodes_attr = ndex_to_sim_matrix(
                params['ndex_uuid'],
                params['ndex_server'],
                params['ndex_user'],
                params['ndex_pass'],
                similarity=params['similarity'],
                input_fmt=input_fmt,
                output_fmt='sparse')
            
            print 'Similarity'
            print G_df.head()

            features = [params['similarity']]
            if params['all_features']:
                features.extend(G_df.columns.values.tolist())

            gene_names = np.unique(np.append(
                G_df['Gene1'].values,
                G_df['Gene2'].values))
            if verbose:
                print 'Gene names:', gene_names

            if params['features']:
                to_concat = []

                for f in params['features']:
                    if f in G_df.columns:
                        in_G.append(f)
                    else:
                        tmp = ndex_to_sim_matrix(
                            params['ndex_uuid'],
                            params['ndex_server'],
                            params['ndex_user'],
                            params['ndex_pass'],
                            similarity=params['similarity'],
                            input_fmt=input_fmt,
                            output_fmt='sparse',
                            subset=None)
                        tmp.set_index(['Gene1', 'Gene2'], inplace=True)
                        to_concat.append(tmp)
                to_concat.append(G_df)
                G_df = pd.concat(to_concat, axis=1)

            if verbose:
                print 'Features:', features

            errors = [e for e in errors if e.message != "Error decoding token from stream, EOF"]
            if len(errors) == 0:

                ## Run CLIXO
                graph = G_df[['Gene1', 'Gene2', params['similarity']]]
                ont = Ontology.run_clixo(
                    graph,
                    params['alpha'],
                    params['beta'],
                    min_dt=params['min_dt'],
                    timeout=params['timeout']
                )                

                # Create an NDEX network for every term's subnetwork
                term_2_uuid = ont.upload_subnets_ndex(
                    G_df,
                    features,
                    params['ndex_server'],
                    params['ndex_user'],
                    params['ndex_pass'],
                    params['name'],
                    propagate=True
                )

                if params['output_fmt']=='cx':
                    # Use CXmate to stream to NDEX
                    for elt in self.stream_ontology(ont, term_sizes, term_2_uuid, tree_edges):
                        yield elt

                elif params['output_fmt']=='ndex':
                    description = (
                        'Data-driven ontology created by CLIXO '
                        '(parameters: alpha={alpha}, beta={beta}). '
                        'Created from similarity network '
                        'at {ndex_server}/{ndex_uuid}').format(**params)

                    # nx_nodes_to_pandas(G)

                    ont_ndex = ont.to_NdexGraph(
                        name=params['name'],
                        description=description,
                        term_2_uuid=term_2_uuid,
                        gene_attr=nodes_attr)

                    ont_url = ont_ndex.upload_to(
                        params['ndex_server'],
                        params['ndex_user'],
                        params['ndex_pass'])
                    print 'ontology_url:', ont_url

                    for elt in yield_ndex(ont_url):
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
