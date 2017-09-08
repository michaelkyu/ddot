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

from ddot.utils import expand_seed, sim_matrix_to_NdexGraph, ndex_to_sim_matrix
from ddot.cx_services.cx_utils import yield_ndex, required_params, cast_params
from ddot.config import default_params

_ONE_DAY_IN_SECONDS = 60 * 60 * 24

verbose = True

class CyServiceServicer(cx_pb2_grpc.CyServiceServicer):

    def format_params(self, params):
        required = [
            'seed',
            'ndex_uuid',
            'ndex_user',
            'ndex_pass',
            'ndex_server',
            'input_fmt',
            'output_fmt'
        ]
        required_params(params, required)

        cast = [
            ('filter_percentile', float),
            ('seed_percentile', float),
            ('aggregation_percentile', float),
            ('min_similarity', float),
            ('expand_size', int)
        ]
        cast_params(params, cast)

        params['seed'] = [x.strip() for x in params['seed'].split(',')]

        assert params['input_fmt'] in ['cx', 'cx_matrix', 'ndex', 'ndex_matrix']
        assert params['output_fmt'] in ['cx', 'cx_matrix', 'ndex', 'ndex_matrix']

    # def expand_seed(self,
    #                 ndex_uuid,
    #                 name,
    #                 similarity,
    #                 aggregation,
    #                 filter_percentile,
    #                 seed_percentile,
    #                 aggregation_percentile,
    #                 expand_size,
    #                 min_similarity,
    #                 input_fmt,
    #                 output_fmt,
    #                 include_subnet,
    #                 ndex_server):

    #     if input_fmt=='ndex_matrix':
    #         sim, sim_names, sim_names_col = load_edgeMatrix(
    #             ndex_uuid, ndex_server, ndex_user, ndex_pass)
    #         assert sim_names == sim_names_col
    #     elif input_fmt=='ndex':
    #         sim, sim_names = read_ndex_similarity(params)

    #     expand, expand_idx, expand_sim = expand_seed(
    #         seed,
    #         sim,
    #         sim_names,
    #         agg=aggregation
    #         min_sim=min_similarity,
    #         filter_perc=filter_percentile,
    #         seed_perc=seed_percentile,
    #         agg_perc=aggregation_percentile,
    #         expand_size=expand_size
    #     )
    #     expand = list(expand)

    #     if verbose:
    #         print('Expanded set of genes:', zip(expand, expand_sim))

    #     if output_fmt in ['ndex_matrix', 'ndex']:
    #         # Take sub-array of expand-by-expand
    #         subnet = sim[:, expand_idx][expand_idx, :]
    #         subnet = np.array(subnet,
    #                           order='C',
    #                           dtype=np.float32,
    #                           copy=True)

    #         gene_attr = pd.DataFrame({similarity : expand_sim}, index=expand)

    #         G = NdexGraph_from_sim_matrix(
    #             subnet,
    #             expand,
    #             similarity,
    #             output_fmt,
    #             node_attr=node_attr
    #         )
    #         G.set_name(name)
    #         G.set_network_attribute('Description', self.get_description(params))

    #         ont_url = G.upload_to(ndex_server, ndex_user, ndex_pass)

    #         if verbose:
    #             print 'ontology_url:', ont_url

    #         for elt in yield_ndex(ont_url):
    #             yield elt

    #     elif output_fmt == 'cx':
    #         # Use CXmate to stream to NDEX
    #         for node_id, (gene, sim) in enumerate(zip(expand, expand_sim)):
    #             yield self.create_node(node_id, gene)
    #             yield self.create_nodeAttribute(node_id, similarity, unicode(sim))

    def StreamElements(self, element_iterator, context):
        try:
            params = {'ndex_uuid' : None,
                      'name' : 'Expanded Gene Set',
                      'similarity' : 'similarity',
                      'aggregation' : 'mean',
                      'filter_percentile' : None,
                      'seed_percentile' : None,
                      'aggregation_percentile' : 0,
                      'expand_size' : 100,
                      'min_similarity' : -float('inf'),
                      'input_fmt' : 'ndex_matrix',
                      'output_fmt' : 'ndex',
                      'include_subnet' : True,
            }
            params.update(default_params)
            G, params, errors = self.read_element_stream(element_iterator, params)
            self.format_params(params)
            
            if verbose:
                print('Parameters:', params)

            ###############################

            input_fmt = params['input_fmt'].replace('ndex', 'cx')

            sim, sim_names = ndex_to_sim_matrix(
                params['ndex_uuid'],
                params['ndex_server'],
                params['ndex_user'],
                params['ndex_pass'],
                input_fmt=input_fmt,
                output_fmt='matrix')

            expand, expand_idx, expand_sim = expand_seed(
                params['seed'],
                sim,
                sim_names,
                agg=params['aggregation'],
                min_sim=params['min_similarity'],
                filter_perc=params['filter_percentile'],
                seed_perc=params['seed_percentile'],
                agg_perc=params['aggregation_percentile'],
                expand_size=params['expand_size']
            )
            expand = list(expand)
            
            if verbose:
                print('Expanded set of genes:', zip(expand, expand_sim))

            if params['output_fmt'] in ['ndex_matrix', 'ndex']:
                # Take sub-array of expand-by-expand
                subnet = sim[:, expand_idx][expand_idx, :]
                subnet = np.array(subnet,
                                  order='C',
                                  dtype=np.float32,
                                  copy=True)

                node_attr = pd.DataFrame({params['similarity'] : expand_sim},
                                         index=expand)

                output_fmt = params['output_fmt'].replace('ndex', 'cx')

                G = sim_matrix_to_NdexGraph(
                    subnet,
                    expand,
                    params['similarity'],
                    output_fmt,
                    node_attr=node_attr
                )
                G.set_name(params['name'])
                G.set_network_attribute('Description', self.get_description(params))
                
                ont_url = G.upload_to(params['ndex_server'], params['ndex_user'], params['ndex_pass'])

                if verbose:
                    print 'ontology_url:', ont_url

                for elt in yield_ndex(ont_url):
                    yield elt

            elif params['output_fmt'] == 'cx':
                # Use CXmate to stream to NDEX
                for node_id, (gene, sim) in enumerate(zip(expand, expand_sim)):
                    yield self.create_node(node_id, gene)
                    yield self.create_nodeAttribute(node_id, params['similarity'], unicode(sim))
            
#            if len(errors) == 0:
            if True:
                pass
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

    def get_description(self, params):
        description = (
            'Genes related to a seed set. Created from\n'
            '\nSimilarity network: {ndex_server}/{ndex_uuid}'
            '\nSeed set: {seed}'
            '\naggregation: {aggregation}'
            '\nmin_similarity: {min_similarity}'
            '\nexpand_size: {expand_size}'
            '\nfilter_percentile: {filter_percentile}'
            '\nseed_percentile: {seed_percentile}').format(**params)
        return description

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
    server.add_insecure_port('0.0.0.0:8082')
    server.start()
    try:
        while True:
            time.sleep(_ONE_DAY_IN_SECONDS)
    except KeyboardInterrupt:
        server.stop(0)

if __name__ == '__main__':
    logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
    log_info("Listening for requests on '0.0.0.0:8082'")
    serve()
