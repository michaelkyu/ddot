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
from Ontology import Ontology, networkx_to_NdexGraph
from utilities import pivot_2_square

_ONE_DAY_IN_SECONDS = 60 * 60 * 24

class CyServiceServicer(cx_pb2_grpc.CyServiceServicer):

    def StreamElements(self, element_iterator, context):
        try:
            align_params = {'ndex_uuid' : None,
                            'ndex_server' : 'http://public.ndexbio.org/',
                            'name' : 'Data-Driven Ontology',
                            'ont1_ndex_uuid' : None,
                            'ont2_ndex_uuid' : None,
                            'verbose': True}

            input_G, align_params, errors = self.read_element_stream(element_iterator, align_params)

            ###############################

            if isinstance(align_params['ont1_ndex_uuid'], (str, unicode)) and \
               isinstance(align_params['ont2_ndex_uuid'], (str, unicode)):
                # Read graph using NDEx client
                hier1 = NdexGraph(server=align_params['ndex_server'],
                                  username=align_params['ndex_user'],
                                  password=align_params['ndex_pass'],
                                  uuid=align_params['ont1_ndex_uuid'])

                hier2 = NdexGraph(server=align_params['ndex_server'],
                                  username=align_params['ndex_user'],
                                  password=align_params['ndex_pass'],
                                  uuid=align_params['ont2_ndex_uuid'])                
            else:
                raise Exception()

            if len(errors) == 0:               
                graph = 

                hier1_edges = [(hier1.node[u]['name'], hier1.node[v]['name'], attr['EdgeType']) for u, v, attr in hier1.edges_iter(data=True)]
                hier2_edges = [(hier2.node[u]['name'], hier2.node[v]['name'], attr['EdgeType']) for u, v, attr in hier2.edges_iter(data=True)]
                align_hierarchies(hier1, hier2,
                                  align_params['output'],
                                  align_params['iterations'],
                                  align_params['threads'])
                
                alignment = pd.read_csv(align_params['output'],
                                        columns=['Term_1', 'Term_2', 'Similarity', 'FDR', 'Term_1_Size'],
                                        header=None,
                                        sep='\t')
                for t in hier1.nodes:
                    hier1.node[t]['Aligned_Term'] = alignment[t]['Term_2']
                    hier1.node[t]['Aligned_Similarity'] = alignment[t]['Term_2']
                    hier1.node[t]['Aligned_FDR'] = alignment[t]['Term_2']
                #yield hier1

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

    def stream_ontology(self, ontology, term_sizes, term_2_uuid):
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
                edge_id += 1
        for p, c_list in ontology.term_2_terms.iteritems():                        
            for c in c_list:
                yield self.create_edge(edge_id, node_2_id[c], node_2_id[p])
                yield self.create_edgeAttribute(edge_id, 'EdgeType', 'Child-Parent Hierarchical Relation')
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
        error.code = 'cy://align-hierarchies/' + str(status)
        error.message = message
        error.link = 'http://align-hierarchies'
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
