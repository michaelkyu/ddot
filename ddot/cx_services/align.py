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

from ddot import Ontology, align_hierarchies
from ddot.utils import update_nx_with_alignment
from ddot.cx_services.cx_utils import yield_ndex, required_params, cast_params
from ddot.config import default_params

_ONE_DAY_IN_SECONDS = 60 * 60 * 24

verbose = True

class CyServiceServicer(cx_pb2_grpc.CyServiceServicer):

    def format_params(self, params):
        required = [
            'ndex_user',
            'ndex_pass',
            'ndex_server',
            'ont1_ndex_uuid',
            'ont2_ndex_uuid',
        ]
        required_params(params, required)

        cast = [
            ('iterations', int),
            ('threads', int),
            ('ont1_ndex_uuid', str),
            ('ont2_ndex_uuid', str)
        ]
        cast_params(params, cast)

    def StreamElements(self, element_iterator, context):
        try:
            params = {
                'name' : 'Data-Driven Ontology',
                'ont1_ndex_uuid' : None,
                'ont2_ndex_uuid' : None,
                'iterations' : 3,
                'threads' : 4
            }
            params.update(default_params)
            G, params, errors = self.read_element_stream(element_iterator, params)
            self.format_params(params)

            if verbose:
                print('Parameters:', params)

            ## Read graphs using NDEx client
            hier1, hier2 = self.read_hierarchies(params)

            if True:
                hier1_ont = Ontology.from_NdexGraph(hier1)
                hier2_ont = Ontology.from_NdexGraph(hier2)
                print('Summary of hier1_ont:', hier1_ont.summary())
                print('Summary of hier2_ont:', hier2_ont.summary())

                hier1_collapsed, hier2_collapsed = Ontology.mutual_collapse(hier1_ont, hier2_ont, verbose=True)
                assert len(hier1_collapsed.terms) < 3000, len(hier1_collapsed.terms)
                assert len(hier2_collapsed.terms) < 3000, len(hier2_collapsed.terms)

                if verbose:
                    print 'Aligning hierarchies'

                alignment = align_hierarchies(
                    hier1_collapsed,
                    hier2_collapsed,
                    params['iterations'],
                    params['threads'])
                
                if verbose:
                    print('One-to-one term alignments:', alignment.shape[0])
                    print(alignment.iloc[:30,:])

                update_nx_with_alignment(hier1, alignment)
                ont_url = hier1.upload_to(params['ndex_server'],
                                           params['ndex_user'],
                                           params['ndex_pass'])
                if verbose:
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

    def read_hierarchies(self, params):
        # Read hierarchy 1
        hier1 = NdexGraph(server=params['ndex_server'],
                          username=params['ndex_user'],
                          password=params['ndex_pass'],
                          uuid=params['ont1_ndex_uuid'])

        # Read hierarchy 2
        hier2 = NdexGraph(server=params['ndex_server'],
                          username=params['ndex_user'],
                          password=params['ndex_pass'],
                          uuid=params['ont2_ndex_uuid'])

        return hier1, hier2

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
                yield self.create_edgeAttribute(edge_id, 'Relation', 'Gene-Term Annotation')
                edge_id += 1
        for p, c_list in ontology.term_2_terms.iteritems():                        
            for c in c_list:
                yield self.create_edge(edge_id, node_2_id[c], node_2_id[p])
                yield self.create_edgeAttribute(edge_id, 'Relation', 'Child-Parent Hierarchical Relation')
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
    server.add_insecure_port('0.0.0.0:8081')
    server.start()
    try:
        while True:
            time.sleep(_ONE_DAY_IN_SECONDS)
    except KeyboardInterrupt:
        server.stop(0)

if __name__ == '__main__':
    logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
    log_info("Listening for requests on '0.0.0.0:8081'")
    serve()
