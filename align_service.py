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

from align_hierarchies import align_hierarchies
from Ontology import Ontology, networkx_to_NdexGraph, NdexGraph_to_Ontology
from utilities import pivot_2_square

_ONE_DAY_IN_SECONDS = 60 * 60 * 24

class CyServiceServicer(cx_pb2_grpc.CyServiceServicer):

    def StreamElements(self, element_iterator, context):
        try:
            assert 'ALIGN_ONTOLOGY' in os.environ

            align_params = {'ndex_uuid' : None,
                            'ndex_server' : 'http://public.ndexbio.org/',
                            'name' : 'Data-Driven Ontology',
                            'ont1_ndex_uuid' : None,
                            'ont2_ndex_uuid' : None,
                            'calculateFDRs_cmd': os.path.join(os.getenv('ALIGN_ONTOLOGY'), 'calculateFDRs'),
                            'iterations' : 3,
                            'threads' : 4,
                            'verbose': True}

            input_G, align_params, errors = self.read_element_stream(element_iterator, align_params)

            print 'Parameters:'
            print align_params

            ###############################

            if isinstance(align_params['ont1_ndex_uuid'], (str, unicode)) and \
               isinstance(align_params['ont2_ndex_uuid'], (str, unicode)):
                ## Read graph using NDEx client

                print 'Reading hierarchy 1'
                hier1 = NdexGraph(server=align_params['ndex_server'],
                                  username=align_params['ndex_user'],
                                  password=align_params['ndex_pass'],
                                  uuid=align_params['ont1_ndex_uuid'])

                print 'Reading hierarchy 2'
                hier2 = NdexGraph(server=align_params['ndex_server'],
                                  username=align_params['ndex_user'],
                                  password=align_params['ndex_pass'],
                                  uuid=align_params['ont2_ndex_uuid'])
            else:
                raise Exception()

#            if len(errors) == 0:
            if True:
                print 'Creating list representation of hierarchies'
                tr = lambda x : 'gene' if x=='Gene-Term Annotation' else 'default'

                hier1_ont = NdexGraph_to_Ontology(hier1, gene_term='Gene-Term Annotation')
                hier2_ont = NdexGraph_to_Ontology(hier2, gene_term='Gene-Term Annotation')
                print 'Summary of hier1_ont:'
                print hier1_ont.summary()
                print 'Summary of hier2_ont:'
                print hier2_ont.summary()

                common_genes = set(hier1_ont.genes) & set(hier2_ont.genes)

                hier1_ont.delete_genes(set(hier1_ont.genes) - common_genes)
                hier1_collapsed = hier1_ont.collapse_ontology(method='mhkramer')
                print 'Summary of hier1_collapsed:'
                print hier1_collapsed.summary()

                hier2_ont.delete_genes(set(hier2_ont.genes) - common_genes)
                hier2_collapsed = hier2_ont.collapse_ontology(method='mhkramer')
                print 'Summary of hier2_collapsed:'
                print hier2_collapsed.summary()

                assert len(hier1_collapsed.terms) < 3000, len(hier1_collapsed.terms)
                assert len(hier2_collapsed.terms) < 3000, len(hier2_collapsed.terms)

                hier1_collapsed_nx = hier1_collapsed.to_networkx()
                hier2_collapsed_nx = hier2_collapsed.to_networkx()
                
                def nx_to_edges(G):
                    return [(G.node[v]['name'], G.node[u]['name'], tr(attr['EdgeType'])) \
                               for u, v, attr in G.edges_iter(data=True)]

                hier1_edges = nx_to_edges(hier1_collapsed_nx)
                hier2_edges = nx_to_edges(hier2_collapsed_nx)
                
#                with tempfile.NamedTemporaryFile('w', delete=True) as output_file:
                with tempfile.NamedTemporaryFile('w', delete=False) as output_file:
                    output = output_file.name
                    print 'output:', output

                    print 'Aligning hierarchies'
                    align_hierarchies(hier1_edges, hier2_edges,
                                      output,
                                      align_params['iterations'],
                                      align_params['threads'])
                
                    alignment = pd.read_csv(output,
                                            names=['Term_1', 'Term_2', 'Similarity', 'FDR', 'Term_1_Size'],
                                            header=None,
                                            sep='\t')

                # import cPickle
                # f = open('/tmp/tmp.pickle', 'wb')
                # print 'Writing pickle'
                # cPickle.dump({'alignment':alignment, 'hier1':hier1}, f)
                # f.close()

                # hier1_names = nx.get_node_attributes(hier1, 'name')
                # assert len(set(hier1_names.values())) == len(hier1_names.values())
                # hier1_names_idx = {b : a for a, b in hier1_names.items()}
                # hier1_nodes = hier1.nodes(data=True)

                print 'One-to-one term alignments:', alignment.shape[0]
                print alignment.iloc[:30,:]

#                G.get_node_ids("Brown","Color")
#                term1_ids = G.get_node_ids(*alignment['Term_1'])

                # for col in ['Term_2', 'Similarity', 'FDR']:
                #     nx.set_node_attributes(hier1, col, {hier1_names_idx[t1] : v for t1, v in zip(alignment['Term_1'], alignment[col])})

                for index, row in alignment.iterrows():
#                    term1 = hier1_nodes[hier1_names_idx[row['Term_1']]][1]
                    term1 = hier1.node[hier1.get_node_ids(row['Term_1'])[0]]
                    print 'Term1:', term1
                    term1['Aligned_Term'] = row['Term_2']
                    term1['Aligned_Similarity'] = row['Similarity']
                    term1['Aligned_FDR'] = row['FDR']

                ndex_url = hier1.upload_to(align_params['ndex_server'], align_params['ndex_user'], align_params['ndex_pass'])
                ndex_uuid = ndex_url.split('v2/network/')[1]
                print 'ndex_uuid:', ndex_uuid

                element = cx_pb2.Element()
                netAttr = element.networkAttribute
                netAttr.name = 'ndex_uuid'
                netAttr.value = ndex_uuid
                yield element
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
