import cxmate
import sys

sys.path.insert(0, '/cellar/users/mikeyu/DeepTranslate/ddot')
import ddot

class MyService(cxmate.Service):

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

    def process(self, params, input_stream):
        """process is a required method, if it's not implemented, cxmate.service will throw an error
        this process implementation will echo the received network back to the sender

        :param input_stream: a python generator that returns CX
        elements :returns: a python generator that returns CX elements
        """

        default = {
            'similarity' : 'similarity',
            'agg' : 'mean',
            'min_sim' : -np.inf,
            'filter_perc' : None,
            'seed_perc' : None,
            'agg_perc' : 0.5,
            'expand_size' : None,
            'include_seed' : True,
            'figure' : False
        }

        # params = self.format_params(params)
        default.update(params)
        params = default
        print 'Parameters:', params
        
        print 'Reading network...'
        network = cxmate.Adapter.to_networkx(input_stream)
        network = network[0]
        network = ddot.nx_to_NdexGraph(network)
        print 'Finished reading network.'

        print 'Type:', type(network)
        
        similarity = ddot.nx_edges_to_pandas(network, [params['similarity']])
        similarity = melt_square(similarity.reset_index(),
                                 columns=['Node1', 'Node2'],
                                 similarity=params['similarity'])
        print similarity.head()
        
        0 / asdf

        expand, expand_idx, sim_2_seed, fig = ddot.expand_seed(
            params['seed'],
            similarity.values,
            similarity.index.values,
            agg=params['agg'],
            min_sim=params
        )
            
            
        # edges = 0
        # nodes = 0
        # for elt in input_stream:
        #     print elt
        #     elt_type = elt.WhichOneof('element')
        #     if elt_type == 'edge':
        #         edges += 1
        #     elif elt_type == 'node':
        #         nodes += 1
        # print 'Edges:', edges
        # print 'Nodes:', nodes
        
#        print 'Nodes/Edges:', len(network.nodes()), len(network.edges())
        
        yield
        
if __name__=='__main__':
    myService = MyService()
    print 'Starting run'
    myService.run('localhost:8080') #run starts the service listening for requests from cxMate
    print 'Running'
