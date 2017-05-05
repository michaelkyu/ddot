import sklearn
from ndex.networkn import NdexGraph

def rf_integrate(feature_list, feature_fmt=None, labels):
    assert feature_fmt in ['UUID', 'numpy']
    
    if feature_fmt == 'UUID':
        network_list = [NdexGraph(f) for f in feature_list]
        all_edges = set([net.edges() for net in network_list])
        
        # Integrate everything into one feature matrix
        features = None
    elif feature_fmt == 'numpy':
        features = np.hstack([np.loadtxt(f) for f in feature_list])

    if label_fmt == 'UUID':
        NdexGraph(labels)
    elif label_fmt == 'NDEX':
        labels = np.loadtxt(label)

    rf = sklearn.ensemble.RandomForestRegressor()
    rf.fit(features, labels)
    pred = rf.predict(labels)

    return pred

