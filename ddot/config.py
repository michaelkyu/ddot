from __future__ import absolute_import

import ddot
import os, inspect
from ndex.networkn import NdexGraph

clixo = '/cellar/users/mikeyu/clixo'
alignOntology = '/cellar/users/mikeyu/alignOntology'

ndex_server = 'http://test.ndexbio.org'
ndex_user = 'scratch'
ndex_pass = 'scratch'

top_level = os.path.dirname(os.path.abspath(inspect.getfile(ddot)))
import json
with open(os.path.join(top_level, 'ontology_style.cx')) as f:
    ontology_style = NdexGraph(json.load(f))

with open(os.path.join(top_level, 'passthrough_style.cx')) as f:
    passthrough_style = NdexGraph(json.load(f))

GO_HUMAN_URL = 'http://dev2.ndexbio.org/v2/network/36e6e5b2-c2aa-11e7-9379-0660b7976219'
PHAROS_URL = 'http://dev2.ndexbio.org/v2/network/a94f1c0f-789a-11e7-a1d1-0660b7976219'
MONARCH_DISEASE_GENE_URL = 'http://dev2.ndexbio.org/v2/network/07749a5f-7956-11e7-a1d1-0660b7976219'
