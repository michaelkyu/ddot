import ddot
import os, inspect
from ndex.networkn import NdexGraph

clixo = '/cellar/users/mikeyu/clixo'
alignOntology = '/cellar/users/mikeyu/alignOntology'

# ndex_server = 'http://test.ndexbio.org'
# ndex_user = 'scratch'
# ndex_pass = 'scratch'

ndex_server = 'http://test.ndexbio.org'
ndex_user = 'mikeyu_testacct4'
ndex_pass = 'GoHejVeg8'

top_level = os.path.dirname(os.path.abspath(inspect.getfile(ddot)))
import json
with open(os.path.join(top_level, 'ontology_style.cx')) as f:
    ontology_style = NdexGraph(json.load(f))

GO_HUMAN_URL = 'http://test.ndexbio.org/v2/network/503df69b-81f8-11e7-9743-0660b7976219'
#GO_HUMAN_URL = 'http://dev2.ndexbio.org/v2/network/af343c36-818c-11e7-9743-0660b7976219'
#GO_HUMAN_URL = 'http://dev2.ndexbio.org/v2/network/8bfa8318-55ed-11e7-a2e2-0660b7976219'

PHAROS_URL = 'http://dev2.ndexbio.org/v2/network/a94f1c0f-789a-11e7-a1d1-0660b7976219'
MONARCH_DISEASE_GENE_URL = 'http://dev2.ndexbio.org/v2/network/07749a5f-7956-11e7-a1d1-0660b7976219'
