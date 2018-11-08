from __future__ import absolute_import

import os, inspect, io
import json

from ndex.networkn import NdexGraph

import ddot

###########################################
# Default NDEx server, username, password #
###########################################

ndex_server = 'http://test.ndexbio.org'
ndex_user, ndex_pass = 'scratch', 'scratch'

#########################
# Read CX Visual styles #
#########################

passthrough_style = None

def get_passthrough_style():
    global passthrough_style
    if passthrough_style is None:
        top_level = os.path.dirname(os.path.abspath(inspect.getfile(ddot)))
        with io.open(os.path.join(top_level, 'passthrough_style.cx')) as f:
            passthrough_style = NdexGraph(json.load(f))
    return passthrough_style        

##################################
# NDEx URLs for example networks #
##################################

# mikeyu_testacct
HUMAN_GENE_SIMILARITIES_URL = 'http://dev2.ndexbio.org/v2/2f8ab89d-e396-11e8-a260-0660b7976219'
GO_HUMAN_URL = 'http://dev2.ndexbio.org/v2/network/ccf12fb2-e399-11e8-a260-0660b7976219'
MONARCH_DISEASE_GENE_URL = 'http://dev2.ndexbio.org/v2/network/afecb56c-e398-11e8-a260-0660b7976219'
MONARCH_DISEASE_GENE_SLIM_URL = 'http://dev2.ndexbio.org/v2/network/85611be5-e39b-11e8-a260-0660b7976219'
