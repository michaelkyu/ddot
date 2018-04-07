from __future__ import absolute_import

import os, inspect, io
import json

from ndex.networkn import NdexGraph

import ddot

###########################################
# Default NDEx server, username, password #
###########################################

ndex_server = 'http://test.ndexbio.org'
ndex_user = 'scratch'
ndex_pass = 'scratch'

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

GO_HUMAN_URL = 'http://dev2.ndexbio.org/v2/network/3305f6f3-10f9-11e8-84e4-0660b7976219'
PHAROS_URL = 'http://dev2.ndexbio.org/v2/network/a94f1c0f-789a-11e7-a1d1-0660b7976219'
MONARCH_DISEASE_GENE_URL = 'http://dev2.ndexbio.org/v2/network/07749a5f-7956-11e7-a1d1-0660b7976219'
