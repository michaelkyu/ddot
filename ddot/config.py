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

# mikeyu (public server)
HUMAN_GENE_SIMILARITIES_URL = "http://public.ndexbio.org/v2/network/cbd481aa-f73c-11e8-aaa6-0ac135e8bacf"
GO_HUMAN_URL = 'http://public.ndexbio.org/v2/network/1a0917db-f739-11e8-aaa6-0ac135e8bacf'
MONARCH_DISEASE_GENE_URL = 'http://public.ndexbio.org/v2/network/bce041f4-f739-11e8-aaa6-0ac135e8bacf'
MONARCH_DISEASE_GENE_SLIM_URL = 'http://public.ndexbio.org/v2/network/aba9e837-f73b-11e8-aaa6-0ac135e8bacf'
FANGO_URL = "http://public.ndexbio.org/v2/network/0fb9fec3-f772-11e8-aaa6-0ac135e8bacf"
FANGO_DATA_URL = "http://public.ndexbio.org/#/network/84898ca0-fa98-11e8-8438-0ac135e8bacf"
