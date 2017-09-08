import cx_pb2
from ddot.utils import parse_ndex_uuid

def yield_ndex(ndex_url):
    element = cx_pb2.Element()
    netAttr = element.networkAttribute
    netAttr.name = 'ndex_uuid'
    netAttr.value = parse_ndex_uuid(ndex_url)
    yield element

    element = cx_pb2.Element()
    netAttr = element.networkAttribute
    netAttr.name = 'ndex_url'
    netAttr.value = ndex_url
    yield element


def required_params(params, required):
    """Check that required parameters were set"""

    for x in required:
        assert params.has_key(x), 'Parameter %s not set' % x

def cast_params(params, cast):
    for p, c in cast:
        try:
            params[p] = c(params[p])
        except:
            pass
