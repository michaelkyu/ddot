import grpc
import cx_pb2
import cx_pb2_grpc

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

def parse_ndex_uuid(ndex_url):
    # print 'ndex_url:', ndex_url
    # print 'ndex uuid:', ndex_url.split('v2/network/')[1]
    return ndex_url.split('v2/network/')[1]
