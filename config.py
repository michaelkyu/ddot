import os

params = {
    'ndex_server' : 'http://test.ndexbio.org',
    'clixo_folder' : os.getenv('CLIXO'),
    'calculateFDRs_cmd': os.path.join(os.getenv('ALIGN_ONTOLOGY'), 'calculateFDRs'),
    'verbose' : False,
    'output_fmt' : 'ndex',
}
