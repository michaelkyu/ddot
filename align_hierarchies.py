import sys, shutil, os, argparse, tempfile
from subprocess import Popen, PIPE, STDOUT
import tempfile
from utilities import time_print

def align_hierarchies(hier1, hier2,
                      output,
                      iterations, threads,
                      calculateFDRs_cmd='/cellar/users/mikeyu/alignOntology/calculateFDRs'):

    if not isinstance(hier1, (str, unicode)):        
        # Write to file
        with open('/tmp/tmp1.txt', 'w') as f:
#        with tempfile.NamedTemporaryFile('w', delete=False) as f:
            f.write('\n'.join(['\t'.join([str(x[0]),str(x[1]),str(x[2])]) for x in hier1]) + '\n')
        hier1 = f.name
    else:
        assert os.path.exists(hier1)

    if not isinstance(hier2, (str, unicode)):
        with open('/tmp/tmp2.txt', 'w') as g:
#        with tempfile.NamedTemporaryFile('w', delete=False) as g:
            g.write('\n'.join(['\t'.join([str(x[0]),str(x[1]),str(x[2])]) for x in hier2]) + '\n')
        hier2 = g.name
    else:
        assert os.path.exists(hier2)

    output_dir = tempfile.mkdtemp(prefix='tmp')

    cmd = '{5} {0} {1} 0.05 criss_cross {2} {3} {4} gene'.format(
           hier1, hier2, output_dir, iterations, threads, calculateFDRs_cmd)
    time_print(cmd)

    try:
        p = Popen(cmd, shell=True)
        p.wait()

        # Five columns in the output file
        # 1) Term from first ontology (Kramer calls it the "computed" ontology)
        # 2) Term from second ontology (Kramer calls it the "reference" ontology)
        # 3) Similarity value
        # 4) FDR
        # 5) Size of the term in the first ontology
        shutil.copy(os.path.join(output_dir, 'alignments_FDR_0.1_t_0.1'), output)
    finally:
        # if os.path.isdir(output + '_dir'):
        #     shutil.rmtree(output + '_dir')
#        shutil.move(output_dir, output + '_dir')

        if os.path.isdir(output_dir):
            shutil.rmtree(output_dir)
        
        if p.poll() is None:
            if verbose: time_print('Killing alignment process %s. Output: %s' % (p.pid, output))
            p.kill()  # Kill the process

    return

if __name__=='__main__':
    parser = argparse.ArgumentParser('Align two hierarchies')
    parser.add_argument('--hier1', required=True)
    parser.add_argument('--hier2', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--iterations', type=int, default=100)
    parser.add_argument('--threads', type=int, default=30)
    args = parser.parse_args()

    align_hierarchies(args.hier1, args.hier2, args.output, args.iterations, args.threads)
