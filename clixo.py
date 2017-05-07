from time import sleep
import os, argparse, time, tempfile, shutil
import sys
import time

import numpy as np
from itertools import combinations
from subprocess import Popen, PIPE, STDOUT
from utilities import time_print

def run_clixo(graph, alpha, beta, dt_thresh, max_time,
              output, clixo_folder, output_log=None, verbose=True):

    if not isinstance(graph, str):
        # Write graph into a temporary file.
        # Assumes that <graph> is a list of 3-tuples (parent, child, score)
        with tempfile.NamedTemporaryFile('w', delete=False) as f:
            f.write('\n'.join(['\t'.join([str(x[0]),str(x[1]),str(x[2])]) for x in graph]) + '\n')

        graph = f.name
        delete_file = f.name
        print 'Temporary file for graph:', f.name
    else:
        delete_file = False

    if verbose: time_print('\t'.join(map(str, [graph, alpha, beta, dt_thresh])))

    # '/cellar/users/mikeyu/mhk7-clixo_0.3-cec3674'
    clixo_cmd = os.path.join(clixo_folder, 'clixo')
    extract_cmd = os.path.join(clixo_folder, 'extractOnt')

    if not isinstance(output_log, str):
        output_log_file = tempfile.NamedTemporaryFile('w', delete=True)
        output_log = output_log_file.name
        delete_output_log = True
    else:
        # print 'output_log', output_log
        # print 'output_log dirname:', os.path.dirname(output_log)
        # print 'output_log dirname exists:', os.path.isdir(os.path.dirname(output_log))
        assert os.path.isdir(os.path.dirname(output_log))
        delete_output_log = False

    # For timestamping everyline: awk '{ print strftime("%Y-%m-%d %H:%M:%S"), $0; fflush(); }'
    cmd = """{0} {1} {2} {3} | awk""".format(clixo_cmd, graph, alpha, beta) + \
          """ '{if ( $1 ~ /^#/ ) {print "\#", strftime("%Y-%m-%d %H:%M:%S"), $0 ; fflush() } else {print $0}}'""" + \
          """ | tee {}""".format(output_log)
    print >>sys.stderr, cmd

    p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, bufsize=1)

    curr_dt = None
    start = time.time()

    # Asynchronous readout of the command's output
    while p.poll() is None:

        # Break if passed the maximum processing time
        if time.time() - start > max_time:
            if verbose: time_print('Killing process %s (OUT OF TIME). Current dt: %s: Output: %s' % (p.pid, curr_dt, output_log))
            break

        line = p.stdout.readline()

        # Remove newline character
        line = line[:-1]

        # Break if the dt_threshold has been met
        if '# dt: ' in line:
            curr_dt = float(line.split('# dt: ')[1])
            if curr_dt < dt_thresh:
                if verbose: time_print('Killing process %s (BEYOND MIN THRESHOLD). Current dt: %s, dt_thresh: %s' % (p.pid, curr_dt, dt_thresh))
                break

        # If line was empty, then sleep a bit
        if line=='':  time.sleep(0.1)

    if p.poll() is None:
        if verbose: time_print('Killing process %s. Output: %s' % (p.pid, output_log))
        p.kill()  # Kill the process

        # Extract ontology with extractOnt
        p_ext = Popen('%s %s 0 0 %s' % (extract_cmd, output_log, output),
                      shell=True, stdout=PIPE, stderr=STDOUT)
        p_ext.communicate()
    else:
        if verbose: time_print('Extracting by grep -v #')

        # Extract ontology with grep -v '#'
        p_ext = Popen("grep -v '#' %s | grep -v '@' > %s" % (output_log, output), shell=True)
        p_ext.communicate()

    if delete_output_log:
        output_log_file.close()

    if delete_file and os.path.isfile(delete_file):
        os.remove(delete_file)

    time_print('Elapsed time (sec): %s' % (time.time() - start))

    return output, output_log

# if __name__=='__main__':
#     parser = argparse.ArgumentParser('Convert contact maps into a clixo readable format')
#     parser.add_argument('--alpha', type=float, nargs='*', help='List of alpha parameters')
#     parser.add_argument('--beta', type=float, nargs='*', help='List of beta parameters')
#     parser.add_argument('--graphs', type=str, nargs='*', help='List of pairwise similarity files')
#     parser.add_argument('--clixo_folder', type=str, default='/cellar/users/mikeyu/mhk7-clixo_0.3-cec3674')
#     parser.add_argument('--output_suffix', type=str, nargs='?', const='', default='')
#     parser.add_argument('--dt_thresh', type=float,
#                         help="Terminate CliXO when dt falls below this threshold. If not specified, then run to completion.")
#     parser.add_argument('--max_time', type=float,
#                         help="Max allowable time before CliXO is terminated.  If not specified, then run to completion")
#     parser.add_argument('--n_jobs', type=int, help="Number of parallel jobs")
#     parser.add_argument('--iteration_chunk', type=int, default=1,
#                         help="In combination with the arg <iteration>, this divide the parameter settings into equal-sized chunks")
#     parser.add_argument('--iteration', type=int, default=None,
#                         help="In combination with the arg <iteration_chunk>, this is the 0-based index of which chunk of parameter settings to use")
#     parser.add_argument('--warm_start', action='store_true', help="Skip parameters if the output file already exists")
#     parser.add_argument('--output', default=None)
#     args = parser.parse_args()

#     pprint.pprint(vars(args))

#     param_list = [(g, a, b, args.dt_thresh, args.max_time, args.warm_start, args.output, args.clixo_folder, args.output_suffix) \
#                   for g in args.graphs for a in args.alpha for b in args.beta]
#     param_list = [x + (i, ) for i, x in enumerate(param_list)]

# #def run_clixo(graph, alpha, beta, dt_thresh, max_time, warm_start, iteration, verbose=True):

#     # If an iteration index is specified, then run only one setting of parameters
#     if args.iteration is not None:
#         param_list = param_list[args.iteration_chunk*args.iteration: args.iteration_chunk*(args.iteration+1)]

#     try:
#         pool = Pool(args.n_jobs, maxtasksperchild=1)
#         start = time.time()
#         val_list = [run_clixo_star(x) for x in param_list]        
# #        time_print('Running pool')        
# #        val_list = pool.map(run_clixo_star, param_list, chunksize=1)
#     finally:
#         pool.close()
#         print 'Elapsed time (sec)', time.time() - start
