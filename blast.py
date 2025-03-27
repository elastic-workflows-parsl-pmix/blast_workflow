#!/usr/bin/env python
''' Simple parallelized blast mockup.
# The basic workflow is structured as follows :
#           Fastq_file
#              |
#         fast_q_split()
#             / \
#         [splitfiles ...]
#           |   |       |
#         blast() blast()  blast()
#           |   |       |
#          [bam_files ...]
#            \  |      /
#             merge()
#               |
#            check_merged()
'''
from os.path import abspath
import argparse
import random
import parsl
from parsl.app.app import python_app, bash_app
from parsl.data_provider.files import File
import os
import time
import sys
from parsl.config import Config
from parsl.providers import PMIxProvider
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SimplePMIxLauncher
from parsl.addresses import address_by_hostname

path = "/home/rbhattara/europar25/tests/test/blast/"

#Here we will use the bash split utility to do the task of splitting the fast_q_file
#into n_chunks. We take the abs_path to the input fast_q_file and return a bunch of split
#files whose names are already in the [outputs] array.
@bash_app
def fast_q_split(fast_q_file, n_chunks, outputs=[],
                 stdout=abspath('test/blast/split.out'), stderr=abspath('test/blast/split.err')):
    cmd_line = '''

    echo "Received file: {0}, splitting into :{1} pieces"
    lines=$(wc -l {0} | awk '{{print $1}}')
    split -l $(($lines / {1})) -d {0} /home/rbhattara/europar25/tests/test/blast/sequence_collagen.fasta.split.

    '''.format(fast_q_file, n_chunks)
    return cmd_line

@bash_app
def blast(inputs=[], outputs=[], stdout=abspath('test/blast/blast.out'), stderr=abspath('test/blast/blast.err')):
    cmd_line = '''

    prun --map-by node --bind-to none --display bind -np 1 /home/rbhattara/europar25/blast_workflow/mpiblastpfix {0} {1} /home/rbhattara/blast_sandbox/swissprot/swissprot

    '''.format(inputs[0].filepath, outputs[0].filepath)
    return cmd_line

# The merge file simply concatenated the array of files passed in through the inputs list
# into the single file specified though the outputs keyword_arg.
# Note: On the commandline, {inputs} expands to the python list of files expressed as a
# string, eg "['.../bwa.01', '.../bwa.02', ...] which does not make sense in bash.
# So we translate this string with tr to ('.../bwa.01' '.../bwa.02' ...)
@bash_app
def merge(inputs=[], outputs=[], stdout=abspath('test/blast/merge.out'), stderr=abspath('test/blast/merge.err')):
    cmd_line = '''

    echo "Merging {0} -> {1}"
    input_array=($(echo {0} | tr -d '[],' ))
    cat ${{input_array[*]}} &> {1}

    '''.format(inputs, outputs[0].filepath)
    return cmd_line


if __name__ == "__main__" :

    parser   = argparse.ArgumentParser()
    parser.add_argument("-s", "--split", default='10', help="Number of files to split the fastq file to")
    parser.add_argument("-f", "--fastq", default='mock.fastq', help="File path to the fastq file")
    parser.add_argument("-N", "--nnodes", default='6', help="Number of nodes")
    parser.add_argument("-L", "--nodelist", default='', help="List of nodes")
    parser.add_argument('-J', "--jobid", default=1, help="Job id from scheduler")
    parser.add_argument('-wt', "--wt", default='00:30:00', help="Walltime")
    parser.add_argument("-v", "--verbose", dest='verbose', action='store_true', help="Verbose output")
    args   = parser.parse_args()

    # Handle command-line options
    if args.verbose:
        parsl.set_stream_logger()

    n_chunks   = int(args.split)
    fastq_file = args.fastq
    n_nodes = int(args.nnodes)
    nodelist = args.nodelist
    job_id = int(args.jobid)
    walltime = args.wt
    # config goes here
    config = Config(
            executors=[
                HighThroughputExecutor(
                    label="blast_new",
                    address=address_by_hostname(),
                    prefetch_capacity=40,
                    heartbeat_period=15,
                    max_workers_per_node=n_nodes,
                    provider=PMIxProvider(
                        cmd_timeout=20,
                        cores_per_node=1,
                        nodes_per_block=n_nodes,
                        node_list=nodelist,
                        walltime=walltime,
                        init_blocks=1,
                        job_id=job_id,
                        max_blocks = 1,
                        worker_init_env=f"""~/spack/opt/spack/linux-almalinux8-thunderx2/gcc-8.5.0/python-venv-1.0-kpcvdxxd3vvp25x32hkihc2yhnw2yaso""",
                        launcher=SimplePMIxLauncher(),
                    ),
                )
                ],
                strategy="pmix_scale_simple",
                strategy_policy_file="/home/rbhattara/europar25/tests/policy.json",
            )

    start = time.time()

    parsl.load(config)

    # Call the first stage where we split the fastq file
    splitfiles=['sequence_collagen.fasta.split.{0:02d}'.format(i) for i in range(0, n_chunks)]
    fast_q_split(abspath(fastq_file), n_chunks).result()
    # We pass each of the split fastq file to bwa, and these run in parallel
    stage_start = time.time()

    bam_files = []
    for i in range(n_chunks):
        in_file = f"{path}{splitfiles[i]}"
        out_file = f"{path}{splitfiles[i].replace('fasta.split', 'bam')}"
        bam_files.append(blast(inputs=[File(in_file)], outputs=[File(out_file)]))

    print("Waiting for results ....")
    for i in range(n_chunks):
        bam_files[i].result(timeout=240)
    
    end = time.time()
    delta_stage = end - start

    bam_outputs = [r.outputs[0].filepath for r in bam_files]

    # The results from the blast are in bam_files, and these are passed to the merge app
    status = merge(inputs=bam_outputs, outputs=[File(abspath('/tmp/merged.result'))])
    status.result()

    delta = time.time() - start

    f = open(f"/home/rbhattara/europar25/tests/test/blast/blast1_runtime.txt", "w+")
    f.write(f"Runtime:{delta} Stage Runtime:{delta_stage}\n")
    f.close()

    if status.done():
        print("Cleaning up DFK")
        try:
            parsl.dfk().cleanup()
            parsl.clear()
            sys.exit(0)
        except Exception as e:
            print(f"Workflow failed: {e}", file=sys.stderr)
            sys.exit(1)  # Return nonzero on failure
