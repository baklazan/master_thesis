import sys
import os
import numpy as np
from nadavca.genome import Genome
from nadavca.read import Read

if len(sys.argv) != 6:
    print("usage: {} reference repeats_file reads alignments output".format(sys.argv[0]))
    sys.exit(0)
reference_filename = sys.argv[1]
repeats_filename = sys.argv[2]
reads_dir = sys.argv[3]
alignments_dir = sys.argv[4]
output_dir = sys.argv[5]

references = Genome.load_from_fasta(reference_filename)
ref_lengths = {r.description[1:]: len(r.bases) for r in references}

reference_maps = {contig_name : np.zeros(length, dtype=bool) for contig_name, length in ref_lengths.items()}

with open(repeats_filename, 'r') as file:
    for line in file:
        tokens = line.split()
        contig = tokens[0]
        start = int(tokens[1])
        #start -= 1
        end = int(tokens[2])
        reference_maps[contig][start:end] = True

for read_filename in os.listdir(reads_dir):
    read_path = os.path.join(reads_dir, read_filename)
    read = Read.load_from_fast5(read_path, 'Analyses/Basecall_1D_000')
    length = len(read.raw_signal)
    read_map = np.zeros(length, dtype=bool)
    
    alignment_filename = os.path.splitext(read_filename)[0] + '.txt'
    alignment_path = os.path.join(alignments_dir, alignment_filename)
    with open(alignment_path, 'r') as file:
        contig = file.readline().rstrip()
        reference_map = reference_maps[contig]
        file.readline()
        last_event_end = None
        for line in file:
            reference_position, event_start, event_end = map(int, line.split())
            value = reference_map[reference_position]
            for i in range(event_start, event_end):
                read_map[i] = value
            if last_event_end is not None and last_event_end < event_start and value:
                for i in range(last_event_end, event_start):
                    read_map[i] = True
            last_event_end = event_end
    
    output_filename = alignment_filename
    output_path = os.path.join(output_dir, output_filename)
    with open(output_path, 'w') as file:
        file.write('\n'.join(map(str, map(int,read_map))))
        file.write('\n')
