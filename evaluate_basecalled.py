import sys
import nadavca
import os
import numpy as np

if len(sys.argv) < 6 or len(sys.argv) > 7:
    print('usage: {} reads_dir basecalls_dir repeats_dir alignments_dir repeat_maps_dir [log_dir]'.format(sys.argv[0]))
    sys.exit(0)

reads_dir = sys.argv[1]
basecalls_dir = sys.argv[2]
repeats_dir = sys.argv[3]
alignments_dir = sys.argv[4]
repeat_maps_dir = sys.argv[5]
log_dir = None
if len(sys.argv) == 7:
    log_dir = sys.argv[6]

total_gt_repeat_positions, total_computed_repeat_positions = 0, 0
total_intersection, total_union = 0, 0
total_length = 0

for read_filename in os.listdir(reads_dir):
    print("processing {}".format(read_filename))
    
    basename = os.path.splitext(read_filename)[0]
    basecall_filename = basename + '.fasta'
    other_filename = basename + '.txt'
    
    read = nadavca.read.Read.load_from_fast5(os.path.join(reads_dir, read_filename), 'Analyses/Basecall_1D_000')
    
    basecall_path = os.path.join(basecalls_dir, basecall_filename)
    nadavca_result = nadavca.align_signal(basecall_path, [read])[0]
    if nadavca_result is None:
        continue
    _, basecall_alignment = nadavca_result
    basecall = nadavca.genome.Genome.load_from_fasta(basecall_path)[0]
    start_in_signal = np.full(len(basecall.bases), -1, dtype=int)
    end_in_signal = np.full(len(basecall.bases), -1, dtype=int)
    
    for ref_pos, event_start, event_end in basecall_alignment:
        start_in_signal[ref_pos] = event_start
        end_in_signal[ref_pos] = event_end
    
    computed_repeat_map = np.zeros(len(read.raw_signal), dtype=bool)
    with open(os.path.join(repeats_dir, other_filename), 'r') as file:
        for line in file:
            contig, rep_start, rep_end = line.split()[:3]
            rep_start = int(rep_start)
            rep_end = int(rep_end)
            computed_repeat_map[start_in_signal[rep_start] : end_in_signal[rep_end-1]] = 1

    ground_truth_repeat_map = np.zeros(len(read.raw_signal), dtype=bool)
    with open(os.path.join(repeat_maps_dir, other_filename), 'r') as file:
        for i, line in enumerate(file):
            ground_truth_repeat_map[i] = bool(int(line))
        
    
    aligned_start, aligned_end = None, None
    with open(os.path.join(alignments_dir, other_filename), 'r') as file:
        file.readline()
        file.readline()
        for line in file:
            ref_position, event_start, event_end = map(int, line.split())
            if aligned_start is None or event_start < aligned_start:
                aligned_start = event_start
            if aligned_end is None or event_end > aligned_end:
                aligned_end = event_end
    
    print("aligned part is {} - {}".format(aligned_start, aligned_end))
    
   
    ground_truth = ground_truth_repeat_map[aligned_start : aligned_end]
    computed = computed_repeat_map[aligned_start : aligned_end]
    
    length = aligned_end - aligned_start
    
    gt_count = np.sum(ground_truth)
    computed_count = np.sum(computed)
    union = np.sum(np.logical_or(ground_truth, computed))
    intersection = np.sum(np.logical_and(ground_truth, computed))
    iou = intersection / union
    
    print("iou: {} [computed: {}, ground: {}, intersection: {}, union: {}]".format(iou, np.sum(computed), np.sum(ground_truth), intersection, union))
    
    sensitivity = intersection / gt_count
    specificity = (length - union) / (length - gt_count)
    print("sensitivity: {}".format(sensitivity))
    print("specificity: {}".format(specificity))
    total_intersection += intersection
    total_union += union
    total_computed_repeat_positions += computed_count
    total_gt_repeat_positions += gt_count
    total_length += length
    print("cumulative iou: {}, cumulative sensitivity: {}, cumulative specificity: {}".format(total_intersection / total_union,
                                                                                              total_intersection / total_gt_repeat_positions,
                                                                                              (total_length - total_union) / (total_length - total_gt_repeat_positions)))
    if log_dir is not None:
        with open(os.path.join(log_dir, other_filename), 'w') as file:
            file.write('{}\n'.format(''.join(map(str, map(int, ground_truth)))))
            file.write('{}\n'.format(''.join(map(str, map(int, computed)))))
