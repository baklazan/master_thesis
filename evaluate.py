import sys
import os
from nadavca.read import Read
from scipy.special import softmax
import numpy as np
import signal_dtw

def log_stats():
    with open('stats.txt', 'w') as file:
        file.write('{} {} {} {} {}\n'.format(total_gt_repeat_positions, 
                                             total_computed_repeat_positions,
                                             total_intersection,
                                             total_union,
                                             total_length))

def log_read(basename):
    with open('done.txt', 'a') as file:
        file.write('{}\n'.format(basename))

if len(sys.argv) < 5 or len(sys.argv) > 6:
    print("usage: {} reads_dir chiron_predictions_dir alignments_dir repeat_maps_dir [logdir]".format(sys.argv[0]))
    sys.exit(0)


done_reads = set()

if os.path.exists('done.txt'):
    with open('done.txt', 'r') as file:
        for line in file:
            done_reads.add(line.rstrip())

reads_dir = sys.argv[1]
chiron_dir = sys.argv[2]
alignments_dir = sys.argv[3]
repeat_maps_dir = sys.argv[4]
logdir = None
if len(sys.argv) == 6:
    logdir = sys.argv[5]


total_gt_repeat_positions, total_computed_repeat_positions = 0, 0
total_intersection, total_union = 0, 0
total_length = 0

if os.path.exists('stats.txt'):
    with open('stats.txt', 'r') as file:
        tokens = file.readline().split()
        total_gt_repeat_positions, total_computed_repeat_positions, total_intersection, total_union, total_length = map(int, tokens)

log_stats()

for read_filename in os.listdir(reads_dir):
    basename = os.path.splitext(read_filename)[0]
    if basename in done_reads:
        continue
    other_filename = basename + '.txt'
    chiron_filename = basename + '.signal'
    
    print("processing {}".format(basename))
    
    read = Read.load_from_fast5(os.path.join(reads_dir, read_filename), 'Analyses/Basecall_1D_000')
    Read.normalize_reads([read])
    
    chiron_logits = []
    with open(os.path.join(chiron_dir, chiron_filename), 'r') as file:
        for line in file:
            chiron_logits.append(list(map(float, line.split())))
    chiron_logits = np.array(chiron_logits)
    chiron_probs = softmax(chiron_logits, axis=1)
    
    
    min_dist, max_dist = 3, 60
    bonus_for_moving = 0.13
    alignments = signal_dtw.local_alignment(read.normalized_signal, chiron_probs, min_dist, max_dist, bonus_for_moving, 4, 0.5, 250)
    
    computed_repeat_map = np.zeros(len(read.normalized_signal), dtype=bool)
    for path in alignments:
        computed_repeat_map[path[0][0] : path[-1][1]] = True
    
    ground_truth_repeat_map = np.zeros(len(read.normalized_signal), dtype=bool)
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
    if logdir is not None:
        with open(os.path.join(logdir, other_filename), 'w') as file:
            file.write('{}\n'.format(''.join(map(str, map(int, ground_truth)))))
            file.write('{}\n'.format(''.join(map(str, map(int, computed)))))
    log_stats()
    log_read(basename)

os.remove('done.txt')
os.remove('stats.txt')
