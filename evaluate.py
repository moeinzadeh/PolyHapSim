import sys
from itertools import permutations

ploidy = int(sys.argv[1])
true_haps_file = sys.argv[2]
methods_to_evaluate = sys.argv[3]




def check_true_map(true_map,scaf_len):
    for i in true_map.values():
        for j in i:
            if len(j) != scaf_len:
                raise Exception('Grand truth haps are not correct')


def read_grand_truth(true_haps_file):
    true_haps_file = open(true_haps_file)
    lines = true_haps_file.readlines()

    true_map = {}
    inx = 0
    _, scaf_name, scaf_len = lines[inx].split()
    true_map[scaf_name] = []
    for i in range(ploidy):
        inx += 1
        true_map[scaf_name] += [lines[inx].split()[3]]
    check_true_map(true_map, int(scaf_len))
    return true_map


def read_a_methods_result(methods_haps_file):
    fi = open(methods_haps_file)
    lines = fi.readlines()

    i = 0
    block = 0
    map_method = {}
    while i < len(lines):
        if lines[i][0]==">":
            _, scaf , _ = lines[i].split()
            map_method [scaf] = {}
            i += 1
        block_index, _, start, hap = lines[i].split()
        base_block_index = block_index
        map_method[scaf][block_index] = []
        while base_block_index == block_index:
            map_method[scaf][block_index] += [(start, hap)]
            i += 1
            if i >= len(lines):
                break
            block_index, _, start, hap = lines[i].split()
    return map_method

def report_m_mm(h1,h2_start, h2):
    #m = match
    #mm = mismatch
    m,mm,gap = 0,0,0

    for i in range(h2_start,h2_start + len(h2)):
        if h1[i] == '-':
            pass
        elif h2[i-h2_start] == '-':
            gap += 1
        elif h1[i]==h2[i-h2_start]:
            m += 1
        else:
            mm +=1
    return m,mm,gap

def assign_haps(true_haps,method_haps):

    def score(m,mm):
        return m-mm


    #print true_haps,method_haps
    perm = permutations(range(len(true_haps)), len(method_haps))
    main_assignment_match,main_assignment_mismatch,main_assignment_score,main_assignment_gap = -100000,-100000,-1000000,-1000000
    for i in list(perm):
        block_match,block_mismatch,block_gap = 0,0,0
        for c,j in enumerate(i):
            m,mm,gap = report_m_mm(true_haps[j],int(method_haps[c][0]),method_haps[c][1])
            #print m,mm
            block_match += m
            block_mismatch += mm
            block_gap += gap
        #print block_match,block_mismatch
        #print
        if score(block_match,block_mismatch) > main_assignment_score:
            #for c, j in enumerate(i):
                #print true_haps[j]
                #print int(method_haps[c][0]) * '-' + method_haps[c][1],'   ',block_match, block_mismatch
            main_assignment_match, main_assignment_mismatch, main_assignment_score, main_assignment_gap = block_match,block_mismatch,score (block_match,block_mismatch),block_gap
    #print '---'
    return main_assignment_match,main_assignment_mismatch,main_assignment_score, main_assignment_gap

def reconstruction_rate(tm,mm,method):
    #tm true map
    #mm method map

    for m in tm:
        for j in mm[m]:
            res = assign_haps(tm[m],mm[m][j])
            print '----> result', method, res[0],res[1],res[2],res[3]





true_map = read_grand_truth (true_haps_file)

for method in methods_to_evaluate.split(','):
    file_path = true_haps_file
    method_file = file_path.replace ("true_haps",method+"_output")
    method_map = read_a_methods_result(method_file)
    reconstruction_rate(true_map,method_map,method)

