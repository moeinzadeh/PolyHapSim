import sys

true_haps_matrix_path = sys.argv[1]
ploidy = int(sys.argv[2])
ref_name = sys.argv[3]
out_path = sys.argv[4]

true_haps = [[] for n in range(ploidy)]
with open(true_haps_matrix_path) as true_haps_matrix:
    true_haps_matrix.readline()
    for line in true_haps_matrix:
        line = line.strip().split('\t')
        for hap_n, allele in enumerate(line[5:]):
            true_haps[hap_n].append(allele)

with open(out_path, 'w') as out:
    out.write('\t'.join(['>>>', ref_name, str(len(true_haps[0]))]) + '\n')
    for hap_n, hap in enumerate(true_haps):
        out.write('\t'.join(['0', 'true_hap_' + str(hap_n), '0', ''.join(hap)]) + '\n')
