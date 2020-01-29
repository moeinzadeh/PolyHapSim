import sys

hapcompass_out_path = sys.argv[1]
ploidy = int(sys.argv[2])
out_path = sys.argv[3]

hapcompass_haps = {}
with open(hapcompass_out_path) as hapcompass_out:
    for line in hapcompass_out:
        if line == '\n':
            continue
        if 'BLOCK' in line:
            block_start_snp = int(line.split('\t')[3]) - 1  # MHM format is 0-based, so -1
            hapcompass_haps[block_start_snp] = [[] for n in range(ploidy)]
        else:
            line = line.strip().split('\t')
            for hap_n, allele in enumerate(line[3:]):
                hapcompass_haps[block_start_snp][hap_n].append(allele)

with open(out_path, 'a') as out:
    for block_n, block_start in enumerate(hapcompass_haps):
        for hap_n, hap in enumerate(hapcompass_haps[block_start]):
            out.write('\t'.join([str(block_n), '_'.join(['hapcompass', 'block', str(block_n), 'hap', str(hap_n)]),
                                 str(block_start), ''.join(hap)]) + '\n')
