import sys

hpop_out_path = sys.argv[1]
ploidy = int(sys.argv[2])
out_path = sys.argv[3]

hpop_haps = {}
with open(hpop_out_path) as hpop_out:
    for line in hpop_out:
        if '********' in line:
            continue
        elif 'BLOCK' in line:
            block_start_snp = int(line.split(' ')[2]) - 1  # MHM format is 0-based, so -1
            hpop_haps[block_start_snp] = [[] for n in range(ploidy)]
        else:
            line = line.strip().split('\t')
            for hap_n, allele in enumerate(line[1:]):
                hpop_haps[block_start_snp][hap_n].append(allele)

with open(out_path, 'a') as out:
    for block_n, block_start in enumerate(hpop_haps):
        for hap_n, hap in enumerate(hpop_haps[block_start]):
            out.write('\t'.join([str(block_n), '_'.join(['hpop', 'block', str(block_n), 'hap', str(hap_n)]),
                                 str(block_start), ''.join(hap)]) + '\n')
