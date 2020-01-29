import sys

sdhap_out_path = sys.argv[1]
ploidy = int(sys.argv[2])
out_path = sys.argv[3]

sdhap_haps = {}
with open(sdhap_out_path) as sdhap_out:
    for line in sdhap_out:
        if line == '\n':
            continue
        elif 'Block' in line:
            block_start = True
        else:
            if block_start:
                block_start_snp = int(line.split('\t')[0]) - 1  # MHM format is 0-based, so -1
                sdhap_haps[block_start_snp] = [[] for n in range(ploidy)]
                block_start = False
            line = line.strip().split('\t')
            for hap_n, allele in enumerate(line[2:]):
                sdhap_haps[block_start_snp][hap_n].append(allele)

with open(out_path, 'a') as out:
    for block_n, block_start in enumerate(sdhap_haps):
        for hap_n, hap in enumerate(sdhap_haps[block_start]):
            out.write('\t'.join([str(block_n), '_'.join(['sdhap', 'block', str(block_n), 'hap', str(hap_n)]),
                                 str(block_start), ''.join(hap)]) + '\n')
