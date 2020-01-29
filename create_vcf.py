import sys

true_haps_matrix_path = sys.argv[1]
out_vcf_path = sys.argv[2]

with open(true_haps_matrix_path) as true_haps_matrix, open(out_vcf_path, 'w') as out_vcf:
    out_vcf.write('##fileformat=VCFv4.2\n')
    out_vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tL1P1\n')

    true_haps_matrix.readline()
    for line in true_haps_matrix:
        line = line.split('\t')
        cigar = '1X'
        alt_alleles = line[4]
        n_alt_alleles = len(alt_alleles.split(','))
        if n_alt_alleles > 1:
            cigar = ','.join([cigar]*n_alt_alleles)
        out_vcf.write('\t'.join([line[1], line[2], '.', line[3], alt_alleles, '.', '.', 'AB=0;ABP=0;AC=;AF=0;AN=0;AO=0;CIGAR=' + cigar + ';DP=0;DPB=0;', 'GT', '/'.join(line[5:])]))