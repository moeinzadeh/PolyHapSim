class InputArguments:
    def __init__(self, **kwargs):
        """
        :param kwargs: dictionary with command line arguments
        """
        self.args = kwargs
        print ['='.join([str(arg), str(self.args[arg])]) for arg in self.args]

    def __getattr__(self, item):
        return self.args.get(item)


class Polymorphism:
    def __init__(self, vcf_line):
        """
        :param vcf_line: line from VCF file containing polymorphic site information
        """
        vcf_line = vcf_line.split('\t')
        self.attributes = {'ref_allele': vcf_line[3],
                           'alt_alleles': vcf_line[4].split(','),
                           'start_pos': int(vcf_line[1]) - 1}  # -1 for compatibility with pysam
        self.attributes['end_pos'] = self.attributes['start_pos'] + len(self.ref_allele) - 1
        self.attributes['span'] = range(self.attributes['start_pos'], self.attributes['end_pos'])
        if len(self.attributes['ref_allele']) == 1:
            self.attributes['type'] = 'SNP'
        elif set([len(allele) for allele in self.attributes['alt_alleles']]) == {len(self.attributes['ref_allele'])}:
            self.attributes['type'] = 'MNP'
        else:
            self.attributes['type'] = 'INDEL'

    def __getattr__(self, attribute):
        return self.attributes[attribute]

    def get_alleles(self):
        return [self.attributes['ref_allele']] + self.attributes['alt_alleles']

    def get_allele_index(self, allele):
        alleles = [self.attributes['ref_allele']] + self.attributes['alt_alleles']
        return alleles.index(allele)

    def is_multiallelic(self):
        return True if len([self.attributes['ref_allele']] + self.attributes['alt_alleles']) > 2 else False


class Polymorphisms:
    def __init__(self):
        self.polymorphisms = {}

    def __getitem__(self, pos):
        return self.polymorphisms[pos]

    def __len__(self):
        return len(self.polymorphisms)

    def keys(self):
        return self.polymorphisms.keys()

    def add_polymorphism(self, polymorphism):
        self.polymorphisms[polymorphism.start_pos] = polymorphism

    def get_polymorphisms_number(self, start, end):
        return len([pos for pos in self.polymorphisms.keys() if pos in range(start, end + 1)])


class Fragments:
    def __init__(self, vcf_path, solution_path=None, input_args=None, vcf_position=None, region=None, real_haps=False):
        """
        :param vcf_path: path to VCF file which was used for SNP matrix creation
        :param solution_path: path to phased haplotypes
        :param input_args: InputArgument class instance
        :param vcf_position: INT, region start position in vcf file
        :param region: STR, genomic region
        :param real_haps: bool, if input vcf file should be used as a template to store real haplotypes
        """
        if real_haps:
            self.fragments = {}
            self.fragments_map = {}

            with open(vcf_path, 'r') as vcf:
                vcf.seek(vcf_position)

                haps = [{} for n in range(input_args.ploidy)]
                for line in vcf:
                    if line.split('\t')[0] != region:
                        break

                    polymorphism = Polymorphism(vcf_line=line)

                    genotype = [int(allele) for allele in line.split('\t')[9].split(':')[0].split('/')]

                    for hap_num, hap in enumerate(haps):
                        hap[polymorphism.start_pos] = polymorphism.get_alleles()[genotype[hap_num]]
                    self.fragments_map[polymorphism.start_pos] = set(['hap_' + str(n)
                                                                      for n in range(input_args.ploidy)])

            for hap_num, hap in enumerate(haps):
                hap_name = '_'.join([region, 'hap', str(hap_num)])
                self.fragments[hap_name] = HapFragment(hap, name=hap_name)

        else:
            if input_args.haplotyper in ['SDhaP', 'HapCompass']:
                import subprocess

                bases_encoding = {'1': 'A', '2': 'C', '3': 'G', '4': 'T'}
                with open(solution_path, 'r') as solution, \
                     open(solution_path + '.mod', 'w') as solution_mod:
                    line = solution.readline()
                    solution_mod.write(line)
                    for line in solution:
                        if not line:
                            continue
                        if line[0] == 'B':
                            solution_mod.write('********\n')
                            solution_mod.write(line)
                            continue
                        if input_args.haplotyper == 'SDhaP':
                            solution_mod.write(line)
                        if input_args.haplotyper == 'HapCompass':
                            solution_mod.write('\t'.join(line.split('\t')[2:]))
                    solution_mod.write('********\n')
                solution_path += '.mod'

            index_to_snp = {}
            with open(vcf_path, 'r') as vcf:
                for line in vcf:
                    if line[:2] != '##':
                        break
                start = 1 if input_args.haplotyper in ['H-PoPG', 'HapCompass'] else 0
                for line_num, line in enumerate(vcf, start=start):
                    if input_args.haplotyper == 'HapCompass':
                        snp_pos = int(line.split('\t')[1]) - 1
                        #snp_pos = int(line.split('\t')[2][8:]) - 1 #Hossein change
                    else:
                        snp_pos = int(line.split('\t')[1]) - 1
                    index_to_snp[line_num] = snp_pos
            self.fragments = {}
            self.fragments_map = {snp_pos: set() for snp_pos in index_to_snp.values()}
            with open(solution_path, 'r') as solution:
                print solution_path
                for line in solution:
                    if 'BLOCK' in line:
                        block = [{} for n in range(input_args.ploidy)]
                    elif line[0] == '*' or line == "":
                        for hap in block:
                            if len(hap) > 1:
                                fragment = HapFragment(hap)
                                fragment_name = fragment.get_name()
                                self.fragments[fragment_name] = fragment
                                for snp_pos in fragment.get_positions():
                                    self.fragments_map[snp_pos].add(fragment_name)
                    else:
                        line = line.strip().split('\t')
                        print line
                        snp_pos = index_to_snp[int(line[0])]
                        for hap_num, allele in enumerate(line[1:]):
                            if allele != '-':
                                if input_args.haplotyper == 'SDhaP':
                                    block[hap_num][snp_pos] = bases_encoding[allele]
                                else:
                                    block[hap_num][snp_pos] = int(allele)

            if input_args.haplotyper in ['SDhaP', 'HapCompass']:
                subprocess.check_call(['rm', solution_path])

    def fetch_fragments(self, start, end):
        fragment_names = set()
        for snp_pos in set(range(start, end + 1)) & set(self.fragments_map.keys()):
            fragment_names.update(self.fragments_map[snp_pos])
        return [self.fragments[fragment_name] for fragment_name in list(fragment_names)]

    def fetch_all_fragments(self):
        return self.fragments.values()


class HapFragment:
    def __init__(self, d, name=None):
        """
        :param d: dictionary, {snp_pos: allele}
        :param name: STR, set haplotype name
        """
        self.fragment = d
        self.name = hash(frozenset(self.fragment.items())) if name is None else name

    def __getitem__(self, snp_pos):
        return self.fragment[snp_pos]

    def __len__(self):
        return len(self.fragment)

    def get_name(self):
        return self.name

    def get_positions(self):
        return self.fragment.keys()

    def get_leftmost_position(self):
        return min(self.fragment.keys())

    def get_rightmost_position(self):
        return max(self.fragment.keys())

    def get_allele(self, position):
        return self.fragment[position]

    def __str__(self):
        return str(self.fragment)


class Fragment:
    def __init__(self, alignment, polymorphisms, store_qualities=False, store_non_informative_alleles=False):
        """
        :param alignment: pysam AlignedSegment object
        :param polymorphisms: Polymorphisms class instance
        :param store_qualities: bool
        :param store_non_informative_alleles: bool
        """

        alignment_positions = alignment.get_reference_positions()
        if not alignment_positions:
            self.fragment = None
            return

        if not store_non_informative_alleles:
            alignment_snp_positions = list(set(alignment_positions) & set(polymorphisms.keys()))
        else:
            alignment_snp_positions = list(set(range(alignment_positions[0], alignment_positions[-1] + 1))
                                           & set(polymorphisms.keys()))

        # if there is less then 2 polymorphic sites in fragment, set fragment to None
        if len(alignment_snp_positions) < 2:
            self.fragment = None
            return

        # create the read map with {reference_pos: abs_read_pos}
        alignment_map = {ref_pos: read_pos
                         for read_pos, ref_pos in alignment.get_aligned_pairs()
                         if ref_pos is not None and read_pos is not None}

        # fill the fragment dictionary with alleles: {polymorphic_site: allele}
        query_sequence = alignment.query_sequence
        query_qualities = alignment.query_qualities
        self.fragment = {}
        self.qualities = {}

        for snp_pos in alignment_snp_positions:
            if polymorphisms[snp_pos].type == 'SNP':
                if snp_pos in alignment_map:
                    try:
                        allele = query_sequence[alignment_map[snp_pos]]
                        qual = query_qualities[alignment_map[snp_pos]]
                    except IndexError:
                        continue
                else:
                    allele = None
            else:
                if snp_pos in alignment_map and polymorphisms[snp_pos].end_pos in alignment_map:
                    if polymorphisms[snp_pos].end_pos + 1 not in alignment_map:
                        allele = query_sequence[alignment_map[snp_pos]:]
                        quals = query_qualities[alignment_map[snp_pos]:]
                    else:
                        allele = query_sequence[alignment_map[snp_pos]:alignment_map[polymorphisms[snp_pos].end_pos + 1]]
                        quals = query_qualities[alignment_map[snp_pos]:alignment_map[polymorphisms[snp_pos].end_pos + 1]]
                    qual = int(round(sum(quals) / float(len(quals))))
                else:
                    allele = None

            if allele in polymorphisms[snp_pos].get_alleles():
                self.fragment[snp_pos] = allele
                if store_qualities:
                    self.qualities[snp_pos] = chr(qual + 33)
            elif store_non_informative_alleles:
                self.fragment[snp_pos] = None

        self.name = alignment.query_name

        # If there are less then 2 informative alleles, set fragment to None
        if store_non_informative_alleles:
            if len([allele for allele in self.fragment.values() if allele is not None]) < 2:
                self.fragment = None
                return
        else:
            if len(self.fragment) < 2:
                self.fragment = None
                return

        self.leftmost_position = min(self.fragment.keys())
        self.rightmost_position = max(self.fragment.keys())

    def get_allele(self, position):
        return self.fragment[position] if self.fragment is not None else None

    def get_positions(self):
        return sorted(self.fragment.keys()) if self.fragment is not None else None

    def get_leftmost_position(self):
        return self.leftmost_position if self.fragment is not None else None

    def get_rightmost_position(self):
        return self.rightmost_position if self.fragment is not None else None

    def is_informative(self):
        return True if self.fragment is not None else False

    def get_name(self):
        return self.name

    def get_fragment(self):
        return self.fragment

    def get_qualities(self):
        return self.qualities


def evaluate_haplotypes(**kwargs):
    """Main haplotypes evaluation script"""
    import pysam
    import sympy
    import time
    import sys
    from datetime import timedelta
    import subprocess

    input_args = InputArguments(**kwargs)

    eval_bam_path = os.path.abspath(kwargs.get('eval_bam'))
    with pysam.AlignmentFile(eval_bam_path) as eval_bam:
        try:
            eval_bam.check_index()
        except ValueError:
            print '\033[93m' + eval_bam_path + '.bai not found, indexing bam file...\033[0m'
            subprocess.check_call(['samtools', 'index', eval_bam_path])

    vcf_path = os.path.abspath(kwargs.get('vcf'))
    vcf_index = index_vcf_file(vcf_path=vcf_path)

    haps_file_path = os.path.abspath(kwargs.get('haps_file'))
    if input_args.haplotyper == 'RANBOW':
        with pysam.AlignmentFile(haps_file_path) as haps_bam:
            try:
                haps_bam.check_index()
            except ValueError:
                print '\033[93m' + haps_file_path + '.bai not found, indexing bam file...\033[0m'
                subprocess.check_call(['samtools', 'index', haps_file_path])

    outfile_path = os.path.abspath(kwargs.get('outfile'))
    if not os.path.isdir(os.path.dirname(outfile_path)):
        os.makedirs(os.path.dirname(outfile_path))
    outfile = open(outfile_path, input_args.append)

    scoring_formula = sympy.sympify(kwargs.get('scoring_formula').lower())

    best_score_function = eval(kwargs.get('best_score_function').lower())

    # load regions for evaluation; if not specified, take all regions from evaluation bam file
    if input_args.region is None:
        if input_args.regions_file is None:
            with pysam.AlignmentFile(eval_bam_path) as eval_bam:
                regions = eval_bam.references
        else:
            with open(input_args.regions_file, 'r') as regions_file:
                regions = [r for r in regions_file.read().split()]
    else:
        regions = input_args.region

    start_time = time.time()
    for region_num, region in enumerate(regions):

        # update progress line
        elapsed_time = time.time() - start_time
        sys.stdout.write('\r')
        sys.stdout.write('%d/%d regions processed, elapsed time: %s' %
                         (region_num, len(regions), str(timedelta(seconds=elapsed_time))))
        sys.stdout.flush()

        # Load polymorphic sites in given region
        if input_args.haplotyper == 'RANBOW':
            region_polymorphisms = get_region_snp_positions(region=region, vcf_path=vcf_path,
                                                            vcf_position=vcf_index[region],
                                                            bam_path=haps_file_path, input_args=input_args)
        else:
            region_polymorphisms = get_region_snp_positions(region=region, vcf_path=vcf_path,
                                                            vcf_position=vcf_index[region], input_args=input_args)

        # Load phased solution if not RANBOW
        if input_args.haplotyper in ['SDhaP', 'H-PoPG', 'HapCompass']:
            haplotypes = Fragments(vcf_path=vcf_path, solution_path=haps_file_path, input_args=input_args)

        # Fetch all evaluation alignments from given region and store them as fragments
        with pysam.AlignmentFile(eval_bam_path) as eval_bam:
            alignments = eval_bam.fetch(region=region)
            eval_fragments = []

            for alignment in alignments:
                # Check if alignment passes MAPQ threshold and has at least 2 polymorphic sites
                if alignment.mapping_quality < input_args.mapq_threshold:
                    continue

                fragment = Fragment(alignment=alignment, polymorphisms=region_polymorphisms)
                eval_fragments.append(fragment)

        for eval_fragment in eval_fragments:
            # Discard non informative fragments (<2 polymorphic sites or <2 informative alleles)
            if not eval_fragment.is_informative():
                continue

            # Fetch haplotypes that overlap with region from first to last evaluation alignment SNP
            # and store them as a fragments
            if input_args.haplotyper == 'RANBOW':
                with pysam.AlignmentFile(haps_file_path) as haps_bam:
                    haplotypes = haps_bam.fetch(region=region,
                                                start=eval_fragment.get_leftmost_position(),
                                                end=eval_fragment.get_rightmost_position())
                    haplotype_fragments = []

                    for haplotype in haplotypes:
                        fragment = Fragment(alignment=haplotype, polymorphisms=region_polymorphisms)
                        haplotype_fragments.append(fragment)

            if input_args.haplotyper in ['SDhaP', 'H-PoPG', 'HapCompass']:
                haplotype_fragments = haplotypes.fetch_fragments(start=eval_fragment.get_leftmost_position(),
                                                                 end=eval_fragment.get_rightmost_position())

            # find the best-scoring haplotype and store #MATCH and #MISMATCH for it
            best_hap_name, matches, mismatches = find_best_score(evaluation_fragment=eval_fragment,
                                                                 haplotype_fragments=haplotype_fragments,
                                                                 scoring_formula=scoring_formula,
                                                                 best_score_function=best_score_function,
                                                                 input_args=input_args,
                                                                 polymorphic_sites=region_polymorphisms)
            # if not None, write evaluation read name, hap_name, #MATCH and #MISMATCH out
            if matches is None or mismatches is None:
                continue
            outfile.write(eval_fragment.get_name() + '\t'
                          + str(best_hap_name) + '\t'
                          + str(matches) + '\t'
                          + str(mismatches) + '\n')

    # update progress line
    elapsed_time = time.time() - start_time
    sys.stdout.write('\r')
    sys.stdout.write('%d/%d regions processed, elapsed time: %s\n' %
                     (len(regions), len(regions), str(timedelta(seconds=elapsed_time))))
    sys.stdout.flush()

    outfile.close()

    if input_args.plot_stats:
        subprocess.check_call(['Rscript',
                               os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'plot_hap_evaluation.R'),
                               outfile_path,
                               haps_file_path.replace('/', '_'),
                               kwargs.get('scoring_formula').upper(),
                               kwargs.get('best_score_function').upper(),
                               os.path.join(os.path.dirname(outfile_path),
                                            os.path.basename(haps_file_path) + '_eval_plot.png')])


def create_snp_matrix(**kwargs):
    """Main script to create input SNP matrices for given haplotyper"""
    import os
    import subprocess
    import pysam
    import time
    import sys
    from datetime import timedelta

    input_args = InputArguments(**kwargs)

    bam_path = os.path.abspath(kwargs.get('bam'))
    with pysam.AlignmentFile(bam_path) as bam_file:
        try:
            bam_file.check_index()
        except ValueError:
            print '\033[93m' + bam_path + '.bai not found, indexing bam file...\033[0m'
            subprocess.check_call(['samtools', 'index', bam_path])
    bam_file = pysam.AlignmentFile(bam_path)

    vcf_path = os.path.abspath(kwargs.get('vcf'))
    vcf_index = index_vcf_file(vcf_path=vcf_path)

    output_dir = os.path.abspath(kwargs.get('output_dir'))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    output_files_prefix = kwargs.get('output_files_prefix')
    if output_files_prefix is None:
        output_files_prefix = input_args.haplotyper

    if input_args.region is None:
        if input_args.regions_file is None:
            with pysam.AlignmentFile(bam_path) as bam:
                regions = bam.references
        else:
            with open(input_args.regions_file, 'r') as regions_file:
                regions = [r for r in regions_file.read().split()]
    else:
        regions = input_args.region

    start_time = time.time()
    for region_num, region in enumerate(regions):
        # update progress line
        elapsed_time = time.time() - start_time
        sys.stdout.write('\r')
        sys.stdout.write('%d/%d regions processed, elapsed time: %s' %
                         (region_num, len(regions), str(timedelta(seconds=elapsed_time))))
        sys.stdout.flush()

        region_vcf_file_path = os.path.join(output_dir, '_'.join([output_files_prefix, region]) + '.vcf')

        region_polymorphisms = get_region_snp_positions(region=region, vcf_path=vcf_path, vcf_position=vcf_index[region],
                                                        region_vcf_file_path=region_vcf_file_path, input_args=input_args)

        if input_args.create_vcf_only:
            continue

        # fetch all alignments in given region
        alignments = bam_file.fetch(region=region,
                                    start=sorted(region_polymorphisms.keys())[0],
                                    end=sorted(region_polymorphisms.keys())[-1])

        # run matrix creation script for given haplotyper
        if input_args.haplotyper == 'HapTree':
            output_file_path = os.path.join(output_dir, '_'.join([output_files_prefix, region]))
            create_haptree_matrices(alignments=alignments, polymorphic_sites=region_polymorphisms,
                                    output_file_path=output_file_path, input_args=input_args)

        else:
            output_matrix_path = os.path.join(output_dir, '_'.join([output_files_prefix, input_args.haplotyper.replace('-', '').lower(), 'snp_matrix']) + '.tsv')
            create_hpopg_or_sdhap_matrices(alignments=alignments, polymorphisms=region_polymorphisms,
                                           output_matrix_path=output_matrix_path, input_args=input_args)

    # update progress line
    elapsed_time = time.time() - start_time
    sys.stdout.write('\r')
    sys.stdout.write('%d/%d regions processed, elapsed time: %s\n' %
                     (len(regions), len(regions), str(timedelta(seconds=elapsed_time))))
    sys.stdout.flush()

    bam_file.close()


def evaluate_haplotypes_sim(**kwargs):
    """
    Main script to simulated data evaluation
    Output table:
    1 column: assembled haplotype name
    2 column: real haplotype name, having best score with assembled haplotype
    3 column: # match
    4 column: # mismatch
    5 column: # multiallelic sites match
    6 column: # multiallelic sites mismatch
    7 column: total polymorphic sites in assembled haplotype range
    8 column: position of first polymorphic site in assembled haplotype
    9 column: position of last polymorphic site in assembled haplotype
    """
    import pysam
    import sympy
    import time
    import sys
    from datetime import timedelta
    import subprocess

    input_args = InputArguments(**kwargs)

    origin_vcf_path = os.path.abspath(kwargs.get('origin_vcf'))
    hapl_vcf_path = os.path.abspath(kwargs.get('vcf'))
    vcf_index = index_vcf_file(vcf_path=origin_vcf_path)

    haps_file_path = os.path.abspath(kwargs.get('haps_file'))
    if input_args.haplotyper == 'RANBOW':
        with pysam.AlignmentFile(haps_file_path) as haps_bam:
            try:
                haps_bam.check_index()
            except ValueError:
                print '\033[93m' + haps_file_path + '.bai not found, indexing bam file...\033[0m'
                subprocess.check_call(['samtools', 'index', haps_file_path])

    outfile_path = os.path.abspath(kwargs.get('outfile'))
    if not os.path.isdir(os.path.dirname(outfile_path)):
        os.makedirs(os.path.dirname(outfile_path))
    outfile = open(outfile_path, input_args.append)

    scoring_formula = sympy.sympify(kwargs.get('scoring_formula').lower())

    best_score_function = eval(kwargs.get('best_score_function').lower())

    # load regions for evaluation
    if input_args.region is None:
        if input_args.regions_file is None:
            sys.exit('Need --region or --regions_file argument.')
        else:
            with open(input_args.regions_file, 'r') as regions_file:
                regions = [r for r in regions_file.read().split()]
    else:
        regions = input_args.region

    supplementary_outfile = open(os.path.join(os.path.dirname(outfile_path),
                                              'supplementary_' + os.path.basename(outfile_path)), 'a')

    start_time = time.time()
    for region_num, region in enumerate(regions):

        # update progress line
        elapsed_time = time.time() - start_time
        sys.stdout.write('\r')
        sys.stdout.write('%d/%d regions processed, elapsed time: %s' %
                         (region_num, len(regions), str(timedelta(seconds=elapsed_time))))
        sys.stdout.flush()

        # Load polymorphic sites in given region
        region_polymorphisms = get_region_snp_positions(region=region, vcf_path=origin_vcf_path,
                                                        vcf_position=vcf_index[region], input_args=input_args)

        # Write real haplotypes length to supplementary output file
        supplementary_outfile.write('\t'.join([region, str(len(region_polymorphisms))]) + '\n')

        # Load real haplotypes
        real_haplotypes = Fragments(vcf_path=origin_vcf_path, input_args=input_args,
                                    vcf_position=vcf_index[region], region=region, real_haps=True)

        # Load assembled haplotypes if not RANBOW
        if input_args.haplotyper in ['SDhaP', 'H-PoPG', 'HapCompass']:
            haplotypes = Fragments(vcf_path=hapl_vcf_path, solution_path=haps_file_path, input_args=input_args)
            haplotypes = haplotypes.fetch_all_fragments()

        # Fetch all assembled haplotypes from given region and store them as fragments if RANBOW
        else:
            with pysam.AlignmentFile(haps_file_path) as haps_bam:
                haplotypes = [Fragment(alignment=haplotype, polymorphisms=region_polymorphisms)
                              for haplotype in haps_bam.fetch(region=region)]

        real_haplotypes_list = real_haplotypes.fetch_all_fragments()
        for haplotype in haplotypes:
            # find the best-scoring haplotype and store #MATCH, #MISMATCH, #MULTALLELE_MATCH, #MULTALLELE_MISMATCH
            match_hap_name, matches, mismatches, multallele_matches, multallele_mismatches =\
                find_best_score(evaluation_fragment=haplotype, haplotype_fragments=real_haplotypes_list,
                                scoring_formula=scoring_formula, best_score_function=best_score_function,
                                input_args=input_args, polymorphic_sites=region_polymorphisms)
            if matches is None or mismatches is None:
                continue
            outfile.write('\t'.join([str(haplotype.get_name()), str(match_hap_name), str(matches), str(mismatches),
                                     str(multallele_matches), str(multallele_mismatches),
                                     str(0),
                                     str(0), str(0)])
                          + '\n')

    # update progress line
    elapsed_time = time.time() - start_time
    sys.stdout.write('\r')
    sys.stdout.write('%d/%d regions processed, elapsed time: %s\n' %
                     (len(regions), len(regions), str(timedelta(seconds=elapsed_time))))
    sys.stdout.flush()

    outfile.close()

    supplementary_outfile.close()

    if input_args.plot_stats:
        subprocess.check_call(['Rscript',
                               os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'plot_hap_evaluation.R'),
                               outfile_path,
                               haps_file_path.replace('/', '_'),
                               kwargs.get('scoring_formula').upper(),
                               kwargs.get('best_score_function').upper(),
                               os.path.join(os.path.dirname(outfile_path),
                                            os.path.basename(haps_file_path) + '_eval_plot.png')])


def get_region_snp_positions(region, vcf_path, vcf_position, input_args, bam_path=None, region_vcf_file_path=None):
    """
    :param region: STR, genomic region name
    :param vcf_path: STR, path to the vcf file
    :param vcf_position: INT, starting position of region in vcf file
    :param input_args: InputArguments class instance
    :param bam_path: STR, path to the bam file
    :param region_vcf_file_path: STR, path to the output vcf file with selected region SNP's
    :return: Polymorphisms class instance

    Return a Polymorphisms class instance containing all polymorphisms informations from given region,
    add bam_path to not include SNP positions with no alignments on it.
    If the region_vcf_name is definded, the selected vcf lines will be stored as region_vcf_file as well.
    """
    import pysam

    polymorphisms = Polymorphisms()

    if bam_path is not None:
        bam_file = pysam.AlignmentFile(bam_path)

    if region_vcf_file_path is not None:
        output_vcf_file = open(region_vcf_file_path, 'w')
        output_vcf_file.write('##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tunknown\n')

    with open(vcf_path, 'r') as vcf_file:
        vcf_file.seek(vcf_position)
        for line in vcf_file:
            if line.split('\t')[0] != region:
                break

            if input_args.utility_name == 'evaluate_haplotypes' and input_args.haplotyper == 'HapCompass':
                line = line.split('\t')
                line = '\t'.join([line[0], line[2][8:], line[1]] + line[3:])
                polymorphism = Polymorphism(vcf_line=line)
            else:
                polymorphism = Polymorphism(vcf_line=line)

            if input_args.haplotyper == 'SDhaP' and polymorphism.type != 'SNP':
                continue

            if bam_path is not None:
                # Check if the depth on polymorphic site > 0
                if not bam_file.count(region=region, start=polymorphism.start_pos, end=polymorphism.end_pos):
                    continue

            polymorphisms.add_polymorphism(polymorphism)

            if region_vcf_file_path is not None:
                output_vcf_file.write(line)

    if bam_path is not None:
        bam_file.close()

    if region_vcf_file_path is not None:
        output_vcf_file.close()
    return polymorphisms


def index_vcf_file(vcf_path):
    """
    Index a vcf file, write it to vcf_file.index
        region name         position in vcf file
    and return it also as a dictionary:
        {region name: position in vcf file}
    """

    index_path = vcf_path + '.index'
    vcf_index = {}
    with open(vcf_path, 'r') as vcf_file, open(index_path, 'w') as index_file:
        #  skip headers and store pointer position
        while True:
            line = vcf_file.readline()
            if line[:2] != '##':
                vcf_position = vcf_file.tell()
                line = vcf_file.readline()
                break

        while line:
            region_name = line.split('\t')[0]
            index_file.write(region_name + '\t' + str(vcf_position) + '\n')
            vcf_index[region_name] = vcf_position

            while line and region_name == line.split('\t')[0]:
                vcf_position = vcf_file.tell()
                line = vcf_file.readline()
    return vcf_index


def find_best_score(evaluation_fragment, haplotype_fragments, scoring_formula, best_score_function, input_args,
                    polymorphic_sites):
    """
    Find haplotype with highest score and return evaluation read name, haplotype name, #MATCH and #MISMATCH
    Score in computed using the scoring_formula and best score is choosen using the best_score_function argument
    Only SNP's are used to compute the score

    :param evaluation_fragment: Fragment object used to evaluation
    :param haplotype_fragments: list containing Fragment or HapFragment instances - haplotypes fragments
    :param scoring_formula: sympy object containing formula to compute the score
    :param best_score_function: function to choose between two scores; currently implemented: max, min
    :param input_args: InputArguments instance representing command line arguments
    :param polymorphic_sites: Polymorphisms class instance
    :return: STR, INT, INT, INT, INT
    """

    best_score = None
    best_score_matches = None
    best_score_mismatches = None
    best_score_haplotype = None

    for hap_fragment in haplotype_fragments:
        if input_args.haplotyper == 'RANBOW' and input_args.utility_name != 'evaluate_haplotypes_sim':
            # Discard non informative fragments (<2 polymorphic sites or <2 informative alleles)
            if not hap_fragment.is_informative():
                continue

        try:
            common_polymorphisms = list(set(hap_fragment.get_positions()) & set(evaluation_fragment.get_positions()))
        except TypeError:
            continue

        # Discard haplotypes if #common polymorphic sites with alignment < 2
        if len(common_polymorphisms) < 2:
            continue

        # Count #MATCH and #MISMATCH
        matches = 0
        mismatches = 0
        for snp in common_polymorphisms:
            if input_args.utility_name == 'evaluate_haplotypes_sim':
                alignment_allele = evaluation_fragment.get_allele(position=snp)
                if input_args.haplotyper in ['RANBOW', 'SDhaP']:
                    haplotype_allele = hap_fragment.get_allele(position=snp)
                else:
                    haplotype_allele = polymorphic_sites[snp].get_allele_index(hap_fragment.get_allele(snp))
            else:
                haplotype_allele = hap_fragment.get_allele(position=snp)
                if input_args.haplotyper in ['RANBOW', 'SDhaP']:
                    alignment_allele = evaluation_fragment.get_allele(position=snp)
                else:
                    alignment_allele = polymorphic_sites[snp].get_allele_index(evaluation_fragment.get_allele(snp))

            if haplotype_allele == alignment_allele:
                matches += 1
            else:
                mismatches += 1

        # Discard score if #MATCH + #MISMATCH < 2
        if matches + mismatches < 2:
            continue

        score = float(scoring_formula.evalf(subs={'m': matches, 'mm': mismatches}))
        best_score = best_score_function([s for s in [score, best_score] if s is not None])
        if score == best_score:
            best_haplotype = hap_fragment
            best_score_matches = matches
            best_score_mismatches = mismatches
            best_score_haplotype = hap_fragment.get_name()

    if input_args.utility_name == 'evaluate_haplotypes_simm':
        # Get number of matched and mismatched multiallelic polymorphisms
        common_polymorphisms = list(set(best_haplotype.get_positions()) & set(evaluation_fragment.get_positions()))
        best_score_ma_matches = 0
        best_score_ma_mismatches = 0
    
        if input_args.haplotyper in ['RANBOW', 'SDhaP']:
            for position in common_polymorphisms:
                if not polymorphic_sites[position].is_multiallelic():
                    continue

                alignment_allele = evaluation_fragment.get_allele(position=position)
                haplotype_allele = best_haplotype.get_allele(position=position)

                if haplotype_allele == alignment_allele:
                    best_score_ma_matches += 1
                else:
                    best_score_ma_mismatches += 1

        return best_score_haplotype, best_score_matches, best_score_mismatches,\
               best_score_ma_matches, best_score_ma_mismatches
    else:
        return best_score_haplotype, best_score_matches, best_score_mismatches, 0, 0


def create_haptree_matrices(alignments, polymorphic_sites, output_file_path, input_args):
    """
    :param alignments: an iterator over a collection of pysam ALignedSegment objects
    :param polymorphic_sites: Polymorphisms class instance
    :param output_file_path: STR, path to write the output file
    :param input_args: InputArguments instance representing command line arguments

    Produce SNP input file to use with HapTree v0.1
    alignment -> {SNP_number: 0/1}, 0 if alignment allele = reference allele, 1 otherwise.
    """
    # Convert alignments to fragments and store them
    fragments = []
    for alignment in alignments:
        if alignment.mapping_quality < input_args.mapq_threshold:
            continue
        if input_args.filter_supplementary and alignment.is_supplementary:
            continue
        if input_args.filter_secondary and alignment.is_secondary:
            continue
        fragments.append(Fragment(alignment=alignment, polymorphisms=polymorphic_sites))

    # if reads are pair-ended, combine the mates by the read name
    if input_args.pair_end:
        reads = {}
        for fragment in fragments:
            if not fragment.is_informative():
                continue

            if fragment.get_name() not in reads:
                reads[fragment.get_name()] = [fragment]
            else:
                reads[fragment.get_name()].append(fragment)

    # Index the snp positions
    snp_pos_map = {snp_pos: snp_index for snp_index, snp_pos in enumerate(sorted(polymorphic_sites.keys()))}

    # Create dictionaries to divide polymorphic sites and reads into connected blocks
    snp_blocks = {}  # {polymorphic_site: block_num}
    read_blocks = {}  # {block_number: Fragment}

    # Combine pair-ended reads
    if input_args.pair_end:
        fragments = []
        for read_name in reads:
            mates_dict = {}
            for fragment in reads[read_name]:
                fragment_dict = fragment.get_fragment()
                for snp_pos in fragment_dict:
                    # Conflict handling for overlapping mate reads
                    if snp_pos in mates_dict:
                        if fragment_dict[snp_pos] == mates_dict[snp_pos]:
                            continue
                        else:
                            mates_dict[snp_pos] = None
                    else:
                        mates_dict[snp_pos] = fragment_dict[snp_pos]
            fragments.append(HapFragment({snp_pos: allele for snp_pos, allele in mates_dict.items()
                                          if allele is not None}))

    # Single-end reads and combined pair-end reads
    for fragment in fragments:
        if not input_args.pair_end:
            if not fragment.is_informative():
                continue

        fragment_snps = fragment.get_positions()
        connected_blocks = []
        for snp_pos in fragment_snps:
            if snp_pos not in snp_blocks:
                continue
            else:
                connected_blocks.append(snp_blocks[snp_pos])

        if not connected_blocks:
            # create new block if no snp position overlaps with known blocks
            block_num = max(snp_blocks.values()) + 1 if snp_blocks else 0
            for snp_pos in fragment_snps:
                snp_blocks[snp_pos] = block_num
            read_blocks[block_num] = [fragment]
        else:
            block_num = min(connected_blocks)
            for snp_pos in fragment_snps:
                snp_blocks[snp_pos] = block_num
            read_blocks[block_num].append(fragment)
            if len(connected_blocks) > 1:
                blocks_to_delete = set(connected_blocks) - {block_num}
                for snp_pos in snp_blocks:
                    if snp_blocks[snp_pos] in blocks_to_delete:
                        snp_blocks[snp_pos] = block_num
                for block_to_delete in blocks_to_delete:
                    read_blocks[block_num] = read_blocks[block_num] + read_blocks[block_to_delete]
                    del read_blocks[block_to_delete]

    snp_blocks_map = {}
    for block_num, block in enumerate(read_blocks.keys()):
        with open(output_file_path + '.vcf') as vcf, \
             open(''.join([output_file_path, '_block', str(block_num), '.vcf']), 'w') as out_vcf:
            lines_to_write = [snp_pos_map[snp_pos] for snp_pos in snp_blocks if snp_blocks[snp_pos] == block]
            line_abs_num = 0

            for line in vcf:
                if line[:2] != '##':
                    break

            for line_num, line in enumerate(vcf):
                if line_num in lines_to_write:
                    out_vcf.write(line)
                    pos = int(line.split('\t')[1]) - 1
                    snp_blocks_map[pos] = line_abs_num
                    line_abs_num += 1
        with open(''.join([output_file_path, '_block', str(block_num), '.reads']), 'w') as output_file:
            for fragment in read_blocks[block]:
                fragment_str = str({snp_blocks_map[snp_pos]: int(fragment.get_allele(snp_pos)
                                                                 != polymorphic_sites[snp_pos].ref_allele)
                                    for snp_pos in fragment.get_positions()})
                output_file.write(fragment_str + '\n')


def create_hpopg_or_sdhap_matrices(alignments, polymorphisms, output_matrix_path, input_args):
    """
    :param alignments: an iterator over a collection of pysam ALignedSegment objects
    :param polymorphisms: Polymorphisms class instance
    :param output_matrix_path: STR, path to write the output SNP matrix
    :param input_args: InputArguments instance representing command line arguments

    Produce SNP input matrix to work with H-PoPG or with SDhaP polyploid haplotyper.
    H-PoPG:
    read blocks number \t read name \t first block start SNP index \t first block 0/1 sequence \t quality values
    If there are more than 1 blocks, additional columns 3 and 4 are needed

    SDhaP:
    Number of reads
    Number of SNP sites in region
    read blocks number \t read name \t first block start SNP index \t first block bases \t quality values
    If there are more than 1 blocks, additional columns 3 and 4 are needed
    """
    if input_args.haplotyper == 'SDhaP':
        from string import maketrans
        import subprocess
        bases_encoding = maketrans('ACGT', '1234')
        reads_counter = 0

    output_matrix = open(output_matrix_path, 'w')

    # Index the region snp positions
    sorted_snp_positions = sorted(polymorphisms.keys())
    snp_to_index = {snp_pos: snp_index for snp_index, snp_pos in enumerate(sorted_snp_positions, start=1)}

    # Convert alignments to fragments and store them
    if input_args.pair_end:
        reads = {}
        for alignment in alignments:
            if alignment.mapping_quality < input_args.mapq_threshold:
                continue
            if input_args.filter_supplementary and alignment.is_supplementary:
                continue
            if input_args.filter_secondary and alignment.is_secondary:
                continue

            read_name = alignment.query_name

            if read_name in reads:
                reads[read_name].append(Fragment(alignment=alignment, polymorphisms=polymorphisms,
                                                 store_qualities=True, store_non_informative_alleles=True))
            else:
                reads[read_name] = [Fragment(alignment=alignment, polymorphisms=polymorphisms,
                                             store_qualities=True, store_non_informative_alleles=True)]

        for read_name in reads:
            mates = reads[read_name]
            mates_dict = {}
            mates_quals = {}

            for fragment in mates:
                if not fragment.is_informative():
                    continue

                fragment_dict = fragment.get_fragment()
                fragment_quals = fragment.get_qualities()
                for snp_pos in fragment_dict:
                    # Conflict handling by mates, aligned to same place, or by overlapping mate reads
                    fragment_allele = fragment_dict[snp_pos]
                    if snp_pos in mates_dict:
                        if fragment_allele == mates_dict[snp_pos]:
                            try:
                                if fragment_quals[snp_pos] == mates_quals[snp_pos]:
                                    continue
                                else:
                                    mates_quals[snp_pos] = max(fragment_quals[snp_pos], mates_quals[snp_pos])
                                    continue
                            except KeyError:
                                continue
                        else:
                            mates_dict[snp_pos] = None

                    else:
                        mates_dict[snp_pos] = fragment_allele
                        if fragment_allele is not None:
                            mates_quals[snp_pos] = fragment_quals[snp_pos]

            if not mates_dict:
                continue

            if len([allele for allele in mates_dict.values() if allele is not None]) < 2:
                continue

            # add None alleles between the mates to indicate gap
            mates_leftmost_pos = snp_to_index[min(mates_dict.keys())]
            mates_rightmost_pos = snp_to_index[max(mates_dict.keys())]
            mates_snp_interval = sorted_snp_positions[mates_leftmost_pos : mates_rightmost_pos+1]

            for snp_pos in mates_snp_interval:
                if snp_pos not in mates_dict:
                    mates_dict[snp_pos] = None

            blocks = []
            block = []
            block_start_positions = []

            for snp_pos in sorted(mates_dict.keys()):
                if not block:
                    allele = mates_dict[snp_pos]
                    if allele is None:
                        continue
                    block_start_positions.append(snp_pos)
                    if input_args.haplotyper == 'H-PoPG':
                        block.append(str(polymorphisms[snp_pos].get_allele_index(allele)))
                    if input_args.haplotyper == 'SDhaP':
                        block.append(allele)
                else:
                    allele = mates_dict[snp_pos]
                    if allele is None:
                        blocks.append(''.join(block))
                        block = []
                        continue
                    if input_args.haplotyper == 'H-PoPG':
                        block.append(str(polymorphisms[snp_pos].get_allele_index(allele)))
                    if input_args.haplotyper == 'SDhaP':
                        block.append(allele)

            if block:
                blocks.append(''.join(block))

            if input_args.haplotyper == 'SDhaP':
                for block_index, seq in enumerate(blocks):
                    blocks[block_index] = seq.translate(bases_encoding)

            block_start_positions = [str(snp_to_index[pos]) for pos in block_start_positions]
            blocks_str = '\t'.join(['\t'.join([start_pos, seq])
                                    for (start_pos, seq) in zip(block_start_positions, blocks)])
            block_quals = ''.join([mates_quals[pos] for pos in sorted(mates_quals.keys())])
            output_matrix.write('\t'.join([str(len(blocks)), read_name, blocks_str, block_quals]) + '\n')

            if input_args.haplotyper == 'SDhaP':
                reads_counter += 1

    if not input_args.pair_end:
        fragments = []
        for alignment in alignments:
            if alignment.mapping_quality < input_args.mapq_threshold:
                continue
            if input_args.filter_supplementary and alignment.is_supplementary:
                continue
            if input_args.filter_secondary and alignment.is_secondary:
                continue

            fragments.append(Fragment(alignment=alignment, polymorphisms=polymorphisms,
                                      store_qualities=True, store_non_informative_alleles=True))

        for fragment in fragments:
            if not fragment.is_informative():
                continue

            fragment_dict = fragment.get_fragment()
            blocks = []
            block = []
            block_start_positions = []

            for snp_pos in sorted(fragment_dict.keys()):
                if not block:
                    allele = fragment_dict[snp_pos]
                    if allele is None:
                        continue
                    block_start_positions.append(snp_pos)
                    if input_args.haplotyper == 'H-PoPG':
                        block.append(str(polymorphisms[snp_pos].get_allele_index(allele)))
                    if input_args.haplotyper == 'SDhaP':
                        block.append(allele)
                else:
                    allele = fragment_dict[snp_pos]
                    if allele is None:
                        blocks.append(''.join(block))
                        block = []
                        continue
                    if input_args.haplotyper == 'H-PoPG':
                        block.append(str(polymorphisms[snp_pos].get_allele_index(allele)))
                    if input_args.haplotyper == 'SDhaP':
                        block.append(allele)

            if block:
                blocks.append(''.join(block))

            if input_args.haplotyper == 'SDhaP':
                for block_index, seq in enumerate(blocks):
                    blocks[block_index] = seq.translate(bases_encoding)

            block_start_positions = [str(snp_to_index[pos]) for pos in block_start_positions]
            blocks_str = '\t'.join(['\t'.join([start_pos, seq])
                                    for (start_pos, seq) in zip(block_start_positions, blocks)])
            block_quals = ''.join([fragment.get_qualities()[pos] for pos in sorted(fragment.get_qualities().keys())])
            read_name = fragment.get_name()
            output_matrix.write('\t'.join([str(len(blocks)), read_name, blocks_str, block_quals]) + '\n')

            if input_args.haplotyper == 'SDhaP':
                reads_counter += 1

    output_matrix.close()
    if input_args.haplotyper == 'SDhaP':
        temp_matrix_path = output_matrix_path + '.temp'
        with open(temp_matrix_path, 'w') as temp_matrix:
            temp_matrix.write(str(reads_counter) + '\n' + str(len(polymorphisms)) + '\n')
            with open(output_matrix_path, 'r') as output_matrix:
                for line in output_matrix:
                    temp_matrix.write(line)
        subprocess.check_call(['rm', output_matrix_path])
        subprocess.check_call(['mv', temp_matrix_path, output_matrix_path])

if __name__ == '__main__':
    import argparse
    import os

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='Utilities kit for haplotype assembly.')
    subparsers = parser.add_subparsers(dest='utility_name',
                                       help='python utils.py utility_name to see the utility arguments.')

    parser_eval = subparsers.add_parser('evaluate_haplotypes', formatter_class=argparse.RawTextHelpFormatter,
                                        help='Evaluate haplotypes, write statistics to file.')
    parser_eval.add_argument('haplotyper', choices=['RANBOW', 'H-PoPG', 'SDhaP', 'HapCompass'],
                             help='Output of which haplotyper should be evaluated.')
    parser_eval.add_argument('eval_bam', help='Path to alignment file with evaluation reads.')
    parser_eval.add_argument('haps_file', help='Path to bam file (text file) with haplotypes you want to evaluate.')
    parser_eval.add_argument('vcf', help='Path to vcf file.')
    parser_eval.add_argument('outfile', help='Path to output file.')
    parser_eval.add_argument('--append', action='store_const', const='a', default='w',
                             help='Set this flag to append statistics to outfile instead of creating a new one.')
    parser_eval.add_argument('--ploidy', metavar='<INT>', type=int,
                             help='Specify ploidy number of organism,'
                                  'if the output file produced with H-PoPG, SDhaP or HapCompass')
    parser_eval.add_argument('--region', '-r', metavar='<STR>', action='append', default=None,
                             help='Specify genomic region. Multiple --region arguments are possible:\n'
                                  '\t--region chr1 --region chr2')
    parser_eval.add_argument('--regions_file', metavar='<FILENAME>', default=None,
                             help='Path to text file with region names. Used only if --region argument not given. '
                                  'Use all regions from bam file header if neither --region '
                                  'nor --regions_file options are not specified.')
    parser_eval.add_argument('--mapq_threshold', metavar='<INT>', type=int, default=0,
                             help='Discard evaluation alignments with MAPQ < INT (default: %(default)s).')
    parser_eval.add_argument('--scoring_formula', metavar='<STR>', type=str,
                             help='Specify the formula for haplotype scoring, where\n'
                                  '\tm is #match\n\tmm is #mismatch\n(default: "%(default)s").', default='m - mm**2')
    parser_eval.add_argument('--best_score_function', type=str, choices=['max', 'min'], default='max',
                             help='Specify the function to choose the best score (default: %(default)s).')
    parser_eval.add_argument('--plot_stats', action='store_true',
                             help='Set this flag to plot #Match ~ #Mismatch distribution.'
                                  'The plot will be stored in the same directory with output statistics file.')

    parser_snp_matrix = subparsers.add_parser('create_snp_matrix', formatter_class=argparse.RawTextHelpFormatter,
                                              help='Create input SNP matrices and distinct vcf from bam and vcf file '
                                                   'for different haplotypers.')
    parser_snp_matrix.add_argument('haplotyper', choices=['SDhaP', 'HapTree', 'H-PoPG'],
                                   help='Set a haplotyper for which you want to create SNP matrix. '
                                        'Currently implemented:\n'
                                        '\tSDhaP (polyploid only): https://doi.org/10.1186/s12864-015-1408-5\n'
                                        '\tHapTree v0.1: http://cb.csail.mit.edu/cb/haptree\n'
                                        '\tH-PoPG: https://github.com/MinzhuXie/H-PoPG\n')
    parser_snp_matrix.add_argument('bam', type=str, help='Path to bam file with alignments.')
    parser_snp_matrix.add_argument('vcf', type=str, help='Path to vcf file.')
    parser_snp_matrix.add_argument('output_dir', type=str, help='Path to directory to store the output matrices.')
    parser_snp_matrix.add_argument('--output_files_prefix', '-prefix', metavar='<STR>', type=str,
                                   help='Prefix for output SNP matrices and vcf files: prefix_region.tsv/vcf\n'
                                        '(Default: haplotyper_region.tsv/vcf).')
    parser_snp_matrix.add_argument('--region', '-r', metavar='<STR>', action='append', default=None,
                                   help='Specify genomic region. Multiple --region arguments are possible:\n'
                                        '\t--region chr1 --region chr2')
    parser_snp_matrix.add_argument('--regions_file', '-rfile', metavar='<FILENAME>', type=str, default=None,
                                   help='Path to text file with region names. Used only if --region argument not given.'
                                        ' Use all regions from bam file header if neither --region'
                                        'nor --regions_file options are not specified.')
    parser_snp_matrix.add_argument('--pair_end', '-pe', action='store_true', default=False,
                                   help='The reads are pair-end (Default: consider all reads as single-end).')
    parser_snp_matrix.add_argument('--mapq_threshold', '-mapq', metavar='<INT>', type=int, default=0,
                                   help='Mapping quality threshold (default: %(default)s).')
    parser_snp_matrix.add_argument('--filter_supplementary', '-nosup', action='store_true',
                                   help='Set this flag to filter all supplementary alignments out.')
    parser_snp_matrix.add_argument('--filter_secondary', '-nosec', action='store_true',
                                   help='Set this flag to filter all secondary alignments out.')
    parser_snp_matrix.add_argument('--create_vcf_only', '-vcfonly', action='store_true',
                                   help='Set this flag to create vcf files only')

    parser_eval_sim = subparsers.add_parser('evaluate_haplotypes_sim', formatter_class=argparse.RawTextHelpFormatter,
                                            help='Evaluate haplotypes with known origin haplotypes'
                                                 '(e.g. simulated data).')
    parser_eval_sim.add_argument('haplotyper', choices=['RANBOW', 'H-PoPG', 'SDhaP', 'HapCompass'],
                                 help='Output of which haplotyper should be evaluated.')
    parser_eval_sim.add_argument('haps_file', help='Path to bam file (text file) with haplotypes you want to evaluate.')
    parser_eval_sim.add_argument('vcf', help='Path to vcf file which were specifically created for haplotyping.')
    parser_eval_sim.add_argument('origin_vcf', help='Path to vcf file with GT field.')
    parser_eval_sim.add_argument('outfile', help='Path to output file.')
    parser_eval_sim.add_argument('ploidy', type=int, help='Specify ploidy number of organism.')
    parser_eval_sim.add_argument('--append', action='store_const', const='a', default='w',
                                 help='Set this flag to append statistics to outfile instead of creating a new one.')
    parser_eval_sim.add_argument('--region', '-r', metavar='<STR>', action='append', default=None,
                                 help='Specify genomic region. Multiple --region arguments are possible:\n'
                                      '\t--region chr1 --region chr2')
    parser_eval_sim.add_argument('--regions_file', metavar='<FILENAME>', default=None,
                                 help='Path to text file with region names. Used only if --region argument not given. '
                                      '--region or --regions_file is required')
    parser_eval_sim.add_argument('--scoring_formula', metavar='<STR>', type=str,
                                 help='Specify the formula for haplotype scoring, where\n\t'
                                      'm is #match\n\tmm is #mismatch\n(default: "%(default)s").', default='m - mm**2')
    parser_eval_sim.add_argument('--best_score_function', type=str, choices=['max', 'min'], default='max',
                                 help='Specify the function to choose the best score (default: %(default)s).')
    parser_eval_sim.add_argument('--plot_stats', action='store_true',
                                 help='Set this flag to plot #Match ~ #Mismatch distribution.'
                                      'The plot will be stored in the same directory with output statistics file.')
    args = vars(parser.parse_args())
    utility = args.get('utility_name')
    eval(utility)(**args)
