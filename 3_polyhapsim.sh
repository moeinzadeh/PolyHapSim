#!/bin/bash


out_dir="$1"
out_prefix="$2"
contig_len=$3
mutation_rate=$4
ploidy=$5
sequencing_system=$6
# GA1 - GenomeAnalyzer I (36bp,44bp), GA2 - GenomeAnalyzer II (50bp, 75bp)
# HS10 - HiSeq 1000 (100bp),          HS20 - HiSeq 2000 (100bp),      HS25 - HiSeq 2500 (125bp, 150bp)
# HSXn - HiSeqX PCR free (150bp),     HSXt - HiSeqX TruSeq (150bp),   MinS - MiniSeq TruSeq (50bp)
# MSv1 - MiSeq v1 (250bp),            MSv3 - MiSeq v3 (250bp),        NS50 - NextSeq500 v2 (75bp)
read_length=$7
multiple_inserts=${8}
timeout_sec=${9}
scripts_folder=${10}

source ${scripts_folder}/set_env_variables.sh


echo_with_timestamp() {
    date +"[%F %T] $1"
}

echo_with_timestamp "INFO: Starting pipeline. Input arguments:"
echo -e "\tout_dir=${out_dir}"
echo -e "\tout_prefix=${out_prefix}"
echo -e "\tcontig_len=${contig_len}"
echo -e "\tmutation_rate=${mutation_rate}"
echo -e "\tploidy=${ploidy}"
echo -e "\tsequencing_system=${sequencing_system}"
echo -e "\tread_length=${read_length}" 
echo -e "\tlibraries=${multiple_inserts}"
echo -e "\ttimeout_sec=${timeout_sec}"



mkdir -p "${out_dir}"
prefix="${out_dir}/${out_prefix}"



# Generate haplotypes
sim_data_dir="${prefix}_sim_data"
mkdir -p "${sim_data_dir}"




## Make a random DNA sequence template
$FAINX ${SAMPLE_FASTA} random_reference:1-${contig_len} > "${sim_data_dir}"/sample.fasta.tmp
echo ">random_reference_${contig_len}" > "${sim_data_dir}"/sample.fasta
sed '1d' "${sim_data_dir}"/sample.fasta.tmp >> "${sim_data_dir}"/sample.fasta
$SAMTOOLS faidx "${sim_data_dir}"/sample.fasta
rm "${sim_data_dir}"/sample.fasta.tmp

## Run haplogenerator
echo_with_timestamp "INFO: Running haplogenerator.py"
$HAPLOGENERATOR \
    --reference "${sim_data_dir}"/sample.fasta \
    --output "${sim_data_dir}/haplogenerator" \
    --model poisson \
    --subfreq "[${mutation_rate},0,0]" \
    --ploidy ${ploidy} \
|| exit 1


## Map true haplotypes to reference
echo_with_timestamp "INFO: Running minimap2 for mapping true haps to reference"
all_haps="${sim_data_dir}"/"all_true_haps.fa"
all_haps_bam="${sim_data_dir}"/"all_true_haps.bam"
all_haps__sorted_bam="${sim_data_dir}"/"all_true_haps.sorted.bam"
touch $all_haps
for (( i=1; i<=$ploidy; i++ ))
do
    #echo "${sim_data_dir}/haplogenerator_hap"$i".fa"
    cat <(head -1 "${sim_data_dir}/haplogenerator_hap"$i".fa") <(echo "_"$i) | tr -d '\n' >> $all_haps
    echo >> $all_haps
    sed '1d' "${sim_data_dir}/haplogenerator_hap"$i".fa" >> $all_haps
done
$MINIMAP2 -ax asm5 "${sim_data_dir}"/sample.fasta ${all_haps} | $SAMTOOLS view -Sb - > ${all_haps_bam}
$SAMTOOLS sort -o ${all_haps__sorted_bam} -O BAM ${all_haps_bam}
$SAMTOOLS index ${all_haps__sorted_bam}
rm $all_haps $all_haps_bam


# Generate reads
echo_with_timestamp "INFO: Running ART"
touch "${sim_data_dir}"/art.log
for hap in "${sim_data_dir}/haplogenerator_hap"*.fa
do
    IFS='-' read -r -a libraries <<< "$multiple_inserts"
    for library in ${libraries[@]}
    do
        echo ${library}
        IFS='_' read -r -a lib_params <<< "$library"
        fragment_len=${lib_params[0]}
        fragment_len_sd=${lib_params[1]}
        haploid_cov=${lib_params[2]}

        hap_name_tmp=`basename "${hap}" | cut -d"_" -f2 | cut -f1 -d"."`
        echo $hap_name_tmp
        hap_name=${hap_name_tmp}"_"${fragment_len}
        echo $hap_name
        $ART_ILLUMINA \
            -d ${hap_name} \
            -f ${haploid_cov} \
            -i "${hap}" \
            -ir 0.0 \
            -ir2 0.0 \
            -dr 0.0 \
            -dr2 0.0 \
            -m ${fragment_len} \
            -o "${sim_data_dir}/art_${hap_name}" \
            -p \
            -l ${read_length} \
            -s ${fragment_len_sd} \
            -ss ${sequencing_system} \
            -na \
        >> "${sim_data_dir}"/art.log || exit 1
    done
done
# Map reads to the reference
$BWA index "${sim_data_dir}"/sample.fasta
cat "${sim_data_dir}"/art_hap*1.fq > "${sim_data_dir}"/art_1.fq && rm "${sim_data_dir}"/art_hap*1.fq
cat "${sim_data_dir}"/art_hap*2.fq > "${sim_data_dir}"/art_2.fq && rm "${sim_data_dir}"/art_hap*2.fq
echo_with_timestamp "INFO: Running bwa mem"
$BWA mem \
    "${sim_data_dir}"/sample.fasta \
    "${sim_data_dir}"/art_1.fq \
    "${sim_data_dir}"/art_2.fq \
> "${sim_data_dir}"/bwa_mem_aln.sam || exit 1
$SAMTOOLS view -b "${sim_data_dir}"/bwa_mem_aln.sam > "${sim_data_dir}"/bwa_mem_aln.bam
rm "${sim_data_dir}"/bwa_mem_aln.sam
$SAMTOOLS sort -o "${sim_data_dir}"/bwa_mem_aln.bam.sorted "${sim_data_dir}"/bwa_mem_aln.bam
mv "${sim_data_dir}"/bwa_mem_aln.bam.sorted "${sim_data_dir}"/bwa_mem_aln.bam

# Add read group
echo_with_timestamp "INFO: Running AddOrReplaceReadGroups"
java -Xmx8g -jar $PICARD AddOrReplaceReadGroups \
    RGLB="${sequencing_system}_random_reference_${contig_len}_${sim_data_dir}/sample.fasta" \
    RGPL="Illumina" \
    RGID="bwa" \
    RGPU=95 \
    RGSM=L1P1 \
    I="${sim_data_dir}"/bwa_mem_aln.bam \
    O="${sim_data_dir}"/bwa_mem_aln.RG.bam \
    VALIDATION_STRINGENCY=LENIENT \
    QUIET=true || exit 1
mv "${sim_data_dir}"/bwa_mem_aln.RG.bam "${sim_data_dir}"/bwa_mem_aln.bam
$SAMTOOLS index "${sim_data_dir}"/bwa_mem_aln.bam




# Create vcf file with true genotypes
echo_with_timestamp "INFO: Creating VCF file from true haplotype matrix"
python ${scripts_folder}/create_vcf.py "${sim_data_dir}"/haplogenerator_varianthaplos.txt "${sim_data_dir}"/haplogenerator.vcf

# Convert true haps to MHM format
echo_with_timestamp "INFO: Converting true haplotypes matrix into MHM format"
python ${scripts_folder}/true_haps_to_MHM_format.py \
    "${sim_data_dir}"/haplogenerator_varianthaplos.txt \
    ${ploidy} \
    random_reference_${contig_len} \
    "${prefix}"_true_haps_MHM_format.txt



rm "${prefix}"_evaluation.txt
touch "${prefix}"_evaluation.txt



# HapCompass
## Run HapCompass
hapcompass_out_dir="${prefix}_hapcompass_out"
mkdir -p "${hapcompass_out_dir}"
echo_with_timestamp "INFO: Running HapCompass"
/usr/bin/time -v timeout --signal=SIGKILL ${timeout_sec} java -Xmx64g -jar $HAPCOMPASS \
    --ploidy ${ploidy} \
    --iterations 2 \
    --bam "${sim_data_dir}"/bwa_mem_aln.bam \
    --vcf "${sim_data_dir}"/haplogenerator.vcf \
    -o "${hapcompass_out_dir}/${out_prefix}"_hapcompass \
> "${hapcompass_out_dir}/${out_prefix}"_hapcompass.log &&
## Convert hapcompass output to MHM format
{ echo_with_timestamp "INFO: Converting HapCompass output into MHM format";
head -1 "${prefix}"_true_haps_MHM_format.txt > "${prefix}"_hapcompass_output_MHM_format.txt;
python ${scripts_folder}/hapcompass_out_to_MHM_format.py \
    "${hapcompass_out_dir}"/${out_prefix}_hapcompass_MWER_solution.txt \
    ${ploidy} \
    "${prefix}"_hapcompass_output_MHM_format.txt; } &&

{ python ${scripts_folder}/evaluate.py ${ploidy}  ${prefix}"_true_haps_MHM_format.txt" "hapcompass" >> "${prefix}"_evaluation.txt ; }



# H-PoP
hpop_out_dir="${prefix}_hpop_out"
mkdir -p "${hpop_out_dir}"
## Create snp matrix for H-PoP
echo_with_timestamp "INFO: Creating snp matrix for H-PoP"
python ${scripts_folder}/evaluation_scripts.py create_snp_matrix \
    --output_files_prefix "${out_prefix}" \
    --pair_end \
    --filter_supplementary \
    --filter_secondary \
    H-PoPG \
    "${sim_data_dir}"/bwa_mem_aln.bam \
    "${sim_data_dir}"/haplogenerator.vcf \
    "${hpop_out_dir}" \
|| exit 1

## Run H-PoP
echo_with_timestamp "INFO: Running H-PoP"
/usr/bin/time -v timeout --signal=SIGKILL ${timeout_sec} java -Xmx64g -jar $HPOP \
    -p ${ploidy} \
    -f "${hpop_out_dir}/${out_prefix}"_hpopg_snp_matrix.tsv \
    -o "${hpop_out_dir}/${out_prefix}"_hpopg.out \
> "${hpop_out_dir}/${out_prefix}"_hpopg.log &&

## Convert H-PoP output to MHM format
{ echo_with_timestamp "INFO: Converting H-PoP output into MHM format";
head -1 "${prefix}"_true_haps_MHM_format.txt > "${prefix}"_hpop_output_MHM_format.txt;
python ${scripts_folder}/hpop_out_to_MHM_format.py \
    "${hpop_out_dir}/${out_prefix}"_hpopg.out \
    ${ploidy} \
    "${prefix}"_hpop_output_MHM_format.txt; } &&

{ python ${scripts_folder}/evaluate.py ${ploidy}  ${prefix}"_true_haps_MHM_format.txt" "hpop" >> "${prefix}"_evaluation.txt ; }


# SDhaP
sdhap_out_dir="${prefix}_sdhap_out"
mkdir -p "${sdhap_out_dir}"
## Create snp matrix for SDhaP
echo_with_timestamp "INFO: Creating snp matrix for SDhaP (converting H-PoP matrix with FragmentPoly.py)"
${scripts_folder}/FragmentPoly.py \
    -f "${hpop_out_dir}/${out_prefix}"_hpopg_snp_matrix.tsv \
    -o "${sdhap_out_dir}/${out_prefix}"_sdhap_snp_matrix.tsv \
    -x SDhaP

## Run SDhaP
echo_with_timestamp "INFO: Running SDhaP"
/usr/bin/time -v timeout --signal=SIGKILL ${timeout_sec} $SDHAP \
    "${sdhap_out_dir}/${out_prefix}"_sdhap_snp_matrix.tsv \
    "${sdhap_out_dir}/${out_prefix}"_sdhap.out ${ploidy} \
> "${sdhap_out_dir}/${out_prefix}"_sdhap.log &&

## Convert SDhaP output to normal
{ echo_with_timestamp "INFO: Converting SDhaP output back to normal format";
${scripts_folder}/ConvertAllelesSDhaP.py \
    -p "${sdhap_out_dir}/${out_prefix}"_sdhap.out \
    -o "${sdhap_out_dir}/${out_prefix}"_sdhap_norm.out \
    -v "${hpop_out_dir}/${out_prefix}"_random_reference_${contig_len}.vcf; } &&

## Convert normalized SDhaP output to MHM format
{ echo_with_timestamp "INFO: Converting normalized SDhaP output to MHM format";
head -1 "${prefix}"_true_haps_MHM_format.txt > "${prefix}"_sdhap_output_MHM_format.txt;
python ${scripts_folder}/sdhap_out_to_MHM_format.py \
    "${sdhap_out_dir}/${out_prefix}"_sdhap_norm.out \
    ${ploidy} \
    "${prefix}"_sdhap_output_MHM_format.txt; } &&

{ echo_with_timestamp "INFO: Evaluate SdHap results"; python ${scripts_folder}/evaluate.py ${ploidy}  ${prefix}"_true_haps_MHM_format.txt" "sdhap" >> "${prefix}"_evaluation.txt ; }


# Ranbow
ranbow_out_dir="${prefix}_ranbow_out"
mkdir -p "${ranbow_out_dir}"
## Create ranbow parameter file
echo_with_timestamp "INFO: Creating Ranbow parameter file"
ranbow_param_file="${ranbow_out_dir}"/${out_prefix}_ranbow.par
touch "${ranbow_param_file}"
echo "-ploidy ${ploidy}" > "${ranbow_param_file}"
echo "-noProcessor 1" >> "${ranbow_param_file}"
echo "-bamFile ${sim_data_dir}/bwa_mem_aln.bam" >> "${ranbow_param_file}"
echo "-refFile ${sim_data_dir}/sample.fasta" >> "${ranbow_param_file}"
echo "-vcfFile ${sim_data_dir}/haplogenerator.vcf" >> "${ranbow_param_file}"
echo "-outputFolderBase ${ranbow_out_dir}" >> "${ranbow_param_file}"
echo "random_reference_${contig_len}" > "${ranbow_out_dir}"/contigs.list
echo "-selectedScf ${ranbow_out_dir}/contigs.list" >> "${ranbow_param_file}"


## Run Ranbow, Run! Hahahahaha
echo_with_timestamp "INFO: Running Ranbow..."
python ${RANBOW} hap -mode index -par "${ranbow_param_file}" &&
{ /usr/bin/time -v timeout --signal=SIGKILL ${timeout_sec} python ${RANBOW} hap \
    -mode hap \
    -par "${ranbow_param_file}" \
> "${ranbow_out_dir}"/${out_prefix}_ranbow.log;
python $RANBOW hap -mode collect -par "${ranbow_param_file}" ; } &&


#convert Ranbow output to correct format
{ echo_with_timestamp "run convert_ranbow_out.py";
echo ${ranbow_out_dir}'/ranbow.single.hap';
echo ${out_dir}"/test_ranbow_output_MHM_format.txt";
python ${scripts_folder}/convert_ranbow_out.py  ${ranbow_out_dir}'/ranbow.single.hap' ${prefix}"_ranbow_output_MHM_format.txt" ; } &&

{ echo_with_timestamp "evaluate.py";
python ${scripts_folder}/evaluate.py ${ploidy}  ${prefix}"_true_haps_MHM_format.txt" "ranbow" >> "${prefix}"_evaluation.txt ; } &&
{ echo_with_timestamp "RanbowMEC:";
 grep -v ">" ${ranbow_out_dir}'/ranbow.single.hap' | cut -f 5 |  paste -sd+ | bc ; }

echo "--------------------------"
cat "${prefix}"_evaluation.txt
