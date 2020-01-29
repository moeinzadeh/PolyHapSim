working_dir=$1
genomeLen=$2
error=$3
multiple_inserts=$4
timeout=$5
replicate=$6
workingDir_mxqsub=$7
grpidentifier=$8
ploidy=$9
scriptFolder=${10}


rm -rf ${working_dir}
mkdir -p ${working_dir}

prefix=${workingDir_mxqsub}'/'${grpidentifier}"_"${error}"_"${genomeLen}"_"${replicate}"_"${ploidy}

bash ${scriptFolder}/3_polyhapsim.sh ${working_dir} haplosim ${genomeLen} ${error} ${ploidy} MSv3 150 ${multiple_inserts} ${timeout} ${scriptFolder}
cp ${working_dir}'/haplosim_evaluation.txt' ${prefix}"_evaluation.txt"
cp ${working_dir}'/haplosim_ranbow_out/haplosim_ranbow.log'  ${prefix}"_ranbow.txt"
cp ${working_dir}'/haplosim_sdhap_out/haplosim_sdhap.log'  ${prefix}"_sdhap.txt"
cp ${working_dir}'/haplosim_hapcompass_out/haplosim_hapcompass.log'  ${prefix}"_hapcompass.txt"
cp ${working_dir}'/haplosim_hpop_out/haplosim_hpopg.log'  ${prefix}"_hpopg.txt"
#rm -r ${working_dir}
