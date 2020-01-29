#!/usr/bin/env bash
multipleInserts=$1
error=$2
genomeLen=$3
timeout=$4
replicate=$5
working_dir_result=$6
runTimeOnCluster=$7
memoryOnCluster=$8
grpidentifier=$9
ploidy=${10}
scriptFolder=${11}


mkdir -p ${working_dir_result}

id=$9"_"$2"_"$3"_"$4"_"$5"_"$7"_"$8"_"${10}
groupID=$9"_"$2"_"$3"_"$4"_"$7"_"$8"_"${10}
working_dir_haploEvgeny='/scratch/local2/moeinzadeh/'$id
stdOut=${working_dir_result}"/"$id".out"
stdErr=${working_dir_result}"/"$id".err"



bash ${scriptFolder}/2_run_polyhapsim.sh ${working_dir_haploEvgeny} ${genomeLen} ${error} ${multipleInserts} ${timeout} ${replicate} ${working_dir_result} ${grpidentifier} ${ploidy} ${scriptFolder}

