# PolyHapSim: A haplotype simulator for polyploid genomes

PolyHapSim takes genome length, heterozygousity rates, ploidy levels, and set of illumina library insert sizes and their coverage to simulate haplotypes and Illumina reads for evaluation.
It takes paths to installed Ranbow, H-PoP, SDhaP, HapCompass tools in order to evaluates them.

### Set parameters in `main.py`

```
result_folder = '/project/moeinzadeh-data/tmp'
scriptFolder = '/project/moeinzadeh-data/project/polyploid_haplotyping/script_github_main'
genomelen = '100000'
timeout_each_method = '10000'
runTimeOnCluster = '20h'
memoryOnCluster = '120G'
replicates = 5
grp_indentifier='t7'
coverage = '30'
multiple_inserts=""
for insert in [('350','50',coverage),('500','100',coverage),('1000','200',coverage),('2000','400',coverage),('5000','500',coverage)]:
    multiple_inserts+=insert[0]+"_"+insert[1]+"_"+insert[2]+"-"
heterozygousity_rates =  ['0.001','0.005','0.01','0.05','0.1']
ploidies = range(4,9,2)
```

### Set env variables in `set_env_variables.sh`
```

#!/usr/bin/env bash
export FAINX=PATHTO/faidx
export HAPLOGENERATOR=PATHTO/haplogenerator.py
export MINIMAP2=PATHTO/minimap2
export ART_ILLUMINA=/PATHTO/art_illumina
export BWA=PATHTO/bwa
export PICARD=PATHTO/picard.jar
export SAMTOOLS=PATHTO/samtools
export HAPCOMPASS=PATHTO/hapcompass.jar
export HPOP=PATHTO/H-PoPGv0.2.0.jar
export SDHAP=PATHTO/hap
export RANBOW=PATHTO/ranbow.py
export SAMPLE_FASTA=PATHTO/sample.fasta


```


Other contributors to this project:
Eugene Muzychenko (github id: emuzychenko)
