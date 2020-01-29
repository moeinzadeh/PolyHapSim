result_folder = 'pathto???'
scriptFolder = 'pathto???'
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
errors =  ['0.001','0.005','0.01','0.05','0.1']




import subprocess
import os


collect_results = result_folder+'/G'+genomelen+'_t'+timeout_each_method+'s_Ct'+runTimeOnCluster +"_M" + memoryOnCluster




for error in errors[:-1]:
    for ploidy in range(4,9,2):
        for replicate in range(replicates):
            st = "bash "+scriptFolder + "/1_cluster_polyhapsim.sh "+ \
             multiple_inserts[:-1]+ ' '  + error + ' ' + genomelen + ' ' +  timeout_each_method + ' ' + str(replicate) + " "+ \
                 collect_results + " " + runTimeOnCluster + " " + memoryOnCluster + " " + grp_indentifier + " " + str(ploidy) + " "  + scriptFolder
            os.system(st)
os.system("echo " + multiple_inserts + '> ' + collect_results + "/readme.txt")
