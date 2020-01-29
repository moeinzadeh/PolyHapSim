cd $1
##match mismatch
grep "\-\-" *evaluation.txt | tr "_" " "  | cut -f1,2,3,4,5,8,9,10,11,12 -d' ' | tr -s " " "\t" > all_methods_reasults_M_MM.txt


#MEC
grep -A1 Time *sdhap.txt | grep MEC | tr '_' " " | sed "s:sdhap.txt-Total MEC of all blocks=: :g" | tr "," " " | tr -s " " | cut -d' ' -f 1,2,3,4,5,6| tr -s " " "\t" | awk '{print $0 "\t" "sdhap"}' > all_methods_reasults_MEC.txt


#time
grep Time *sdhap.txt | grep haplotype | tr "_" " " | tr -d ":" | sed "s:sdhap.txtTime to haplotype:sdhap:g" | tr " " "\t" > all_methods_reasults_time.txt
grep "Phased\:" *hpopg.txt | tr "_" " " | sed "s:hpopg.txt:hpopg:g" | tr ":" " " |tr -s " " | tr " " "\t" | cut -f 1,2,3,4,5,6,14 >> all_methods_reasults_time.txt
grep readBAM *ranbow.txt | tr "_\t" " " | sed "s:ranbow.txt:Ranbow :g" | cut -d' ' -f1,2,3,4,5,6,14 | tr " " "\t" >> all_methods_reasults_time.txt
grep MWER  *_hapcompass.txt | tr "_" " " | tr -s " " | sed "s:hapcompass.txt:hapcompass :g" | cut -f1,2,3,4,5,6,8 -d' ' | tr " " "\t" >> all_methods_reasults_time.txt
grep RanbowMEC -A 1  *.out | grep -v RanbowMEC | grep -v "\-\-" | tr "_-" " " | cut -f 1,2,3,8,9 -d' ' | tr " " "\t" | sed 's/.out/ /g' > ranbow_error_correction.txt
egrep "Maximum resident set size|Command being timed" *.err | tr "_" " "  | cut -f1,2,3,4,5,6,10,11,12,13,14,15,16,17,18,8  -d' ' | tr -s " " "\t" > all_methods_reasults_memory.txt
