import sys
fi = open(sys.argv[1],'r')
fo = open(sys.argv[2],'w')
lines = fi.readlines()

haplen = int(lines[0].split()[2])
fo.write(lines[0])
inx = 1
while inx < len(lines):
    if lines[inx][0]==">":
        fo.write(lines[inx])
    else:
        a = lines[inx].split()
        fo.write(a[5].split("_")[0]+"\t"+ a[1]+"\t"+ a[2]+"\t"+ a[3]+'\n')
    inx +=1
