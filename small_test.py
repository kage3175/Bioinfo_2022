file=open("../files_bioinfo2022/mature.fa", "r")
count=0
for sLine in file.readlines():
    if(sLine[:4]==">hsa"):
        count+=1
        
print(count)


file.close()