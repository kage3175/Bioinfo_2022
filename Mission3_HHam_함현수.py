chr_list=[str(i) for i in range(1,23)]
chr_list.append('X')
chr_list.append('Y')

class RefSeq:
    def __init__(self):
        self.RefSeqID='NULL'
        self.Gene_Symbol='NULL'
        self.ChrID=0
        self.Strand='NULL'
        self.Num_Exons=0
        self.Txn_start=0
        self.Txn_end=0
        self.Coding_start=0
        self.Coding_end=0
        self.Exon_num=0
        self.Exon_starts=[]
        self.Exon_ends=[]
        self.NUM_RefSeqID=0
        self.temp=[]
    def exon_list_precessing(self,list_exon):
        temp=list_exon.split(",")
        temp.pop()
        return list(map(int,temp))
    def parsing(self, sLine):
        list_sLine=sLine.split('\t')
        self.Gene_Symbol=list_sLine[0]
        self.RefSeqID=list_sLine[1]
        self.ChrID=list_sLine[2][3:]
        self.Strand=list_sLine[3]
        self.Txn_start=int(list_sLine[4])
        self.Txn_end=int(list_sLine[5])
        self.Coding_start=int(list_sLine[6])
        self.Coding_end=int(list_sLine[7])
        self.Exon_num=int(list_sLine[8])
        self.Exon_starts=self.exon_list_precessing(list_exon=list_sLine[9])
        self.Exon_ends=self.exon_list_precessing(list_exon=list_sLine[10])
        self.NUM_RefSeqID=int(self.RefSeqID[3:])
        
def main():
    global chr_list
    #test=RefSeq()
    #temp_RefSeq=RefSeq()
    list_RefSeq_raw=[]
    list_RefSeq_NM=[]
    list_RefSeq_SingleEntry=[]
    dict_ID={}
    file=open("../files_bioinfo2022/refFlat.txt", 'r')
    for sLine in file.readlines():
        temp_RefSeq=RefSeq()
        temp_RefSeq.parsing(sLine)
        list_RefSeq_raw.append(temp_RefSeq)
        if(temp_RefSeq.RefSeqID[:3]=='NM_' and (temp_RefSeq.ChrID in chr_list)):
            list_RefSeq_NM.append(temp_RefSeq)
            try:
                dict_ID[temp_RefSeq.RefSeqID]+=1
            except KeyError:
                dict_ID[temp_RefSeq.RefSeqID]=1
    file.close()
    outfile=open("../files_bioinfo2022/result.txt", 'w')
    for refseq in list_RefSeq_NM:
        if dict_ID[refseq.RefSeqID]==1:
            list_RefSeq_SingleEntry.append(refseq)
    print(len(list_RefSeq_raw), len(list_RefSeq_NM), len(list_RefSeq_SingleEntry))
    list_RefSeq_SingleEntry.sort(key=lambda x:x.NUM_RefSeqID)
    for refseq in list_RefSeq_SingleEntry:
        print(refseq.RefSeqID+'\t'+refseq.Gene_Symbol, file=outfile)
    outfile.close()
    
main()