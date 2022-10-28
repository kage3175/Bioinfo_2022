import time

chr_list=[str(i) for i in range(1,23)]
chr_list.append('X')
chr_list.append('Y')#Chr list를 담기위한 것

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
    # End of __init__
    def exon_list_precessing(self,list_exon):
        temp=list_exon.split(",")
        temp.pop()
        return list(map(int,temp))
    # End of exon_list_processing
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
    # End of parsing
######################################################################## End of RefSeq

def make_dict(list_RefSeq_NM):
    dict={}
    for refseq in list_RefSeq_NM:
        try:
            dict[refseq.RefSeqID]+=1
        except KeyError:
            dict[refseq.RefSeqID]=1
    return dict
def delete_multientry(dict_check, list_RefSeq_NM):
    templist=[]
    for refseq in list_RefSeq_NM:
        if dict_check[refseq.RefSeqID]==1:
            templist.append(refseq)
    return templist

def main():
    start=time.time()
    global chr_list
    list_RefSeq_raw=[]
    list_RefSeq_NM=[]
    file=open("../files_bioinfo2022/refFlat.txt", 'r')
    for sLine in file.readlines():
        temp_RefSeq=RefSeq()
        temp_RefSeq.parsing(sLine)
        list_RefSeq_raw.append(temp_RefSeq)
        if(temp_RefSeq.RefSeqID[:3]=='NM_' and (temp_RefSeq.ChrID in chr_list)):
            list_RefSeq_NM.append(temp_RefSeq)
    dict_ID=make_dict(list_RefSeq_NM)# entry가 1개 이상인지를 확인하기 위한 딕셔너리
    file.close()
    outfile=open("../files_bioinfo2022/result.txt", 'w')
    list_RefSeq_SingleEntry=delete_multientry(dict_ID, list_RefSeq_NM)
    print(len(list_RefSeq_raw), len(list_RefSeq_NM), len(list_RefSeq_SingleEntry))
    list_RefSeq_SingleEntry.sort(key=lambda x:x.NUM_RefSeqID)#각 RefSeqID의 뒤에 숫자에 대해서 sorting
    for refseq in list_RefSeq_SingleEntry:
        print(refseq.RefSeqID+'\t'+refseq.Gene_Symbol, file=outfile)
    outfile.close()
    print("총 걸린 시간은 "+str(time.time()-start)+"초입니다.")
######################################################################## End of main

main()