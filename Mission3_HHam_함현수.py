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

def make_dict(list_RefSeq_NM):# 각 entry가 1개 이상인지를 확인하기 위한 딕셔너리를 만들어서 반환
    dict={}
    for refseq in list_RefSeq_NM:
        try:
            dict[refseq.RefSeqID]+=1
        except KeyError:
            dict[refseq.RefSeqID]=1
    # End of for body for refseq
    return dict
######################################################################## End of make_dict

def delete_multientry(dict_check, list_RefSeq_NM): # entry가 여러 개인 경우를 제외하고 만든 (클래스)리스트를 반환
    templist=[]
    for refseq in list_RefSeq_NM:
        if dict_check[refseq.RefSeqID]==1:
            templist.append(refseq)
    # End of for body for refseq
    return templist
######################################################################## End of delete_multientry

def main():
    global chr_list
    start=time.time()
    list_RefSeq_raw=[] # 클래스들을 담는 리스트. 해당 gene이 어떤 특성을 갖던 일단 모두 넣고 본다.
    list_RefSeq_NM=[] # RefSeqID가 NM_으로 시작하고 1~22, X, Y chromosome에 존재하는 클래스들만 담을 리스트
    file=open("../files_bioinfo2022/refFlat.txt", 'r') ############################## 제출할 때 바꿔야함
    for sLine in file.readlines():
        temp_RefSeq=RefSeq()
        temp_RefSeq.parsing(sLine)
        list_RefSeq_raw.append(temp_RefSeq)
        if(temp_RefSeq.RefSeqID[:3]=='NM_' and (temp_RefSeq.ChrID in chr_list)): #RefSeqID가 NM_으로 시작하고 위치하는 chromosome이 chr_list(1~22,X,Y)에 존재하는 경우만 리스트에 append 해줌
            list_RefSeq_NM.append(temp_RefSeq)
    # End of for body for sLine
    dict_ID=make_dict(list_RefSeq_NM)# entry가 1개 이상인지를 확인하기 위한 딕셔너리
    file.close()
    outfile=open("../files_bioinfo2022/result.txt", 'w') ############################## 제출할 때 바꿔야함
    list_RefSeq_SingleEntry=delete_multientry(dict_ID, list_RefSeq_NM) # single entry gene만 골라 담은, 최종적인 리스트
    print(len(list_RefSeq_raw), len(list_RefSeq_NM), len(list_RefSeq_SingleEntry))
    list_RefSeq_SingleEntry.sort(key=lambda x:x.NUM_RefSeqID)#각 RefSeqID의 뒤에 숫자에 대해서 sorting
    for refseq in list_RefSeq_SingleEntry:
        print(refseq.RefSeqID+'\t'+refseq.Gene_Symbol, file=outfile)
    # End of for body for refseq
    outfile.close()
    print("총 걸린 시간은 "+str(time.time()-start)+"초입니다.")
######################################################################## End of main

main()