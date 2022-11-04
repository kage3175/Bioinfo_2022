import time

chr_list=[str(i) for i in range(1,23)]
chr_list.append('X')
chr_list.append('Y')#Chr list를 담기위한 것. 사람의 유전자 이름은 절대 바뀔 일이 없으므로 전역변수로 설정하였다.

BASE_COMPLE={'A':'T','C':'G','G':'C','T':'A'} #상보적인 염기쌍을 전역변수로서 미리 잡아둔다. A -T, C- G 간 결합을 딕셔너리로 표현한 것
STOP_CODON=['TAA','TGA', 'TAG']

class RefSeq:
    global BASE_COMPLE
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
        self.mRNASeq='NULL'
        self.ORF_start=0
        self.ORF_end=0
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
        self.Length_mRNA=0
        self.Length_ORF=0
    # End of parsing
    def Get_mRNASeq(self,dict_file_chr):
        mRNA_seq=''
        temp_length=0
        chr_to_search=dict_file_chr[self.ChrID]
        for i in range(self.Exon_num):
            mRNA_seq=mRNA_seq+chr_to_search[self.Exon_starts[i]:self.Exon_ends[i]]#각 엑손의 시작~끝부분(시작점-1:끝점)을 문자열 슬라이싱을 통해 붙여나간다. ####주의: refseq의 location 정보는 +1을 시작으로 삼는다.
            if(self.Exon_starts[i]<=self.Coding_start<self.Exon_ends[i]):
                self.ORF_start=temp_length+self.Coding_start-self.Exon_starts[i]#이거 포함. + 기준 이건 무조건 A여야함
            if(self.Exon_starts[i]<=self.Coding_end<self.Exon_ends[i]):
                self.ORF_end=temp_length+self.Coding_end-self.Exon_starts[i]#이거 제외. + 기준 이거 바로 앞 3문자가 stop codod이여야함
            temp_length+=self.Exon_ends[i]-self.Exon_starts[i]
        self.Length_mRNA=temp_length
        #End of for body for i
        if self.Strand=='-': #strand가 '-'일 경우에는 따로 tempstr에서 서열을 뒤집고 상보적인 염기쌍으로 반전시켜주고 바로 리턴한다. '+' strand일 경우에는 이 if 문을 무시하고 mRNA_seq를 반환한다.
            tempstr=''
            for i in range(len(mRNA_seq)-1,-1,-1):
                tempstr=tempstr+BASE_COMPLE[mRNA_seq[i]]
            self.mRNASeq=tempstr
            temp=self.ORF_end
            self.ORF_end=self.Length_mRNA-1-self.ORF_start+1
            self.ORF_start=self.Length_mRNA-temp
        else:
            self.mRNASeq=mRNA_seq
        self.Length_ORF=self.ORF_end-self.ORF_start
            # End of for body for i
        # End of Get_mRNASeq  
                
######################################################################## End of RefSeq

def file_processing(sFile):#문자열 file을 건네받고 문자열 내부의 문자들은 전부 대문자로 바꾸고, 개행문자를 삭제함
    newsFile=sFile.replace('\n', '')
    newsFile=newsFile.upper()
    return newsFile
######################################################################## End of file_processing

'''def return_seq(refseq, dict_file_chr):
    global BASE_COMPLE
    mRNA_seq=''
    chr_to_search=dict_file_chr[refseq.ChrID]
    for i in range(refseq.Exon_num):
        mRNA_seq=mRNA_seq+chr_to_search[refseq.Exon_starts[i]:refseq.Exon_ends[i]]#각 엑손의 시작~끝부분(시작점-1:끝점)을 문자열 슬라이싱을 통해 붙여나간다. ####주의: refseq의 location 정보는 +1을 시작으로 삼는다.
    # End of for body for i
    if refseq.Strand=='-': #strand가 '-'일 경우에는 따로 tempstr에서 서열을 뒤집고 상보적인 염기쌍으로 반전시켜주고 바로 리턴한다. '+' strand일 경우에는 이 if 문을 무시하고 mRNA_seq를 반환한다.
        tempstr=''
        for i in range(len(mRNA_seq)-1,-1,-1):
            tempstr=tempstr+BASE_COMPLE[mRNA_seq[i]]
        # End of for body for i
        return tempstr
    return mRNA_seq
######################################################################## End of return_seq'''

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

def make_dict_2(list_RefSeq_validORF):# 각 entry가 1개 이상인지를 확인하기 위한 딕셔너리를 만들어서 반환
    dict={}
    for refseq in list_RefSeq_validORF:
        try:
            dict[refseq.Gene_Symbol]+=1
        except KeyError:
            dict[refseq.Gene_Symbol]=1
    # End of for body for refseq
    return dict
######################################################################## End of make_dict_2

def delete_multientry(dict_check, list_RefSeq_NM): # entry가 여러 개인 경우를 제외하고 만든 (클래스)리스트를 반환
    templist=[]
    for refseq in list_RefSeq_NM:
        if dict_check[refseq.RefSeqID]==1:
            templist.append(refseq)
    # End of for body for refseq
    return templist
######################################################################## End of delete_multientry
def delete_multientry_2(dict_check, list_RefSeq_NM): # entry가 여러 개인 경우를 제외하고 만든 (클래스)리스트를 반환
    templist=[]
    for refseq in list_RefSeq_NM:
        if dict_check[refseq.Gene_Symbol]==1:
            templist.append(refseq)
            dict_check[refseq.Gene_Symbol]=0
        elif dict_check[refseq.Gene_Symbol]>1:
            templist.append(refseq)
            dict_check[refseq.Gene_Symbol]=0
    # End of for body for refseq
    return templist
######################################################################## End of delete_multientry_2

def print_outfile(list_RefSeq, outfile):
    for refseq in list_RefSeq:
        #mRNA_seq=return_seq(refseq, dict_file_chr)
        print(refseq.RefSeqID, '\n', refseq.mRNASeq,file=outfile)
    #End for body for refseq
######################################################################## End of print_outfile

def check_valid_mRNA(refseq):
    global STOP_CODON
    if(refseq.Length_ORF%3!=0):
        return False
    elif(refseq.mRNASeq[refseq.ORF_start:refseq.ORF_start+3]!='ATG'):
        return False
    elif(refseq.mRNASeq[refseq.ORF_end-3:refseq.ORF_end] not in STOP_CODON):
        return False
    else:
        for i in range(refseq.Length_ORF//3 -1):#중간에 STOP CODON이 존재하는지 검사
            if(refseq.mRNASeq[refseq.ORF_start+i*3:refseq.ORF_start+i*3+3] in STOP_CODON):
                return False
        return True# 다 통과하면 True
######################################################################## End of check_valid_mRNA


def main():
    global chr_list
    start=time.time()
    list_RefSeq_raw=[] # 클래스들을 담는 리스트. 해당 gene이 어떤 특성을 갖던 일단 모두 넣고 본다.
    list_RefSeq_NM=[] # RefSeqID가 NM_으로 시작하고 1~22, X, Y chromosome에 존재하는 클래스들만 담을 리스트
    list_validORF=[]#ANS4에 해당하는 refseq만 남기는 리스트
    #mRNA_gene={}
    file=open("../files_bioinfo2022/refFlat.txt", 'r')############################## 제출할 때 바꿔야함
    for sLine in file.readlines():
        temp_RefSeq=RefSeq()
        temp_RefSeq.parsing(sLine)
        list_RefSeq_raw.append(temp_RefSeq)
        if(temp_RefSeq.RefSeqID[:3]=='NM_' and (temp_RefSeq.ChrID in chr_list)):
            list_RefSeq_NM.append(temp_RefSeq)
    # End of for body for sLine
    dict_ID=make_dict(list_RefSeq_NM)# entry가 1개 이상인지를 확인하기 위한 딕셔너리
    file.close()
    #outfile=open("../files_bioinfo2022/result.txt", 'w')############################## 제출할 때 바꿔야함
    list_RefSeq_SingleEntry=delete_multientry(dict_ID, list_RefSeq_NM)
    list_RefSeq_SingleEntry.sort(key=lambda x:x.NUM_RefSeqID)#각 RefSeqID의 뒤에 숫자에 대해서 sorting
    '''for refseq in list_RefSeq_SingleEntry:
        print(refseq.RefSeqID+'\t'+refseq.Gene_Symbol, file=outfile) #result.txt 파일에 각 RefSeqID에 대한 gene symbol을 탭으로 구분해서 출력
    # End of for body for refseq'''
    #outfile.close()
    #outfile=open("../files_bioinfo2022/result_mRNAs_test.txt", 'w')############################## 제출할 때 바꿔야함
    dict_file_chr={}
    for chr_num in chr_list: # 1번부터 Y까지의 크로모좀 전체 시퀀스를 해당 크로모좀 번호(string 형태)를 키로 하는 딕셔너리에 저장
        temp_file=open("../files_bioinfo2022/hg38ChrFiles/chr"+chr_num+".fa","r")############################## 제출할 때 바꿔야함
        trash=temp_file.readline()
        temp_seq_chr=temp_file.read()
        temp_seq_chr=file_processing(temp_seq_chr)
        dict_file_chr[chr_num]=temp_seq_chr
    for refseq in list_RefSeq_SingleEntry:
        refseq.Get_mRNASeq(dict_file_chr)
        if(check_valid_mRNA(refseq)):
            list_validORF.append(refseq)
    dict_ID_isoform=make_dict_2(list_validORF)
    list_mRNA_ANS5=delete_multientry_2(dict_ID_isoform,list_validORF)
    print(len(list_RefSeq_raw), len(list_RefSeq_NM), len(list_RefSeq_SingleEntry), len(list_validORF),len(list_mRNA_ANS5))
    
    '''for refseq in list_RefSeq_SingleEntry:
        print_outfile(list_RefSeq_SingleEntry,outfile)
        print(refseq.RefSeqID, refseq.ORF_start, refseq.ORF_end, len(refseq.mRNASeq), refseq.Length_mRNA)
        print(refseq.mRNASeq[refseq.ORF_start:refseq.ORF_start+3], refseq.mRNASeq[refseq.ORF_end-3:refseq.ORF_end])'''
    
    #print_outfile(list_RefSeq_SingleEntry,dict_file_chr,outfile)#outfile에 txt 형태로 각 RefSeqID에 해당하는 mRNA를 write하는 함수.
    #outfile.close()
    print("총 걸린 시간은 "+str(time.time()-start)+"초입니다.")
######################################################################## End of main

main()