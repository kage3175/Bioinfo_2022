import time

from scipy import stats


chr_list=[str(i) for i in range(1,23)]
chr_list.append('X')
chr_list.append('Y')#Chr list를 담기위한 것. 사람의 유전자 이름은 절대 바뀔 일이 없으므로 전역변수로 설정하였다.
chr_list.sort()

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
        self.ORF_start=-1
        self.ORF_end=-1    #ORF_end는 3'UTR start와 일치한다. end 자체가 exclusive이기 때문
        self.Length_3UTR=0
        self.Length_mRNA=0
        self.Length_ORF=0
        self.dict_motif={}
    # End of __init__
    def exon_list_processing(self,list_exon):
        temp=list_exon.split(",")
        temp.pop()
        return list(map(int,temp))
    # End of exon_list_processing
    def parsing(self, sLine):
        list_sLine=sLine.split('\t')
        self.Gene_Symbol=list_sLine[0].upper()
        self.RefSeqID=list_sLine[1]
        self.ChrID=list_sLine[2][3:]
        self.Strand=list_sLine[3]
        self.Txn_start=int(list_sLine[4])
        self.Txn_end=int(list_sLine[5])
        self.Coding_start=int(list_sLine[6])
        self.Coding_end=int(list_sLine[7])
        self.Exon_num=int(list_sLine[8])
        self.Exon_starts=self.exon_list_processing(list_exon=list_sLine[9])
        self.Exon_ends=self.exon_list_processing(list_exon=list_sLine[10])
        self.NUM_RefSeqID=int(self.RefSeqID[3:])
    # End of parsing
    def Get_mRNASeq(self,file_chr): #chromosome 전체 sequence를 넣어주면 mRNA 서열을 저장하는 메소드
        mRNA_seq=''
        temp_length=0
        chr_to_search=file_chr
        for i in range(self.Exon_num):
            mRNA_seq=mRNA_seq+chr_to_search[self.Exon_starts[i]:self.Exon_ends[i]]#각 엑손의 시작~끝부분(시작점-1:끝점)을 문자열 슬라이싱을 통해 붙여나간다. ####주의: refseq의 location 정보는 +1을 시작으로 삼는다.
            if(self.Exon_starts[i]<=self.Coding_start<self.Exon_ends[i]):
                self.ORF_start=temp_length+self.Coding_start-self.Exon_starts[i]#이거 포함. + 기준 이건 무조건 A여야함
            if(self.Exon_starts[i]<self.Coding_end<=self.Exon_ends[i]):
                self.ORF_end=temp_length+self.Coding_end-self.Exon_starts[i]#이거 제외. + 기준 이거 바로 앞 3문자가 stop codod이여야함
            temp_length+=self.Exon_ends[i]-self.Exon_starts[i]
        self.Length_mRNA=temp_length#Length mRNA 변수에 길이 저장
        #End of for body for i
        if self.Strand=='-': #strand가 '-'일 경우에는 따로 tempstr에서 서열을 뒤집고 상보적인 염기쌍으로 반전시켜주고 바로 리턴한다. '+' strand일 경우에는 이 if 문을 무시하고 mRNA_seq를 반환한다.
            tempstr=''
            for i in range(len(mRNA_seq)-1,-1,-1):
                tempstr=tempstr+BASE_COMPLE[mRNA_seq[i]]
            self.mRNASeq=tempstr
            self.ORF_end, self.ORF_start=self.Length_mRNA-self.ORF_start, self.Length_mRNA-self.ORF_end # ORF end와 ORF start도 - strand에 맞게 변환시켜준다.
        else:
            self.mRNASeq=mRNA_seq
        self.Length_ORF=self.ORF_end-self.ORF_start#ORF 길이 저장
        self.Length_3UTR=self.Length_mRNA-self.ORF_end
    def getMotif_3UTR(self):
        if(self.Length_3UTR>=7):
            for i in range(self.ORF_end, self.Length_mRNA-6):
                self.dict_motif[self.mRNASeq[i:i+7]]=1
    #End of Get_mRNASeq
######################################################################## End of RefSeq

def file_processing(sFile):#문자열 file을 건네받고 문자열 내부의 문자들은 전부 대문자로 바꾸고, 개행문자를 삭제함
    newsFile=sFile.replace('\n', '')
    newsFile=newsFile.upper()
    return newsFile
######################################################################## End of file_processing

def make_dict_entry(list_RefSeq_NM):# 각 entry가 1개 이상인지를 확인하기 위한 딕셔너리를 만들어서 반환
    dict={}
    for refseq in list_RefSeq_NM:
        try:
            dict[refseq.RefSeqID]+=1
        except KeyError:
            dict[refseq.RefSeqID]=1
    # End of for body for refseq
    return dict
######################################################################## End of make_dict_entry

def make_dict_for_isoform(list_RefSeq_validORF):# 다른 entry여도 Gene_Symbol이 여러 개인 경우를 추려내기 위한 딕션너리를 반환한다. 딕셔너리는 Gene_Symbol을 키로, 개수 상관없이 무조건 1로 저장된다.
    dict={}
    for refseq in list_RefSeq_validORF:
        dict[refseq.Gene_Symbol]=1
    # End of for body for refseq
    return dict
######################################################################## End of make_dict_for_isoform

def delete_multientry(dict_check, list_RefSeq_NM): # entry가 여러 개인 경우를 제외하고 만든 (클래스)리스트를 반환
    templist=[]
    for refseq in list_RefSeq_NM:
        if dict_check[refseq.RefSeqID]==1:
            templist.append(refseq)
    # End of for body for refseq
    return templist
######################################################################## End of delete_multientry

def leave_representitive_isoform(dict_check, list_RefSeq_NM): # entry가 여러 개인 isoform 중 RefSeqID 숫자가 가장 작은 하나만 남겨서 리스트 반환
    templist=[]
    for refseq in list_RefSeq_NM:
        if dict_check[refseq.Gene_Symbol]==1: #만약 해당 Gene_Symbol에 해당하는 refseq class가 이미 리스트에 추가되어 있다면, 아래 코드에 의해 딕셔너리 값이 0이 되므로, 같은 Gene_Symbol을 가지는 refseq class는 2번 이상 리스트에 추가될 수 없다.
            templist.append(refseq)
            dict_check[refseq.Gene_Symbol]=0
    # End of for body for refseq
    return templist
######################################################################## End of leave_representitive_isoform

def print_outfile(list_RefSeq, outfile):#Excel에 적을 답을 txt 파일에 출력해주는 함수
    for refseq in list_RefSeq:
        print(refseq.RefSeqID+'\t'+refseq.Gene_Symbol+'\t'+str(refseq.ORF_start)+'\t'+str(refseq.Length_ORF)+'\t'+str(refseq.Length_mRNA-refseq.ORF_end) ,file=outfile)
    #End for body for refseq
######################################################################## End of print_outfile

def check_valid_mRNA(refseq): # mRNA의 ORF 길이가 3의 배수인지, start codon이나 stop codon이 시작과 끝에 있는지, 중간에 멈춰버리지 않는 지 등을 검사
    global STOP_CODON
    if(refseq.Length_ORF%3!=0): # ORF의 길이가 3의 배수인지
        return False
    elif(refseq.mRNASeq[refseq.ORF_start:refseq.ORF_start+3]!='ATG'): # ORF의 첫 3글자가 ATG인지
        return False
    elif(refseq.mRNASeq[refseq.ORF_end-3:refseq.ORF_end] not in STOP_CODON): # ORF의 마지막 3글자가 종결코돈인지
        return False
    else:
        for i in range(refseq.Length_ORF//3 -1):#중간에 STOP CODON이 존재하는지 검사. 끝은 어차피 stop codon이어야 하므로 검사하지 않는다.
            if(refseq.mRNASeq[refseq.ORF_start+i*3:refseq.ORF_start+i*3+3] in STOP_CODON):
                return False
        #End of for body for i
        return True# 다 통과하면 True
######################################################################## End of check_valid_mRNA

def make_list_valid(list_RefSeq_SingleEntry):
    global chr_list
    templist=[]
    list_RefSeq_SingleEntry.sort(key=lambda x:x.ChrID)#ChrID 순서대로 먼저 정렬. 이러지 않으면 for문을 5만번씩 24번 돌려야해서 비효율 적이다. 정렬하면 for문 5만번 1번만 돌리면 된다.
    length_list=len(list_RefSeq_SingleEntry)
    cnt=0
    for chr_num in chr_list: # 1번부터 Y까지의 크로모좀 전체 시퀀스를 해당 크로모좀 번호(string 형태)를 키로 하는 딕셔너리에 저장
        chr_file=open("hg38ChrFiles/chr"+chr_num+".fa","r")############################## 제출할 때 바꿔야함
        trash=chr_file.readline()
        full_seq_chr=chr_file.read()
        full_seq_chr=file_processing(full_seq_chr)
        chr_file.close()
        while cnt<length_list: # RefSeq이 담긴 리스트의 길이보다 더 검사할 필요도 없고, 하면 오류가 나니까 끊어준다.
            if(list_RefSeq_SingleEntry[cnt].ChrID!=chr_num):# ChrID가 현재 Chr과 다르면. 이미 ChrID 순서대로 sorting 했으므로, 다음 Chr을 검사해야할 시점이라는 뜻이다.
                break
            list_RefSeq_SingleEntry[cnt].Get_mRNASeq(full_seq_chr)#해당 RefSeq class에 mRNA 시퀀스를 넣어준다.
            if(check_valid_mRNA(list_RefSeq_SingleEntry[cnt])):
                templist.append(list_RefSeq_SingleEntry[cnt])
            cnt+=1  
        #End of while body for cnt 
    #End of for body for chr_num
    return templist
######################################################################## End of make_list_valid

def make_list_raw_NM(file):# 전체 RefSeq을 담는 리스트와 NM_만 골라낸 RefSeq을 담는 리스트를 반환하는 함수
    templist=[]
    templist_raw=[]
    for sLine in file.readlines():
        temp_RefSeq=RefSeq()
        temp_RefSeq.parsing(sLine)
        templist_raw.append(temp_RefSeq)
        if(temp_RefSeq.RefSeqID[:3]=='NM_' and (temp_RefSeq.ChrID in chr_list)):
            templist.append(temp_RefSeq)
    # End of for body for sLine
    return templist_raw, templist
######################################################################## End of make_list_NM

def fancy_print(*args):
    length=len(args)
    templist=[str(i+1)+"번 답:\t" for i in range(length)]
    for elements in zip(templist, args):
        for element in elements:
            print(element, end="")
        print()
    #End of for body for element
######################################################################## End of fancy_print

def main():
    global chr_list
    start=time.time()
    file=open("refFlat.txt", 'r')
    list_RefSeq_raw, list_RefSeq_NM=make_list_raw_NM(file)
    dict_ID_entry=make_dict_entry(list_RefSeq_NM)# entry가 1개 이상인지를 확인하기 위한 딕셔너리
    file.close()
    outfile=open("result.txt", 'w')
    list_RefSeq_SingleEntry=delete_multientry(dict_ID_entry, list_RefSeq_NM)
    list_mRNA_valid=make_list_valid(list_RefSeq_SingleEntry) # 4번 answer에 대한 list
    list_mRNA_valid.sort(key=lambda x:x.NUM_RefSeqID)#각 RefSeqID의 뒤에 숫자에 대해서 sorting
    dict_ID_isoform=make_dict_for_isoform(list_mRNA_valid)#isoform인 친구들은 하나만 나타나도록 딕셔너리 만든다
    list_mRNA_final=leave_representitive_isoform(dict_ID_isoform,list_mRNA_valid) # 5번 answer에 대한 list
    print_outfile(list_mRNA_final,outfile) # 엑셀에 쓸 결과물을 txt 파일로 출력한다.
    outfile.close()
    fancy_print(len(list_RefSeq_raw), len(list_RefSeq_NM), len(list_RefSeq_SingleEntry), len(list_mRNA_valid), len(list_mRNA_final))
    print(list_mRNA_final[0].Length_3UTR,list_mRNA_final[0].dict_motif)
    print("총 걸린 시간은 "+str(time.time()-start)+"초입니다.")
    outfile.close()
######################################################################## End of main

main()