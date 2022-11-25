import time
import datetime
from itertools import product
from scipy import stats
import sys


chr_list=[str(i) for i in range(1,23)]
chr_list.append('X')
chr_list.append('Y')#Chr list를 담기위한 것. 사람의 유전자 이름은 절대 바뀔 일이 없으므로 전역변수로 설정하였다.
chr_list.sort()

BASE_COMPLE={'A':'T','C':'G','G':'C','T':'A'} #상보적인 염기쌍을 전역변수로서 미리 잡아둔다. A -T, C- G 간 결합을 딕셔너리로 표현한 것
STOP_CODON=['TAA','TGA', 'TAG']
BASE=['A','C','G','T']
CUTOFF=-0.5
PSUEDO_COUNT=0.00000000001

class RefSeq:
    global BASE_COMPLE
    global CUTOFF
    def __init__(self):
        self.RefSeqID='NULL'
        self.Gene_Symbol='NULL'
        self.ChrID=0
        self.Strand='NULL'
        self.Txn_start=0
        self.Txn_end=0
        self.Coding_start=0
        self.Coding_end=0
        self.Exon_num=0
        self.Exon_starts=[]
        self.Exon_ends=[]
        self.NUM_RefSeqID=0
        self.mRNASeq=''
        self.ORF_start=-1
        self.ORF_end=-1    #ORF_end는 3'UTR start와 일치한다. end 자체가 exclusive이기 때문
        self.Length_3UTR=0
        self.Length_mRNA=0
        self.Length_ORF=0
        self.dict_motif_3UTR={}
        self.dict_motif_ORF={}
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
    def Put_mRNASeq(self,file_chr): #chromosome 전체 sequence를 넣어주면 mRNA 서열을 저장하는 메소드
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
            for i in range(len(mRNA_seq)-1,-1,-1):
                self.mRNASeq=self.mRNASeq+BASE_COMPLE[mRNA_seq[i]]
            self.ORF_end, self.ORF_start=self.Length_mRNA-self.ORF_start, self.Length_mRNA-self.ORF_end # ORF end와 ORF start도 - strand에 맞게 변환시켜준다.
        else:
            self.mRNASeq=mRNA_seq
        self.Length_ORF=self.ORF_end-self.ORF_start#ORF 길이 저장
        self.Length_3UTR=self.Length_mRNA-self.ORF_end
    #End of Put_mRNASeq
    def Cal_motif_3UTR(self):
        if(self.Length_3UTR>=7):
            for i in range(self.ORF_end, self.Length_mRNA-6):
                self.dict_motif_3UTR[self.mRNASeq[i:i+7]]=1
    def Cal_motif_ORF(self):
        if(self.Length_ORF>=7):
            for i in range(self.ORF_start, self.ORF_end-6):
                self.dict_motif_ORF[self.mRNASeq[i:i+7]]=1
    def Get_Length_ORF(self): ### 이하 함수들은 ORF_start 등의 인스턴스 변수들을 반환해주는 함수들이다.
        return self.Length_ORF
    def Get_mRNASeq(self):
        return self.mRNASeq
    def Get_ORF_start(self):
        return self.ORF_start
    def Get_ORF_end(self):
        return self.ORF_end
    def Get_dict_motif_3UTR(self):
        return self.dict_motif_3UTR
    def Get_dict_motif_ORF(self):
        return self.dict_motif_ORF
    def Get_RefSeqID(self):
        return self.RefSeqID
    def Get_NUM_RefSeqID(self):
        return self.NUM_RefSeqID
    def Get_Gene_Symbol(self):
        return self.Gene_Symbol
    def Get_ChrID(self):
        return self.ChrID
    def Get_Length_mRNA(self):
        return self.Length_mRNA
    #End of Get_mRNASeq
######################################################################## End of RefSeq

class RefSeq_Fisher:
    global PSUEDO_COUNT
    def __init__(self):
        self.n1=0
        self.n2=0
        self.n3=0
        self.n4=0
        self.motif='AAAAAAA'
        self.pvalue=1
        self.Relative_Risk=0
    def Setmotif(self, motif_):
        self.motif=motif_
    def n1_plus(self):
        self.n1+=1
    def n2_plus(self):
        self.n2+=1
    def Cal_n3_n4(self,total_down, total_notdown):
        self.n3=total_down-self.n1
        self.n4=total_notdown-self.n2
    def Cal_RelativeRisk(self):
        A=(self.n1+PSUEDO_COUNT)/(self.n1+self.n2+PSUEDO_COUNT)
        B=(self.n3+PSUEDO_COUNT)/(self.n3+self.n4+PSUEDO_COUNT)
        self.Relative_Risk=A/B
    def Cal_pvalue(self):
        list_n=[[self.n1,self.n2], [self.n3,self.n4]]
        try:
            temp, self.pvalue=stats.fisher_exact(list_n)
        except ValueError:
            print(self.motif, list_n)
            sys.exit()
    def Init_ns(self):
        self.n1=0
        self.n2=0
        self.n3=0
        self.n4=0
    def Get_motif(self): #### 이하 함수들은 motif, pvalue 등을 반환해주는 함수. 클래스 내부 요소에 직접 접근하지 않도록 한다.
        return self.motif
    def Get_pvalue(self):
        return self.pvalue
    def Get_n1(self):
        return self.n1
    def Get_n2(self):
        return self.n2
    def Get_n3(self):
        return self.n3
    def Get_n4(self):
        return self.n4
    def Get_Relative_Risk(self):
        return self.Relative_Risk
######################################################################## End of RefSeq_Fisher

def check_valid_mRNA(refseq): # mRNA의 ORF 길이가 3의 배수인지, start codon이나 stop codon이 시작과 끝에 있는지, 중간에 멈춰버리지 않는 지 등을 검사
    global STOP_CODON
    temp_mRNASeq=refseq.Get_mRNASeq()
    ORF_start=refseq.Get_ORF_start()                                                  
    ORF_end=refseq.Get_ORF_end()
    Length_ORF=refseq.Get_Length_ORF()
    if(Length_ORF%3!=0): # ORF의 길이가 3의 배수인지
        return False
    elif(temp_mRNASeq[ORF_start:ORF_start+3]!='ATG'): # ORF의 첫 3글자가 ATG인지
        return False
    elif(temp_mRNASeq[ORF_end-3:ORF_end] not in STOP_CODON): # ORF의 마지막 3글자가 종결코돈인지
        return False
    else:
        for i in range(Length_ORF//3 -1):#중간에 STOP CODON이 존재하는지 검사. 끝은 어차피 stop codon이어야 하므로 검사하지 않는다.
            if(temp_mRNASeq[ORF_start+i*3:ORF_start+i*3+3] in STOP_CODON):
                return False
        #End of for body for i
        return True# 다 통과하면 True
######################################################################## End of check_valid_mRNA

def count_Fisher_variables(dataset, dict_RefSeq_Fisher, dict_RefSeq_Fisher_ORF, dict_RefSeq):
    total_down=0
    total_notdown=0
    for data in dataset:
        if data=="": #Mission5_Dataset의 맨 마지막 줄은 빈 문자열이므로, 마지막 줄에 다다르면 for문을 종료한다.
            break
        templist=list(data.split("\t")) #templist[0]는 Gene_Symbol, templist[1]은 log2(fold change) 값이다.
        templist[1]=float(templist[1])
        templist[0]=templist[0].upper()
        if templist[1]<CUTOFF: # -0.5보다 fold change 값이 낮을 때. 그 gene의 3'UTR에 있는 모티프에 해당하는 n1 값을 증가시켜야 한다.
            try:
                sMotiflist_3UTR=list(dict_RefSeq[templist[0]].Get_dict_motif_3UTR().keys())
                sMotiflist_ORF=list(dict_RefSeq[templist[0]].Get_dict_motif_ORF().keys())
                total_down+=1
                for motif in sMotiflist_3UTR:
                    dict_RefSeq_Fisher[motif].n1_plus() #해당 모티프를 키로 하는 클래스의 n1 값을 1 증가시킨다. n3, n4는 계산된 n1, n2 값을 기반으로 down, notdown gene의 수에서 n1, n2를 빼서 얻을 수 있다.
                for motif in sMotiflist_ORF:
                    dict_RefSeq_Fisher_ORF[motif].n1_plus()
            except KeyError: #dataset에 들어 있는 Gene_Symbol 중에, refflat.txt나 Mission4의 Ans5에 들어있지 않은 Gene_Symbol이 있을 수도 있으므로 안전하게 try except를 사용하였다.
                #print(templist[0])
                continue
        else:
            try:
                sMotiflist_3UTR=list(dict_RefSeq[templist[0]].Get_dict_motif_3UTR().keys())
                sMotiflist_ORF=list(dict_RefSeq[templist[0]].Get_dict_motif_ORF().keys())
                total_notdown+=1
                for motif in sMotiflist_3UTR:
                    dict_RefSeq_Fisher[motif].n2_plus()
                    '''if(motif=='AAAATAA') :
                        print(total_notdown, dict_RefSeq_Fisher[motif].Get_n2(), end=" || ")'''
                for motif in sMotiflist_ORF:
                    dict_RefSeq_Fisher_ORF[motif].n2_plus()
            except KeyError:
                continue
    return total_down, total_notdown
######################################################################## End of count_Fisher_variables

def delete_multientry(list_RefSeq_NM): # entry가 여러 개인 경우를 제외하고 만든 (클래스)리스트를 반환
    dict_check=make_dict_entry(list_RefSeq_NM)# entry가 1개 이상인지를 확인하기 위한 딕셔너리
    templist=[]
    for refseq in list_RefSeq_NM:
        if dict_check[refseq.RefSeqID]==1:
            templist.append(refseq)
    # End of for body for refseq
    return templist
######################################################################## End of delete_multientry

def fancy_print(*args):
    length=len(args)
    templist=[str(i+1)+"번 답:\t" for i in range(length)]
    for elements in zip(templist, args):
        for element in elements:
            print(element, end="")
        print()
    #End of for body for element
######################################################################## End of fancy_print

def file_processing(sFile):#문자열 file을 건네받고 문자열 내부의 문자들은 전부 대문자로 바꾸고, 개행문자를 삭제함
    newsFile=sFile.replace('\n', '')
    newsFile=newsFile.upper()
    return newsFile
######################################################################## End of file_processing

def leave_representitive_isoform(dict_check, list_RefSeq_NM): # entry가 여러 개인 isoform 중 RefSeqID 숫자가 가장 작은 하나만 남겨서 리스트 반환
    templist=[]
    list_RefSeq_NM.sort(key=lambda x:x.Get_NUM_RefSeqID()) # NM 숫자가 낮은 애들을 representitive으로 남길 거기 때문에, NUM_RefSeqID를 키로 sorting 해준다.
    for refseq in list_RefSeq_NM:
        Gene_Symbol=refseq.Get_Gene_Symbol()
        if dict_check[Gene_Symbol]==1: #만약 해당 Gene_Symbol에 해당하는 refseq class가 이미 리스트에 추가되어 있다면, 아래 코드에 의해 딕셔너리 값이 0이 되므로, 같은 Gene_Symbol을 가지는 refseq class는 2번 이상 리스트에 추가될 수 없다.
            refseq.Cal_motif_3UTR()
            refseq.Cal_motif_ORF()
            templist.append(refseq)
            dict_check[Gene_Symbol]=0
    # End of for body for refseq
    return templist
######################################################################## End of leave_representitive_isoform

def make_dict_entry(list_RefSeq_NM):# 각 entry가 1개 이상인지를 확인하기 위한 딕셔너리를 만들어서 반환
    dict={}
    for refseq in list_RefSeq_NM:
        RefSeqID=refseq.Get_RefSeqID()
        try:
            dict[RefSeqID]+=1
        except KeyError:
            dict[RefSeqID]=1
    # End of for body for refseq
    return dict
######################################################################## End of make_dict_entry

def make_dict_Fisher(production_7mer, class_):
    dict_RefSeq_Fisher={}
    dict_RefSeq_Fisher_ORF={}
    for production in production_7mer:#각 모티프를 키로 하는 refseq_fisher 딕셔너리 초기화
        refseq_F=class_()
        refseq_F_ORF=class_()
        refseq_F.Setmotif(production)
        refseq_F_ORF.Setmotif(production)
        dict_RefSeq_Fisher[production]=refseq_F
        dict_RefSeq_Fisher_ORF[production]=refseq_F_ORF
    return dict_RefSeq_Fisher, dict_RefSeq_Fisher_ORF
######################################################################## End of make_dict_Fisher

def make_dict_isoform(list_RefSeq_validORF):# 다른 entry여도 Gene_Symbol이 여러 개인 경우를 추려내기 위한 딕션너리를 반환한다. 딕셔너리는 Gene_Symbol을 키로, 개수 상관없이 무조건 1로 저장된다.
    dict={}
    for refseq in list_RefSeq_validORF:
        Gene_Symbol=refseq.Get_Gene_Symbol()
        dict[Gene_Symbol]=1
    # End of for body for refseq
    return dict
######################################################################## End of make_dict_isoform

def make_list_raw_NM(file):# 전체 RefSeq을 담는 리스트와 NM_만 골라낸 RefSeq을 담는 리스트를 반환하는 함수
    templist=[]
    templist_raw=[]
    for sLine in file.readlines():
        temp_RefSeq=RefSeq()
        temp_RefSeq.parsing(sLine)
        templist_raw.append(temp_RefSeq)
        RefSeqID=temp_RefSeq.Get_RefSeqID()
        ChrID=temp_RefSeq.Get_ChrID()
        if(RefSeqID[:3]=='NM_' and (ChrID in chr_list)):
            templist.append(temp_RefSeq)
    # End of for body for sLine
    return templist_raw, templist
######################################################################## End of make_list_NM

def make_list_valid(list_RefSeq_SingleEntry):
    global chr_list
    templist=[]
    list_RefSeq_SingleEntry.sort(key=lambda x:x.Get_ChrID())#ChrID 순서대로 먼저 정렬. 이러지 않으면 for문을 5만번씩 24번 돌려야해서 비효율 적이다. 정렬하면 for문 5만번 1번만 돌리면 된다.
    length_list=len(list_RefSeq_SingleEntry)
    cnt=0
    for chr_num in chr_list: # 1번부터 Y까지의 크로모좀 전체 시퀀스를 해당 크로모좀 번호(string 형태)를 키로 하는 딕셔너리에 저장
        chr_file=open("../files_bioinfo2022/hg38ChrFiles/chr"+chr_num+".fa","r")############################## 제출할 때 바꿔야함
        trash=chr_file.readline()
        full_seq_chr=chr_file.read()
        full_seq_chr=file_processing(full_seq_chr)
        chr_file.close()
        while cnt<length_list: # RefSeq이 담긴 리스트의 길이보다 더 검사할 필요도 없고, 하면 오류가 나니까 끊어준다.
            ChrID=list_RefSeq_SingleEntry[cnt].Get_ChrID()
            if(ChrID!=chr_num):# ChrID가 현재 Chr과 다르면. 이미 ChrID 순서대로 sorting 했으므로, 다음 Chr을 검사해야할 시점이라는 뜻이다.
                break
            list_RefSeq_SingleEntry[cnt].Put_mRNASeq(full_seq_chr)#해당 RefSeq class에 mRNA 시퀀스를 넣어준다.
            if(check_valid_mRNA(list_RefSeq_SingleEntry[cnt])):
                templist.append(list_RefSeq_SingleEntry[cnt])
            cnt+=1  
        #End of while body for cnt 
    #End of for body for chr_num
    return templist
######################################################################## End of make_list_valid

def make_product(base, num_product):
    temparr=[]
    production=product(base,repeat=num_product)
    for tup in production:
        temp_str=''
        for i in range(num_product):
            temp_str=temp_str + tup[i]
        #End of for body for i
        temparr.append(temp_str)
    #End of for body for tup
    return temparr
######################################################################## End of make_product

def print_outfile(list_RefSeq, outfile):#Excel에 적을 답을 txt 파일에 출력해주는 함수
    for refseq in list_RefSeq:
        print(refseq.Get_RefSeqID()+'\t'+refseq.Get_Gene_Symbol()+'\t'+str(refseq.Get_ORF_start())+'\t'+str(refseq.Get_Length_ORF())+'\t'+str(refseq.Get_Length_mRNA()-refseq.Get_ORF_end()) ,file=outfile)
    #End for body for refseq
######################################################################## End of print_outfile

def print_result_Mission6(list_RefSeq_Fisher, bonferroni):
    cnt=0
    i=0
    while cnt<5: #상위 5개만 출력한다
        Relative_Risk=list_RefSeq_Fisher[i].Get_Relative_Risk()
        if(Relative_Risk>1):
            print(list_RefSeq_Fisher[i].Get_motif()+'\t'+str(list_RefSeq_Fisher[i].Get_pvalue()*bonferroni)+'\t'+str(list_RefSeq_Fisher[i].Get_n1())+'\t'+str(list_RefSeq_Fisher[i].Get_n2())+'\t'+str(list_RefSeq_Fisher[i].Get_n3())+'\t'+str(list_RefSeq_Fisher[i].Get_n4())+'\t'+str(list_RefSeq_Fisher[i].Get_Relative_Risk()))
            cnt+=1
        i+=1
######################################################################## End of print_result_Mission6

def print_time(str_print,start, end):
    print(str_print+str(end-start)+"초 입니다.")
######################################################################## End of fancy_print


def read_mi_Seqs(file):
    flag=0
    list_return=[]
    templist=[1,1]
    tempstr='a'
    for sLine in file.readlines():
        if flag==1:
            flag=0
            list_return.append([tempstr, sLine.replace('\n', '')])
        elif(sLine[:4]=='>hsa'):
            tempstr=sLine.replace('\n', '')
            flag=1
        else :
            flag=0
    return list_return
            
            
            



def main():
    '''global chr_list
    global BASE
    global CUTOFF
    print(datetime.datetime.now())
    count_test_case_3UTR=count_test_case_ORF=0
    start=time.time()
    dict_RefSeq={}
    production_7mer = make_product(BASE,7)
    file=open("../files_bioinfo2022/refFlat.txt", 'r')
    list_RefSeq_raw, list_RefSeq_NM=make_list_raw_NM(file)
    file.close()
    
    list_RefSeq_SingleEntry=delete_multientry(list_RefSeq_NM)
    list_mRNA_valid=make_list_valid(list_RefSeq_SingleEntry) # 4번 answer에 대한 list
    dict_ID_isoform=make_dict_isoform(list_mRNA_valid)#isoform인 친구들은 하나만 나타나도록 딕셔너리 만든다
    list_mRNA_final=leave_representitive_isoform(dict_ID_isoform,list_mRNA_valid) # 5번 answer에 대한 list
    #fancy_print(len(list_RefSeq_raw), len(list_RefSeq_NM), len(list_RefSeq_SingleEntry), len(list_mRNA_valid), len(list_mRNA_final))
    print_time("Mission4(출력 생략)까지 걸린 시간은 ",start, time.time())'''
    file=open("../files_bioinfo2022/mature.fa",'r')
    list_miSeqs=read_mi_Seqs(file)
    file.close()
    print(list_miSeqs[0], list_miSeqs[1])
    
    '''for refseq in list_mRNA_final: #list_mRNA_final에 있는 refseq들을 다 Gene_Symbol을 키로 하는 딕셔너리로 옮겨줌
        Gene_Symbol=refseq.Get_Gene_Symbol()
        dict_RefSeq[Gene_Symbol]=refseq
    dict_RefSeq_Fisher, dict_RefSeq_Fisher_ORF=make_dict_Fisher(production_7mer, RefSeq_Fisher)
    for i in range(3):
        count_test_case_ORF=count_test_case_3UTR=0
        ist_RefSeq_Fisher_3UTR=[]
        list_RefSeq_Fisher_ORF=[]
        for motif in production_7mer:
            dict_RefSeq_Fisher[motif].Init_ns()
            dict_RefSeq_Fisher_ORF[motif].Init_ns()
        file=open("../files_bioinfo2022/Mission6_Dataset"+str(i+1)+"_2022.txt", "r")
        dataset=list(file.read().split("\n"))
        file.close()
        # End of for body for refseq
        total_down, total_notdown = count_Fisher_variables(dataset, dict_RefSeq_Fisher,dict_RefSeq_Fisher_ORF, dict_RefSeq)
        print(total_down, total_notdown)
        for motif in production_7mer:
            dict_RefSeq_Fisher[motif].Cal_n3_n4(total_down, total_notdown)
            dict_RefSeq_Fisher_ORF[motif].Cal_n3_n4(total_down, total_notdown)
            dict_RefSeq_Fisher[motif].Cal_RelativeRisk()
            dict_RefSeq_Fisher_ORF[motif].Cal_RelativeRisk()
            if(dict_RefSeq_Fisher[motif].Get_Relative_Risk()>1): # Relative Risk가 1 초과인 경우만 검사
                count_test_case_3UTR+=1
                dict_RefSeq_Fisher[motif].Cal_pvalue()
                ist_RefSeq_Fisher_3UTR.append(dict_RefSeq_Fisher[motif])
            if(dict_RefSeq_Fisher_ORF[motif].Get_Relative_Risk()>1): # Relative Risk가 1 초과인 경우만 검사
                count_test_case_ORF+=1
                dict_RefSeq_Fisher_ORF[motif].Cal_pvalue()
                list_RefSeq_Fisher_ORF.append(dict_RefSeq_Fisher_ORF[motif])
        # End of for body for motif
        ist_RefSeq_Fisher_3UTR.sort(key=lambda x:x.Get_pvalue()) #pvalue 값을 기준으로 리스트를 정렬한다.
        list_RefSeq_Fisher_ORF.sort(key=lambda x:x.Get_pvalue())
        print("Mission6_dataset"+str(i+1)+" 3UTR 결과")
        print_result_Mission6(ist_RefSeq_Fisher_3UTR, count_test_case_3UTR)
        print("Mission6_dataset"+str(i+1)+" ORF 결과")
        print_result_Mission6(list_RefSeq_Fisher_ORF, count_test_case_ORF)
    
    
    print_time("총 걸린 시간은 ",start, time.time())'''
######################################################################## End of main

main()