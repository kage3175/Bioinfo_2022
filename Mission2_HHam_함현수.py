#NC_045512v2.fa

import time
from itertools import product

def counting_oligonu(species,*args):  # 딕셔너리와 검사할 mer들을 입력받음. 각 mer에 대한 결과를 리스트로 만들어 반환.
    #print(chr_num,end=' ', flush=True)
    list_dict=[]
    file_tmp = open("../files_bioinfo2022/"+species + ".fa", "r")
    trash = file_tmp.readline()
    full_chr_tmp = file_tmp.read()
    full_chr_tmp = file_processing(full_chr_tmp)#대문자로 바꾸고 개행문자 삭제
    for MULTI in args:
        dict_for_oligonu = {}
        full_chr=full_chr_tmp
        multi_bases=full_chr[:MULTI]#full_chr의 가장 앞 2개만 가져오기
        full_chr=full_chr[MULTI:]#이미 읽은 부분 잘라내기
        dict_for_oligonu[multi_bases]=1
        for ch in full_chr:
            multi_bases=multi_bases[1:]+ch
            try:
                dict_for_oligonu[multi_bases] += 1
            except KeyError:
                dict_for_oligonu[multi_bases]=1
        list_dict.append(dict_for_oligonu)
        # End of for body for ch
    #End of for body for MULTI
    file_tmp.close()
    return list_dict
######################################################################## End of counting_oligonu

def file_processing(sFile):#문자열 file을 건네받고 문자열 내부의 문자들은 전부 대문자로 바꾸고, 개행문자를 삭제함
    newsFile=sFile.replace('\n', '')
    newsFile=newsFile.upper()
    return newsFile

def ratio(dict_of_x, base):#base에 들어있는 키에 해당하는 딕셔너리 값들의 비율을 계산하여 딕셔너리 형태로 반환하는 함수
    ratio_of_ACGT = {}
    sum = 0
    for sKey in base:
        try:
            sum += dict_of_x[sKey]
        except KeyError: #혹시나 특정 base가 아예 없을 수도 있으므로
            continue
    # End of for sKey(sum 구하기)
    for sKey in base:
        try:
            ratio_of_ACGT[sKey] = float(dict_of_x[sKey]) / sum
        except KeyError: #특정 base가 아예 없을 경우는 비율에 0.0 대입
            ratio_of_ACGT[sKey] = 0.0
    # End of for sKey(비율 구하기)
    return ratio_of_ACGT
######################################################################## End of ratio

def ideal_prob(ratio_ACGT,combination_multi): #이상적인 비율을 계산하는 함수
    result={}
    for combination in combination_multi:
        temp_ratio=1
        for ch in combination:
            temp_ratio*=ratio_ACGT[ch]
        result[combination]=temp_ratio
    return result
######################################################################## End of ideal_prob

def print_dict_sorted(dict_to_print,list_of_key): #딕셔너리를 키 값의 순서에 맞게 출력해주는 함수. 딕셔너리 순서는 AA ~ TT가 아니라 랜덤 키부터 시작할 수 있기 때문에, 깔끔한 출력을 위해 따로 함수를 지정하였다.
    for sKey in list_of_key:
        try:
            print(sKey,'\t',dict_to_print[sKey])
        except KeyError:
            print(sKey,'\t', 0.0)
            continue
    # End of for body for sKey
######################################################################## End of print_dict_sorted

def fancy_print(name_species,list_seperation,list_dictionary,key_sort): #깔끔한 출력을 위한 함수.
    print("\n>>", name_species,"\n")
    cnt=0
    for seperation in list_seperation:
        print("-------// ", seperation, "//-------")
        print_dict_sorted(list_dictionary[cnt],key_sort[cnt])
        print('SUM:',"\t",dict_sum(list_dictionary[cnt]))
        print('')
        cnt+=1
'''
>>종 이름                                   #name_species

-------// mononucleotide의 개수 //-------   #list_seperation에 여기 들어갈 문장이 적혀있다
A   nnn                                     #list_dictionary, key_sort에 여기 들어갈 요소가 적혀있다.
C   nnn

이런식으로 출력된다.

'''
######################################################################## End of fancy_print

def dict_sum(dict):
    sKeylist=list(dict.keys())
    sum=0
    for sKey in sKeylist:
        sum+=dict[sKey]
    return sum
######################################################################## End of dict_sum

def make_product(base, num_product):
    temparr=[]
    production=product(base,repeat=num_product)
    for tup in production:
        temp_str=''
        for i in range(num_product):
            temp_str=temp_str + tup[i]
        temparr.append(temp_str)
    return temparr

def main():
    start = time.time()
    base = ['A', 'C', 'G', 'T']
    list_seperation = ["mononucleotide의 개수", "dimer nucleotide의 개수", "tetramer nucleotide의 개수",
                       "mononucleotide의 비율", "dimer nucleotide의 비율","tetramer nucleotide의 비율",
                       "이상적인 dimer nucleotide의 비율", "이상적인 tetramer nucleotide의 비율"]  ##총 8개
    combination_dimer = make_product(base,2)
    combination_tetramer = make_product(base,4)
    temp_list_dict=counting_oligonu("NC_045512v2",1,2,4)
    result_mono=temp_list_dict[0]
    result_dimer=temp_list_dict[1]
    result_tetramer=temp_list_dict[2]
    ratio_of_ACGT=ratio(result_mono,base)# 각 mononucleotide의 비율, 딕셔너리
    ratio_of_dimer=ratio(result_dimer,combination_dimer)# 각 dinucleotide의 비율, 딕셔너리
    ratio_of_tetramer=ratio(result_tetramer,combination_tetramer)# 각 tetranucleotide의 비율, 딕셔너리
    ideal_ratio_dimer=ideal_prob(ratio_of_ACGT,combination_dimer) # dimer 이론상 비율 저장
    ideal_ratio_tetramer = ideal_prob(ratio_of_ACGT, combination_tetramer)  # tetramer 이론상 비율 저장
    list_dictionary=[result_mono,result_dimer,result_tetramer,ratio_of_ACGT,ratio_of_dimer,ratio_of_tetramer,ideal_ratio_dimer,ideal_ratio_tetramer]# 8개
    key_sort=[base,combination_dimer,combination_tetramer,base,combination_dimer,combination_tetramer,combination_dimer,combination_tetramer]# 8개
    fancy_print("SARS-CoV-2",list_seperation,list_dictionary,key_sort)

    print("걸린 시간은", end=' ')  # 이하 3줄은 시간출력용
    print(time.time() - start, end='')
    print("초입니다.\n")

main()




