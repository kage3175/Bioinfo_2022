import time, sys
from itertools import product

MULTI=2

def counting(species,chr_num):  # 딕셔너리와 파일 번호를 입력받음
    dict_for_ACGTN = {'A': 0, 'C': 0, 'G': 0, 'T': 0,'N':0}  # 카운팅을 위한 딕셔너리
    file = open("./"+species+"ChrFiles/chr" + str(chr_num) + ".fa", "r")
    trash = file.readline()#fa 파일의 첫줄은 필요 없으므로
    full_chr = file.read()
    full_chr = file_processing(full_chr)
    for ch in full_chr:
        try:
            dict_for_ACGTN[ch] += 1
        except KeyError:
            #sys.exit("A,C,G,T,N 이외의 문자가 있습니다!")
            continue
    # End of for body for ch
    file.close()
    return dict_for_ACGTN
######################################################################## End of counting

def counting_dinu(species,chr_num):  # 딕셔너리와 파일 번호를 입력받음
    dict_for_dinu = {}
    file = open("./"+species+"ChrFiles/chr" + str(chr_num) + ".fa", "r")
    trash = file.readline()
    full_chr = file.read()#대문자로 바꾸고 개행문자 삭제
    full_chr = file_processing(full_chr)
    multi_bases=full_chr[:MULTI]#full_chr의 가장 앞 2개만 가져오기
    full_chr=full_chr[MULTI:]#이미 읽은 부분 잘라내기
    dict_for_dinu[multi_bases]=1
    for ch in full_chr:
        multi_bases=multi_bases[1:]+ch
        try:
            dict_for_dinu[multi_bases] += 1
        except KeyError:
            dict_for_dinu[multi_bases]=1
    # End of for body for ch
    file.close()
    return dict_for_dinu
######################################################################## End of counting_dinu

def file_processing(file):#문자열 file을 건네받고 문자열 내부의 문자들은 전부 대문자로 바꾸고, 개행문자를 삭제함
    newfile=file.replace('\n', '')
    newfile=newfile.upper()
    return newfile
######################################################################## End of file_processing

def ratio(dict_of_x, base):#base에 들어있는 키에 해당하는 딕셔너리 값들의 비율을 계산하여 딕셔너리 형태로 반환하는 함수
    ratio_of_ACGT = {}
    sum = 0
    for sKey in base:
        try:
            sum += dict_of_x[sKey]
        except KeyError:
            continue
    # End of for sKey(sum 구하기)
    for sKey in base:
        try:
            ratio_of_ACGT[sKey] = float(dict_of_x[sKey]) / sum
        except KeyError:
            ratio_of_ACGT[sKey] = 0.0
    # End of for sKey(비율 구하기)
    return ratio_of_ACGT
######################################################################## End of ratio

def ideal_prob(ratio_ACGT,combination_multi):
    result={}
    for combination in combination_multi:
        temp_ratio=1
        for ch in combination:
            temp_ratio*=ratio_ACGT[ch]
        result[combination]=temp_ratio
    return result
######################################################################## End of ideal_prob

def print_dict(dict_to_print):
    list_sKey=list(dict_to_print.keys())
    for sKey in list_sKey:
        print(sKey, dict_to_print[sKey])
    #End of for body for sKey
######################################################################## End of print_dict

def merging(list_of_dict):#리스트에 저장되어 있는 모든 딕셔너리의 값을 통합. [{1:0, 2:1}, {1:1, 2:1, 3:1}] 이라면 {1:1,2:2,3:1}로 통합되어 return 된다.
    result={}
    for dict in list_of_dict:
        list_sKey=list(dict.keys())
        for sKey in list_sKey:
            try:
                result[sKey]+=dict[sKey]
            except KeyError:
                result[sKey]=dict[sKey]
        #End of for body for sKey
    #End of for body for dict
    return result
######################################################################## End of merging

def return_result_counting(list_of_chr, name_species):#검사해야하는 chromosome 번호(또는 기호)와 검사할 종의 문자열(hg38 등)을 전달받아서 counting을 실행. 결과로 나온 딕셔너리들은 result 리스트에 저장된다.
    result=[]
    for chr in list_of_chr:
        result.append(counting(name_species,chr))
    # End of for body for chr
    return result
######################################################################## End of return_result_counting

def return_result_counting_dinu(list_of_chr, name_species):#검사해야하는 chromosome 번호(또는 기호)와 검사할 종의 문자열(hg38 등)을 전달받아서 counting_dinu을 실행. 결과로 나온 딕셔너리들은 result 리스트에 저장된다.
    result=[]
    for chr in list_of_chr:
        result.append(counting_dinu(name_species,chr))
    # End of for body for chr
    return result
######################################################################## End of return_result_counting_dinu

def remove_N(dict):
    list_sKey=list(dict.keys())
    for sKey in list_sKey:
        for ch in sKey:
            if ch=='N':
                del dict[sKey]
                break
    # End of for body for sKey, ch in sKey
    return dict
######################################################################## End of remove_N

def print_dict_sorted(dict_to_print,list_of_key):
    for sKey in list_of_key:
        try:
            print(sKey,dict_to_print[sKey])
        except KeyError:
            print(sKey, 0.0)
            continue
    # End of for body for sKey
######################################################################## End of print_dict_sorted

def main():
    start = time.time()
    base = ['A', 'C', 'G', 'T']
    species = ["hg38", "galGal6", "dm6", "ce11"]
    combination_multi = []
    production = product(base, repeat=MULTI)
    for tup in production:
        temp_str = ''
        for i in range(MULTI):
            temp_str = temp_str + tup[i]
        combination_multi.append(temp_str)
    list_chr_hg = [i for i in range(1, 23)]#이하 11줄은 각 종의 chromosome 넘버 또는 기호를 리스트에 저장하는 코드이다.
    list_chr_hg.append('X')
    list_chr_hg.append('Y')
    list_chr_galGal=[i for i in range(1,29)]
    for i in range(30,34):
        list_chr_galGal.append(i)
    # End of for body for i
    list_chr_galGal.append('W')
    list_chr_galGal.append('Z')
    list_chr_dm = ['2R', '2L', '3R', '3L',4,'X','Y']
    list_chr_ce=['I', 'II','III','IV','V','X']
    dict_list_chr = {"hg38":list_chr_hg,"galGal6":list_chr_galGal,"dm6":list_chr_dm,"ce11":list_chr_ce}#검사할 종의 기호(hg38)과 각 종마다 검사할 chromosome 리스트를 매칭하는 딕셔너리
    for name_species in species:
        result_mononu=return_result_counting(dict_list_chr[name_species],name_species)#현재 리스트에 chromosome의 수 만큼 딕셔너리가 저장되어 있는 상태
        result_dinu=return_result_counting_dinu(dict_list_chr[name_species],name_species)#현재 리스트에 chromosome의 수 만큼 딕셔너리가 저장되어 있는 상태
        result_mono=merging(result_mononu)#리스트에 저장된 딕셔너리를 하나의 딕셔너리로 통합
        result_di=merging(result_dinu)#리스트에 저장된 딕셔너리를 하나의 딕셔너리로 통합
        ratio_of_ACTG = ratio(result_mono, base)  # mononucleotide 비율 구해서 ratio_of_ACTG 딕셔너리에 저장
        ratio_of_multibases = ratio(result_di, combination_multi)# dinucleotide 비율 구해서 ratio_of_multibases 딕셔너리에 저장
        ideal_ratio=ideal_prob(ratio_of_ACTG,combination_multi) # 이론상 비율 저장
        result_di=remove_N(result_di)################### 결과에서 N 들어간거 삭제. 이 줄 없애면 됨
        print("\n>>result of", name_species)
        print("\n// "+str(MULTI)+"-mer "+"nucleotide의 개수 //-----------------------")
        print_dict(result_mono)
        print("\n// "+str(MULTI)+"-mer "+"nucleotide의 개수 //-----------------------")
        print_dict_sorted(result_di,combination_multi)
        print("\n// "+str(MULTI)+"-mer "+"nucleotide의 비율 //-----------------------")
        print_dict(ratio_of_ACTG)
        print("\n// "+str(MULTI)+"-mer "+"nucleotide의 비율 //-----------------------")
        print_dict_sorted(ratio_of_multibases,combination_multi)
        print("\n// 이상적인"+str(MULTI)+"-mer "+" dinucleotide의 비율 //-----------------------")
        print_dict(ideal_ratio)
    # End of for body for name_species
    print("걸린 시간은", end=' ')  # 이하 3줄은 시간출력용
    print(time.time() - start, end='')
    print("초입니다.\n")
######################################################################## End of main

main()