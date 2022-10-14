import time, sys
from joblib import Parallel, delayed
from itertools import product

NUM_CORE=4
MULTI=3

def counting(species,chr_num):  # 딕셔너리와 파일 번호를 입력받음
    dict_for_ACGTN = {'A': 0, 'C': 0, 'G': 0, 'T': 0}  # 카운팅을 위한 딕셔너리
    file = open("./"+species+"Chrfiles/chr" + str(chr_num) + ".fa", "r")
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
# End of counting

def counting_dinu(species,chr_num):  # 딕셔너리와 파일 번호를 입력받음
    dict_for_dinu = {}
    file = open("./"+species+"Chrfiles/chr" + str(chr_num) + ".fa", "r")
    trash = file.readline()
    full_chr = file.read()
    full_chr = file_processing(full_chr)
    two_bases=full_chr[:MULTI]#full_chr의 가장 앞 MULTI개만 가져오기
    full_chr=full_chr[MULTI:]#이미 읽은 부분 잘라내기
    dict_for_dinu[two_bases]=1
    for ch in full_chr:
        two_bases=two_bases[1:]+ch
        try:
            dict_for_dinu[two_bases] += 1
        except KeyError:
            dict_for_dinu[two_bases]=1
    # End of for body for ch
    file.close()
    return dict_for_dinu
# End of counting_dinu

def file_processing(file):
    newfile=file.replace('\n', '')
    newfile=newfile.upper()
    return newfile

def ratio(dict_of_x, base):
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
    # End of ratio

def ideal_prob(ratio_ACGT,combination_multi):
    result={}
    for combination in combination_multi:
        temp_ratio=1
        for ch in combination:
            temp_ratio*=ratio_ACGT[ch]
        result[combination]=temp_ratio
    return result
# End of ideal_prob

def print_dict(dict_to_print):
    list_sKey = list(dict_to_print.keys())
    for sKey in list_sKey:
        print(sKey,"\t",dict_to_print[sKey])
#End of print_dict

def merging(list_of_dict):
    result={}
    for dict in list_of_dict:
        list_sKey=list(dict.keys())
        for sKey in list_sKey:
            try:
                result[sKey]+=dict[sKey]
            except KeyError:
                result[sKey]=dict[sKey]
    return result

def return_result_counting(list_of_chr, name_species):
    result=Parallel(n_jobs=NUM_CORE)(delayed(counting)(name_species,i) for i in list_of_chr)
    return result
#End of return_result_counting

def return_result_counting_dinu(list_of_chr, name_species):
    result=Parallel(n_jobs=NUM_CORE)(delayed(counting_dinu)(name_species,i) for i in list_of_chr)
    print(result)
    return result
#End of return_result_counting_dinu

def remove_N(dict):
    list_sKey=list(dict.keys())
    for sKey in list_sKey:
        for ch in sKey:
            if ch=='N':
                del dict[sKey]
                break
    return dict

def print_dict_sorted(dict_to_print,list_of_key):
    for sKey in list_of_key:
        try:
            print(sKey,"\t",dict_to_print[sKey])
        except KeyError:
            continue
# End of print_dict_sorted

def main():
    start = time.time()
    base = ['A', 'C', 'G', 'T']
    combination_multi=[]
    production=product(base,repeat=MULTI)
    for tup in production:
        temp_str=''
        for i in range(MULTI):
            temp_str=temp_str+tup[i]
        combination_multi.append(temp_str)
    species=["hg38","galGal6","dm6","ce11"]
    list_chr_hg = [i for i in range(1, 23)]
    list_chr_hg.append('X')
    list_chr_hg.append('Y')
    list_chr_galGal = [i for i in range(1, 29)]
    for i in range(30, 34):
        list_chr_galGal.append(i)
    list_chr_galGal.append('W')
    list_chr_galGal.append('Z')
    list_chr_dm = ['2R', '2L', '3R', '3L', 4, 'X', 'Y']
    list_chr_ce = ['I', 'II', 'III', 'IV', 'V', 'X']
    dict_list_chr = {"hg38": list_chr_hg, "galGal6": list_chr_galGal, "dm6": list_chr_dm, "ce11": list_chr_ce}
    for name_species in species:
        print("\nresult of",name_species)
        results_counting=return_result_counting(dict_list_chr[name_species],name_species)
        results_dinu=return_result_counting_dinu(dict_list_chr[name_species],name_species)
        result_mono=merging(results_counting)
        result_di=merging(results_dinu)
        ratio_of_ACTG = ratio(result_mono, base)  # 비율 구해서 ratio_of_ACTG 딕셔너리에 저장
        ratio_of_twobases=ratio(result_di,combination_multi)
        ideal_ratio=ideal_prob(ratio_of_ACTG,combination_multi)
        result_di=remove_N(result_di)################### 결과에서 N 들어간거 삭제. N도 표기하고 싶으면 이 줄 없애면 됨
        print(">>", name_species,"\n")

        print_dict_sorted(result_mono,base)
        print("/////////////////////// 이상 mononucleotide 개수")
        print_dict_sorted(result_di,combination_multi)
        print("/////////////////////// 이상 dinucleotide 개수")
        print_dict_sorted(ratio_of_ACTG,base)
        print("/////////////////////// 이상 mononucleotide 비율")
        print_dict_sorted(ratio_of_twobases,combination_multi)
        print("/////////////////////// 이상 dinucleotide 비율")
        print_dict_sorted(ideal_ratio,combination_multi)
    print("걸린 시간은", end=' ')  # 이하 3줄은 시간출력용
    print(time.time() - start, end='')
    print("초입니다.\n")
    # End of main

main()