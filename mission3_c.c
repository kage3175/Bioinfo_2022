#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include <memory.h>

//int compare(const void*a, const void*b);



typedef struct _mRNA{
	char Gene_Symbol[20];
	char RefSeqID[20];
	char ChrID[7];
	char Strand[2];
	int Txn_start;
	int Txn_end;
	int coding_start;
	int coding_end;
	int Num_Exon;
	long long Exon_Starts[500];
	long long Exon_Ends[500];
	int Num_RefSeqID;
}mRNA;

void str_splice(char *buf, char*s, int start,int size){
	if(strlen(s)>start){
		s+=start;
		while(size-->0&&*s!='\0'){
			*(buf++)=*(s++);
		}
		*(buf)='\0';
	}
}

long long int* split_to_longlong(char* division, long long int results[], char temp_str[]){
	int i=0;
	char temp[5000];
	char* p;
	p=strtok(temp_str,division);
		while(p!=NULL){
			strcpy(temp,p);
			results[i]=atoi(temp);
			i++;
			p=strtok(NULL,",");
		}
	return results;
}

static int compare(const void *a,const void *b)
{
  mRNA *e1 = (mRNA*)a;
  mRNA *e2 = (mRNA*)b;

  if (e1->Num_RefSeqID > e2->Num_RefSeqID) return 1;
  else if (e1->Num_RefSeqID < e2->Num_RefSeqID) return -1;
  else return 0;
}

int main(void){
	int i=0;int j=0;int k=0;
	char line[10000];
	char temp_line[10000];
	char* p;
	char temp_str[11][5000];
	char temp[500];
	long long exon_starts[1000]={0};
	long long exon_ends[1000]={0};
	mRNA * temp_mRNA;
	temp_mRNA = (mRNA *)malloc(sizeof(mRNA)*100000);
	int total_num_NMgenes=0;
	char chr_id[7];
	int temp_num=0;
	bool flag=true;
	
	
	FILE* file=fopen("../files_bioinfo2022/refFlat.txt", "r");
	if(file==NULL){
		printf("파일 열기 실패\n");
		return 1; 
	}
	
	while(fgets(line, sizeof(line), file) != NULL){
		i=0;
		strcpy(temp_line,line);
		for(j=0;j<11;j++){
			k=0;
			for(k=0;k<strlen(temp_str[j]);k++){
				temp_str[j][k]='\0';
			}
		}
		//printf("1\n");
		p=strtok(temp_line,"\t");
		//printf("%s\n", line);
		while(p!=NULL){
			strcpy(temp_str[i],p);
			p=strtok(NULL,"\t");
			i++;
		}
		//printf("%s\n", temp_str[9]);
		//printf("%s\n", line);
		if(strncmp(temp_str[1],"NM_",3)==0){
			if (temp_str[2][4]!='/0') str_splice(chr_id,temp_str[2],0,5);
			else str_splice(chr_id,temp_str[2],0,4);
			//if(!(chr_id[3]=='1'||chr_id[3]=='2'||chr_id[3]=='3'||chr_id[3]=='4'||chr_id[3]=='5'||chr_id[3]=='6'||chr_id[3]=='7'||chr_id[3]=='8'||chr_id[3]=='9'||chr_id[3]=='X'||chr_id[3]=='Y')) continue;
			if(strcmp(temp_str[2],chr_id)==0){
				flag=true;
				str_splice(temp,temp_str[1],3,100);
				temp_num=atoi(temp);
				for(i=0;i<total_num_NMgenes;i++){
					if(temp_num==temp_mRNA[i].Num_RefSeqID){
						//printf("%s %s\n",temp_mRNA[i].RefSeqID, temp_str[1]);
						flag=false;
						break;
					}
				}
				if(flag==true){
					strcpy(temp_mRNA[total_num_NMgenes].Gene_Symbol,temp_str[0]);
					strcpy(temp_mRNA[total_num_NMgenes].RefSeqID,temp_str[1]);
					//strcat(temp_mRNA[total_num_NMgenes].RefSeqID,"\n");
					str_splice(temp,temp_str[1],3,100);
					temp_mRNA[total_num_NMgenes].Num_RefSeqID=atoi(temp);
					strcpy(temp_mRNA[total_num_NMgenes].ChrID,temp_str[2]);
					strcpy(temp_mRNA[total_num_NMgenes].Strand, temp_str[3]);
					temp_mRNA[total_num_NMgenes].Txn_start=atoi(temp_str[4]);
					temp_mRNA[total_num_NMgenes].Txn_end=atoi(temp_str[5]);
					temp_mRNA[total_num_NMgenes].coding_start=atoi(temp_str[6]);
					temp_mRNA[total_num_NMgenes].coding_end=atoi(temp_str[7]);
					temp_mRNA[total_num_NMgenes].Num_Exon=atoi(temp_str[8]);
					split_to_longlong(",",temp_mRNA[total_num_NMgenes].Exon_Starts,temp_str[9]);
					split_to_longlong(",",temp_mRNA[total_num_NMgenes].Exon_Ends, temp_str[10]);
					total_num_NMgenes++;
				}
			}
			else continue;
			//printf("%d\n", temp_mRNA[total_num_NMgenes-1].Num_RefSeqID);
		}
	}
	qsort(temp_mRNA, total_num_NMgenes, sizeof(mRNA), compare);
	fclose(file);
	FILE* outfile=fopen("../files_bioinfo2022/result_c.txt","w");
	printf("%d\n", total_num_NMgenes);
	for(i=0;i<total_num_NMgenes;i++) {
		fprintf(outfile,"%s\n", temp_mRNA[i].RefSeqID);
	}
	//printf("..%d..\n",exon_starts[0]);
	fclose(outfile);
	free(temp_mRNA);
	
	return 0;
}


