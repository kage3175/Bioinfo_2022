#include <stdio.h>
#include <string.h>
#include <stdlib.h>

char* str_splice(char *buf, char*s, int start,int size){
	if(strlen(s)>start){
		s+=start;
		while(size-->0&&*s!='\0'){
			*(buf++)=*(s++);
		}
		*(buf)='\0';
	}
	return buf;
}

int main(void){
	int num=0;
	char test11[100];
	char test22[100]="chr12";
	printf("%d\n", strncmp(test11,test22,-1));
	if(test22[4]!='\0'){
		str_splice(test11,test22,3,2);
		printf("%s\n", test11);
		num=atoi(test11);
	}
	printf("%d\n",num);
	return 0;
}
