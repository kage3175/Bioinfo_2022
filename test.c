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
	char temp[40]="002341";
	int x=0;
	x=atoi(temp);
	printf("%d", x);
	return 0;
}
