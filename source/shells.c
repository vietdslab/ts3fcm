#include<stdio.h>
#include<stdlib.h>

int main( int  argc, char **argv){
    int k,i,j,start = 0,stop = 10, starti = 0, stopi = 30;
	char *run,s[100],*str[50] = {"australian","dermatology","iris","wine","spambase","heart","balance-scale","wdbc","tae","waveform"};	
	if(argc >= 2){
		run = argv[1]; 
		if(argc >= 3){
			start = atoi(argv[2]);
			if(argc >= 4){
				stop = atoi(argv[3]);
				if(argc >= 5){
					starti = atoi(argv[4]);
					if(argc >= 6){
						stopi = atoi(argv[5]);						
					}
				}
			}
		}
	} 
				
	for(k = start;k < stop;k++){
		for(i = starti;i<=stopi;i+=5){
			sprintf(s,"%s.exe %s 20 %d",run,str[k],i);
			printf("running %s\n",s);
			system(s);	
		}		
	}
    return 0;
}
