/* 
 * icc -o how2read how2read.c -L -lmyram -g -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL
 *
 * To Run,  ./how2read [THE FILENAME YOU WANT TO READ]  [STAR/SINK/HCELL/DM]
 * */

#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<sys/stat.h>
#include<math.h>
#include "ramses.h"
#include "params.h"


int main(int argc, char *argv[]){
	size_t i,j,k;
	FILE *fp = fopen(argv[1],"r");

	DmType *dm;
	StarType *star;
	SinkType *sink;
	HydroCellType *hcell;
	GasType *gas;
	struct stat st;
	size_t size;
	void *aa;

	stat(argv[1], &st);
	size = st.st_size;
	size_t dsize;

	if(strcmp(argv[2],"STAR")==0){
		dsize = sizeof(StarType);
		aa = (void*) star = (StarType*)malloc(sizeof(char)*size);
	}
	else if(strcmp(argv[2],"SINK")==0){
		dsize = sizeof(SinkType);
		aa = (void*)sink = (SinkType*)malloc(sizeof(char)*size);
	}
	else if(strcmp(argv[2],"HCELL")==0){
		dsize = sizeof(HydroCellType);
		aa = (void*)hcell = (HydroCellType*)malloc(sizeof(char)*size);
	}
	else if(strcmp(argv[2],"GAS")==0){
		dsize = sizeof(GasType);
		aa = (void*)gas = (GasType*)malloc(sizeof(char)*size);
	}
	else if(strcmp(argv[2],"DM")==0){
		dsize = sizeof(DmType);
		aa = (void*)dm = (DmType*)malloc(sizeof(char)*size);
	}
	else {
		fprintf(stderr,"Oooooooooops. Wrong size in file and type\n");
	}

	if(size%dsize !=0){
		printf("Error in file size,,,, %s\n", argv[1]);
		exit(99);
	}
	fread(aa, sizeof(char), size, fp);
	fclose(fp);

}
