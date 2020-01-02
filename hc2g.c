/* 
icc -o hc2g hc2g.c -L. -lmyram -g -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL -DWGROUPSIZE=4
 *
 * To Run,  ./how2read [THE FILENAME YOU WANT TO READ]  [STAR/SINK/HCELL/DM]
 * */

#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<sys/stat.h>
#include<math.h>
#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ((a)%(b))

#define Wait2Start(WGroupSize) do{}while(0)
#define Wait2Go(WGroupSize) do{}while(0)
#include "ramses.h"
#include "params.h"

#define MIN(a,b) ( (a)<(b) ? (a):(b) )

idtype id=0;

dptype Mpc = 3.0857E24;
dptype Msun = 2E33;
RamsesType simpar;

int hc2g(char *, char*);

int main(int argc, char *argv[]){
	int nfile,nstep;

	nstep=atoi(argv[1]);
	nfile=atoi(argv[2]);
	char infile[190],outfile[190];
	sprintf(infile,"HR5.%.5d.00000.info",nstep);
	FILE *fp = fopen(infile,"r");
	simpar = read_head(fp);
	fread(&simpar, sizeof(RamsesType), 1, fp);
	fclose(fp);
	int i;
	id = 0;
	for(i=375;i<384;i++)
	{
		sprintf(infile,"HR5.%.5d.HCELL.%.5d.dat",nstep,i);
		sprintf(outfile,"HR5.%.5d.GAS.%.5d.dat",nstep,i);
		hc2g(infile,outfile);
	}
}
int hc2g(char *infile, char *outfile){
	size_t i,j,k;

	HydroCellType *hcell;
	struct stat st;
	size_t size;
	void *aa;
	GasType *gas;
	int WGroupSize = 1000;

	WGroupSize = WGROUPSIZE;

	stat(infile, &st);
	size = st.st_size;
	size_t dsize;

	dsize = sizeof(HydroCellType);
	aa = (void*)hcell = (HydroCellType*)malloc(sizeof(char)*size);

	printf("Opening %s to read hcell data\n",infile);fflush(stdout);
	FILE *fp = fopen(infile,"r");
	if(size%dsize !=0){
		printf("Error in file size,,,, %s\n", infile);
		exit(99);
	}
	Wait2Start(WGroupSize);
	fread(aa, sizeof(char), size, fp);
	Wait2Go(WGroupSize);
	size_t np = size/dsize;
	gas = (GasType*)malloc(sizeof(GasType)*np);
	for(i=0;i<np;i++){
		gas[i].x = hcell[i].x;
		gas[i].y = hcell[i].y;
		gas[i].z = hcell[i].z;
		gas[i].cellsize = hcell[i].cellsize;
		gas[i].vx = hcell[i].vx;
		gas[i].vy = hcell[i].vy;
		gas[i].vz = hcell[i].vz;
		gas[i].den = hcell[i].den;
		gas[i].temp = hcell[i].temp;
		gas[i].metallicity = hcell[i].metallicity;
		gas[i].H = hcell[i].H;
		gas[i].O = hcell[i].O;
		gas[i].Fe = hcell[i].Fe;
		gas[i].ilevel = hcell[i].ilevel;
		gas[i].indx = id++;
		dptype volcell = pow(hcell[i].cellsize*Mpc * simpar.aexp/simpar.H0*100L,3.L);
		gas[i].mass = hcell[i].den * simpar.scale_d *volcell/Msun *simpar.H0/100L;
	}
	free(aa);
	fclose(fp);
	printf("Opening %s to write gas data\n",outfile);fflush(stdout);
	FILE *wp = fopen(outfile,"w");
	Wait2Start(WGroupSize);
	fwrite(gas,sizeof(GasType),np,wp);
	Wait2Go(WGroupSize);
	fclose(wp);
	free(gas);
}
