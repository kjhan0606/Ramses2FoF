/* 
 * icc -o corrinfo corrinfo.c -L. -lmyram -g -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL
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

idtype id=0;

dptype Mpc = 3.0857E24;
dptype Msun = 2E33;
RamsesType simpar;
RamsesType read_head(FILE*);

int hcell2gas(char *, char*);

int main(int argc, char *argv[]){
	int nfile,nstep;

	nstep=atoi(argv[1]);
	nfile=atoi(argv[2]);
	char infile[190],outfile[190];
	FILE *fp;
	int i;
	id = 0;
	sprintf(infile,"info_%.5d.txt", nstep);
	RamsesType ram;
	rd_info(&ram, infile);
	ram.boxlen_ini = 717.229040;


	RamsesType ram2;
	sprintf(infile,"HR5.%.5d.00000.info", nstep);
	fp = fopen(infile,"r");
	ram2 = read_head(fp);
	fclose(fp);
	units(&ram);
	ram.nx = ram.ny=ram.nz = 1;
	ram.nstep = 107492;
	ram.ramses_sizeof = sizeof(ram);

	for(i=0;i<nfile;i++){

		ram.icpu = i;
		ram.xmin = ram.boxlen_ini/1024. * i;
		ram.xmax = ram.boxlen_ini/1024. * (i+1);



		sprintf(outfile,"HR5.%.5d.%.5d.info",nstep,  i);
		FILE *wp = fopen(outfile,"w");
		write_head(wp, ram);
		fwrite(&ram,sizeof(RamsesType), 1, wp);
		fclose(wp);
	}
}
