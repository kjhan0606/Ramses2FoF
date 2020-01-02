#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>

#include "ramses.h"
#include "params.h"
#include "Memory.h"



int main(int argc, char **argv){
	size_t i,j,k;
	FILE *fp;
	char infile[190];
	RamsesType ram;
	int istep, icpu;

	istep = atoi(argv[1]);
	icpu = atoi(argv[2]);

	sprintf(infile,"info_%.5d.txt", istep);
	rd_info(&ram, infile);

	sprintf(infile,"amr_%.5d.out%.5d", istep, icpu);
	rd_amr(&ram, infile, NO);

	sprintf(infile,"part_%.5d.out%.5d", istep, icpu);
	rd_part(&ram, infile);

	sprintf(infile,"sink_%.5d.out%.5d", istep, icpu);
	rd_sink(&ram, infile);




	sprintf(infile,"hydro_%.5d.out%.5d", istep, icpu);
	rd_hydro(&ram, infile);

	/* This is to output the leaf cell information */
	sprintf(infile,"LeafCells_%.5d.out%.5d", istep, icpu);

	find_leaf(&ram, icpu, infile);


	qsort(ram.sink, ram.nsink, sizeof(SinkType), sinksortx);

	SplitDump(&ram, ram.sink, ram.nsink, SINK, istep, icpu);

	qsort(ram.hcell, ram.nleafcell, sizeof(HydroCellType), hcellsortx);
	SplitDump(&ram, ram.hcell, ram.nleafcell,HCELL, istep, icpu);
	{
		DmType *dm = (DmType*)Malloc(sizeof(DmType)*ram.npart,PPTR(dm));
		size_t ndm = 0;
		for(i=0;i<ram.npart;i++){
			if(dm[i].tp ==0) {
				dm[ndm] = ((DmType*)(ram.particle))[i];
				ndm ++;
			}
		}
		qsort(dm, ndm, sizeof(DmType), dmsortx);
		SplitDump(&ram, dm, ndm, DM, istep, icpu);


		size_t nstar = 0;
		StarType *star = (StarType*)dm;
		for(i=0;i<ram.npart;i++){
			if(star[i].tp !=0) {
				star[nstar] = ((StarType*)(ram.particle))[i];
				nstar ++;
			}
		}
		qsort(star, nstar, sizeof(StarType), starsortx);
		SplitDump(&ram, star, nstar,STAR, istep, icpu);
		Free(dm);

	}

	if(1){
		FILE *wp = fopen("Test.out","w");
		write_head(wp, ram);
		fclose(wp);
	}

	cleanup_ramses(&ram);
	return 0;
}
