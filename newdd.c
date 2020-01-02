#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>

#include "ramses.h"
#include "params.h"
#include "Memory.h"

int nid=1, myid=0;

#ifdef USE_MPI
#include<mpi.h>

#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ((a)%(b))

#define Wait2Start(WGroupSize) do{\
	if(RANKINGROUP(myid,WGroupSize) != 0 ) {\
		int iget, src,itag=1; \
		src = myid-1;\
		MPI_Status status; \
		MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status);\
	}\
}while(0)

#define Wait2Go(WGroupSize) do{\
	int isend=0,tgt,itag = 1;\
	tgt = myid+1;\
	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) {\
		MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);\
	}\
} while(0)
#define MIN(a,b) ( (a)<(b)? (a):(b) )

#else
#define Wait2Start(WGroupSize) do{}while(0)
#define Wait2Go(WGroupSize) do{}while(0)

#endif



int main(int argc, char **argv){
	size_t i,j,k;
	FILE *fp;
	char infile[190];
	RamsesType ram;
	int istep, icpu;
	int WGroupSize = WGROUPSIZE;
	int sinmul,nsplit;

	istep = atoi(argv[1]);
	nsplit = atoi(argv[2]);
	(void) Make_Total_Memory();

	sprintf(infile,"./output_%.5d/info_%.5d.txt", istep,istep);


#ifdef USE_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);

	sinmul = 1;
	if(myid==0) rd_info(&ram, infile);
	MPI_Bcast(&ram, sizeof(RamsesType), MPI_BYTE, 0, MPI_COMM_WORLD); 
	int mystart, myfinal; 
	int nstep = (ram.ncpu+nid-1)/nid; 
	if(ram.ncpu%nid !=0){ 
		printf("Error %d %d\n", nid, ram.ncpu); 
		MPI_Finalize(); 
		exit(999); 
	} 
	mystart = myid*nstep; 
	myfinal = (myid+1)*nstep; 
	mystart = MIN(ram.ncpu, mystart)+1; 
	myfinal = MIN(ram.ncpu, myfinal)+1;
	printf("P%d has ranges %d : %d\n", myid, mystart, myfinal);
	for(icpu=mystart;icpu<myfinal;icpu++)
#else
	sinmul = 0;
	rd_info(&ram, infile);
	for(icpu=1;icpu<=ram.ncpu;icpu++)
#endif
	{
		/*
		Wait2Start(WGroupSize);
#ifdef USE_MPI
		printf("P%d is at %d  and has a range of %d and %d\n",myid, icpu,mystart,myfinal);
#endif
		Wait2Go(WGroupSize);
		*/
		sprintf(infile,"./output_%.5d/amr_%.5d.out%.5d", istep, istep, icpu);
		Wait2Start(WGroupSize);
		printf("P%d: Opening %s\n", myid, infile);fflush(stdout);
		rd_amr(&ram, infile, NO);
		Wait2Go(WGroupSize);

		sprintf(infile,"./output_%.5d/part_%.5d.out%.5d", istep, istep, icpu);
		Wait2Start(WGroupSize);
		printf("P%d: Opening %s\n", myid, infile);fflush(stdout);
		rd_part(&ram, infile);
		Wait2Go(WGroupSize);


		sprintf(infile,"./output_%.5d/hydro_%.5d.out%.5d", istep, istep, icpu);
		Wait2Start(WGroupSize);
		printf("P%d: Opening %s\n", myid, infile);fflush(stdout);
		rd_hydro(&ram, infile);
		Wait2Go(WGroupSize);

		sprintf(infile,"./output_%.5d/LeafCells_%.5d.out%.5d", istep, istep, icpu);
		/*
		find_leaf(&ram, icpu, infile);
		qsort(ram.hcell, ram.nleafcell, sizeof(HydroCellType), hcellsortx);
		SplitDump(&ram, ram.hcell, ram.nleafcell,HCELL, istep, icpu,sinmul, nsplit);
		*/
		find_leaf_gas(&ram, icpu, infile);
		qsort(ram.gas, ram.nleafcell, sizeof(GasType), gassortx);
		SplitDump(&ram, ram.gas, ram.nleafcell,GAS, istep, icpu,sinmul, nsplit);



		if(icpu==1){
			sprintf(infile,"./output_%.5d/sink_%.5d.out", istep, istep);
			rd_sink(&ram, infile);
			qsort(ram.sink, ram.nsink, sizeof(SinkType), sinksortx);
			SplitDump(&ram, ram.sink, ram.nsink, SINK, istep, icpu,0, nsplit);
			Free(ram.sink);
		}



		{
			DmType *dm = (DmType*)Malloc(sizeof(DmType)*ram.npart, PPTR(dm));
			size_t ndm = 0;
			for(i=0;i<ram.npart;i++){
				if((ram.particle)[i].tp ==0) {
					dm[ndm] = ((DmType*)(ram.particle))[i];
					ndm ++;
				}
			}
			qsort(dm, ndm, sizeof(DmType), dmsortx);
			printf("P%d End of sorting %d dm\n", myid, ndm);
			SplitDump(&ram, dm, ndm, DM, istep, icpu,sinmul, nsplit);


			size_t nstar = 0;
			StarType *star = (StarType*)dm;
			for(i=0;i<ram.npart;i++){
				if((ram.particle)[i].tp !=0) {
					star[nstar] = ((StarType*)(ram.particle))[i];
					nstar ++;
				}
			}
			qsort(star, nstar, sizeof(StarType), starsortx);
			printf("P%d End of sorting %d star\n", myid, ram.nleafcell);
			SplitDump(&ram, star, nstar,STAR, istep, icpu,sinmul, nsplit);
			Free(dm);
	
		}

		if(0){
			FILE *wp = fopen("Test.out","w");
			write_head(wp, ram);
			fclose(wp);
		}
		cleanup_ramses(&ram);
		printf("Current memory stack %d\n", CurMemStack());fflush(stdout);
	}
	return 0;
}
