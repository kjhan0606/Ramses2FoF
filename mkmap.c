/* 
 * icc -o mkmap mkmap.c -L -lmyram -g -DQUADHILBERT -DNENER=0 -DNVAR=11 -DNPRE=8 -DOUTPUT_PARTICLE_POTENTIAL
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

/*
#define x0 317.L
#define x1 395.L
#define y0 317.L
#define y1 395.L
int nx=22480/78.*(x1-x0);
int ny=22480/78.*(x1-x0);
*/
#define x0 330.L
#define x1 370.L
#define y0 330.L
#define y1 370.L
/*
int nx=22480/78.*(x1-x0);
int ny=22480/78.*(x1-x0);
*/
int nx=4096;
int ny=4096;

#define MIN(a,b) ( (a)<(b)? (a):(b) )
#define MAX(a,b) ( (a)>(b)? (a):(b) )


int main(int argc, char *argv[]){
	size_t i,j,k;
	FILE *fp = fopen(argv[1],"r");

	DmType *dm;
	StarType *star;
	SinkType *sink;
	HydroCellType *hcell;
	struct stat st;
	size_t size, rsize;

	stat(argv[1], &st);
	size = st.st_size;
	size_t dsize;

	if(strcmp(argv[2],"STAR")==0){
		dsize = sizeof(StarType);
		if(size%dsize !=0){
			printf("Error in file size,,,, %s\n", argv[1]);
			exit(99);
		}
		star = (StarType*)malloc(sizeof(char)*size);
		rsize = fread(star, sizeof(char), size, fp);
	}
	else if(strcmp(argv[2],"SINK")==0){
		dsize = sizeof(SinkType);
		if(size%dsize !=0){
			printf("Error in file size,,,, %s\n", argv[1]);
			exit(99);
		}
		sink = (SinkType*)malloc(sizeof(char)*size);
		rsize = fread(sink, sizeof(char), size, fp);
	}
	else if(strcmp(argv[2],"HCELL")==0){
		dsize = sizeof(HydroCellType);
		if(size%dsize !=0){
			printf("Error in file size,,,, %s\n", argv[1]);
			exit(99);
		}
		hcell = (HydroCellType*)malloc(sizeof(char)*size);
		rsize = fread(hcell, sizeof(char), size, fp);
	}
	else if(strcmp(argv[2],"DM")==0){
		dsize = sizeof(DmType);
		if(size%dsize !=0){
			printf("Error in file size,,,, %s\n", argv[1]);
			exit(99);
		}
		dm = (DmType*)malloc(sizeof(char)*size);
		rsize = fread(dm, sizeof(char), size, fp);
	}


	fclose(fp);

	if(rsize != size) printf("Strange size in file %s : %ld\n", argv[1], rsize);



	dptype step = (x1-x0)/nx;
	int ibgwidth;
	int np = size/dsize;



	printf("total np = %d\n", np);

	if(strcmp(argv[2], "HCELL")==0){
		float *dmap = (float*)malloc(sizeof(float)*nx*ny);
		float *Tmap = (float*)malloc(sizeof(float)*nx*ny);
		for(i=0;i<nx*ny;i++) dmap[i] = Tmap[i] = 0;
		for(i=0;i<np;i++){
			int j,k;
			if(hcell[i].ilevel < 13) continue;
			dptype cellsize05 = hcell[i].cellsize*0.5;
			ibgwidth = ceil(cellsize05 /step);
			j = (hcell[i].y-x0)/step;
			k = (hcell[i].z-y0)/step;
			if(j>=-ibgwidth && j < nx + ibgwidth&& k >=-ibgwidth && k < ny+ibgwidth) {
				int ix,iy;
				int jx,jy;
				int ii,jj;
				ix = ((hcell[i].y-cellsize05-x0)/step);
				jx = ((hcell[i].y+cellsize05-x0)/step);
				iy = ((hcell[i].z-cellsize05-y0)/step);
				jy = ((hcell[i].z+cellsize05-y0)/step);
				ix = MAX(0, ix);
				jx = MIN(nx-1, jx);
				iy = MAX(0, iy);
				jy = MIN(ny-1, jy);
				for(jj=iy;jj<jy;jj++){
					for(ii=ix;ii<jx;ii++){
						Tmap[ii+nx*jj] = MAX(Tmap[ii+nx*jj], hcell[i].temp);
						dmap[ii+nx*jj] = MAX(dmap[ii+nx*jj], hcell[i].den);
					}
				}
			}
		}
		FILE *wp = fopen("temp.dat","w");
		fwrite(&nx, sizeof(int), 1, wp);
		fwrite(&ny, sizeof(int), 1, wp);
		fwrite(Tmap, sizeof(float), nx*ny, wp);
		fclose(wp);

		wp = fopen("gas.dat","w");
		fwrite(&nx, sizeof(int), 1, wp);
		fwrite(&ny, sizeof(int), 1, wp);
		fwrite(dmap, sizeof(float), nx*ny, wp);
		fclose(wp);
	}
	else if(strcmp(argv[2],"STAR")==0){
		float *dmap = (float*)malloc(sizeof(float)*nx*ny);
		for(i=0;i<nx*ny;i++) dmap[i] = 0;
		for(i=0;i<np;i++){
			int j,k;
			j = (star[i].y-x0)/step;
			k = (star[i].z-y0)/step;
			if(j >=0 && j < nx && k >=0 && k < ny) dmap[j+nx*k] += 1;
		}
		FILE *wp = fopen("star.dat","w");
		fwrite(&nx, sizeof(int), 1, wp);
		fwrite(&ny, sizeof(int), 1, wp);
		fwrite(dmap, sizeof(float), nx*ny, wp);
		fclose(wp);
	}
	else if(strcmp(argv[2],"DM")==0){
		float *dmap = (float*)malloc(sizeof(float)*nx*ny);
		for(i=0;i<nx*ny;i++) dmap[i] = 0;
		for(i=0;i<np;i++){
			int j,k;
			j = (dm[i].y-x0)/step;
			k = (dm[i].z-y0)/step;
			if(j >=0 && j < nx && k >=0 && k < ny) dmap[j+nx*k] += 1;
		}
		FILE *wp = fopen("dm.dat","w");
		fwrite(&nx, sizeof(int), 1, wp);
		fwrite(&ny, sizeof(int), 1, wp);
		fwrite(dmap, sizeof(float), nx*ny, wp);
		fclose(wp);
	}
	else {
		printf("Ooooooooops.... Nothing happens!\n");
	}

}
