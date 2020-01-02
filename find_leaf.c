#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include "ramses.h"
#include "Memory.h"

#define MAX(a,b) ( (a)>(b) ? (a): (b) )



int find_leaf(RamsesType *ram, int myid, char *outfile){

	int i,j,k,ig;

	MeshType *mesh = &(ram->mesh);
	int *numbl = mesh->numbl;
	int *headl = mesh->headl;
	HydroType *hydro = &(ram->hydro);
	int ldim, ncpu, nlevelmax;
	size_t ncell;
	int ngridmax;
	dptype dx;
	int ndim;
	int ncoarse = ram->ncoarse;
	int *next = (mesh->next);

	ndim = ram->ndim;
	ncpu = ram->ncpu;
	nlevelmax = ram->nlevelmax;
	ngridmax = ram->ngridmax;
	dptype boxlen_ini = ram->boxlen_ini;
	dptype hubin = ram->H0 / 100.L;
	int nboundary = ram->nboundary;
	int ilevel,istart;
	int *numbb = mesh->numbb;
	int *headb = mesh->headb;
	int igrid, ind;
	int twotondim = ram->twotondim;
	size_t iskip;
	int *son = mesh->son;
	int CPI[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1}};
	int *cpi;

	ncell = ncoarse + twotondim * ngridmax;


	dptype Lbox = boxlen_ini;
	dptype Lbox05 = Lbox/2;
	dptype smallr = ram->smallr;
	int imetal = ram->imetal;


	dptype *uold = hydro->uold;




	int ntot = 0;
	int ncache = 512;
//	int *ind_leaf = (int*)Malloc(sizeof(int)*ncache, PPTR(ind_leaf));
	int ileaf=0;


	int nleaf = 0;
	int ibound;
	for(ibound=1;ibound <= nboundary+ncpu;ibound++){
		if(ibound == myid){
			for(ilevel=1;ilevel<=nlevelmax;ilevel++){
				if(ibound <=ncpu) {
					ncache = numbl[ibound -1 + (ilevel-1)*ncpu]; 
					istart = headl[ibound -1 + (ilevel-1)*ncpu]; 
				}
				else {
					ncache = numbb[ibound -1 + (ilevel-1)*nboundary]; 
					istart = headb[ibound -1 + (ilevel-1)*nboundary]; 
				}
				if(ncache>0) {
					int *ind_grid = (int*)Malloc(sizeof(int)*ncache, PPTR(ind_grid));
					igrid = istart;
					for(i=0;i<ncache;i++){
						ind_grid[i] = igrid;
						igrid = next[igrid-1];
					}
					for(ind=0;ind<twotondim;ind++){
						iskip = ncoarse + ind*ngridmax;
						for(i=0;i<ncache;i++){
							if(son[ind_grid[i]-1+iskip] ==0){
								nleaf ++;
							}
						}
					}
					Free(ind_grid);
					ntot+=ncache;
				}
			}
		}
	}
	printf("Total nleaf= %d\n", nleaf);
	HydroCellType *hcell;
	hcell = (ram->hcell = (HydroCellType*)Malloc(sizeof(HydroCellType)*nleaf, PPTR(ram->hcell)));

	ileaf = 0;
	for(ibound=1;ibound <= nboundary+ncpu;ibound++){
		if(ibound == myid){
			for(ilevel=1;ilevel<=nlevelmax;ilevel++){
				dptype dx = pow(0.5, ilevel);
				if(ibound <=ncpu) {
					ncache = numbl[ibound -1 + (ilevel-1)*ncpu]; 
					istart = headl[ibound -1 + (ilevel-1)*ncpu]; 
				}
				else {
					ncache = numbb[ibound -1 -ram->ncpu + (ilevel-1)*nboundary]; 
					istart = headb[ibound -1 -ram->ncpu + (ilevel-1)*nboundary]; 
				}
				if(ncache>0) {
					int *ind_grid = (int*)Malloc(sizeof(int)*ncache, PPTR(ind_grid));
					igrid = istart;
					nleaf = 0;
					for(i=0;i<ncache;i++){
						ind_grid[i] = igrid;
						igrid = next[igrid-1];
					}
					for(ind=0;ind<twotondim;ind++){
						iskip = ncoarse + ind*ngridmax;
						for(i=0;i<ncache;i++){
							if(son[ind_grid[i]-1+iskip] ==0){
								nleaf ++;
							}
						}
					}
					if(nleaf >0) {
						j = 0;
						for(ind=0;ind<twotondim;ind++){
							iskip = ncoarse + ind*ngridmax;
							cpi = CPI[ind];
							for(i=0;i<ncache;i++){
								if(son[ind_grid[i]-1+iskip]==0){
									hcell[ileaf].den = uold[ind_grid[i]-1+iskip];
									hcell[ileaf].x  = (mesh->xg[ind_grid[i]-1+ngridmax*0] + (cpi[0]-0.5)*dx)* Lbox;
									hcell[ileaf].y  = (mesh->xg[ind_grid[i]-1+ngridmax*1] + (cpi[1]-0.5)*dx)* Lbox;
									hcell[ileaf].z  = (mesh->xg[ind_grid[i]-1+ngridmax*2] + (cpi[2]-0.5)*dx)* Lbox;
									hcell[ileaf].vx = uold[ind_grid[i]-1+iskip+1*ncell]/MAX(uold[ind_grid[i]-1+iskip],smallr);
									hcell[ileaf].vy = uold[ind_grid[i]-1+iskip+2*ncell]/MAX(uold[ind_grid[i]-1+iskip],smallr);
									hcell[ileaf].vz = uold[ind_grid[i]-1+iskip+3*ncell]/MAX(uold[ind_grid[i]-1+iskip],smallr);
									hcell[ileaf].cellsize = dx*Lbox; /* in cMpc/h */
									hcell[ileaf].metallicity = uold[ind_grid[i]-1+iskip + ncell*imetal]/hcell[ileaf].den/0.02;
									dptype ekk = 0;
									ekk += 0.5*hcell[ileaf].den*hcell[ileaf].vx*hcell[ileaf].vx;
									ekk += 0.5*hcell[ileaf].den*hcell[ileaf].vy*hcell[ileaf].vy;
									ekk += 0.5*hcell[ileaf].den*hcell[ileaf].vz*hcell[ileaf].vz;

									dptype Temp = uold[ind_grid[i]-1+iskip+(ndim+1)*ncell];
									dptype err = 0;
#if NENER>0
									int irad;
									for(irad=0;irad<nener;irad++){ 
										err += uold[ind_grid[i]-1 + iskip + ncell*(inener+irad)]; 
									}
#endif
									dptype emag = 0;
#ifdef SOLVERmhd 
									for(idim=0;idim<ndim;idim++){ 
										emag += 0.125L*TWICE(uold[ind_grid[i]-1+iskip+ncell*(idim+ndim+2)] + 
												uold[ind_grid[i]-1+iskip+ncell*(idim+nvar)]); 
									}
#endif
									hcell[ileaf].temp = (ram->gamma-1)*(uold[ind_grid[i]-1+iskip+ncell*(ndim+1)]-ekk-err-emag);
									hcell[ileaf].temp *= ram->scale_T2; /* in K/mu */
									hcell[ileaf].vx *= ram->kmscale_v; /* in km/sec */
									hcell[ileaf].vy *= ram->kmscale_v; /* in km/sec */
									hcell[ileaf].vz *= ram->kmscale_v; /* in km/sec */

									hcell[ileaf].H  = uold[ind_grid[i]-1+iskip+6*ncell]/hcell[ileaf].den;
									hcell[ileaf].O  = uold[ind_grid[i]-1+iskip+7*ncell]/hcell[ileaf].den;
									hcell[ileaf].Fe = uold[ind_grid[i]-1+iskip+8*ncell]/hcell[ileaf].den;
									hcell[ileaf].ilevel = ilevel;

									ileaf++;
								}
							}
						}
					}
					Free(ind_grid);
				}
			}
		}
	}
	ram->nleafcell = ileaf;
	printf("Total number of leaf cells are %d from %d\n", ileaf,ntot);
	if(0){
		FILE *wp = fopen(outfile,"w");
		for(i=0;i<ileaf;i++){
			/*
			if(hcell[i].metallicity > 1.E-5 && hcell[i].temp < 1.E3)
			*/
			if(hcell[i].temp < 1.E3)
			{
				fprintf(wp, "%d %g %g %g %g %g\n",i, hcell[i].metallicity, hcell[i].temp, hcell[i].H, hcell[i].O, hcell[i].Fe);
			}
		}
		fclose(wp);
	}
}
