/* NAME: RD_PaRT
 * PURPOSE: This proceduer reads particles from a RAMSES PART file. And it is rewritten from the IDL version of rd_part.pro.
 */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>

#include "ramses.h"
#include "Memory.h"




int rd_part(RamsesType *ram, char *infile){
	FILE *fp = fopen(infile,"r");
	size_t i;
	int chip;
	int npartp,mstar;
	F77read(&(ram->ncpu), sizeof(int), 1, fp);
	F77read(&(ram->ndim), sizeof(int), 1, fp);
	F77read(&(ram->npart), sizeof(int), 1, fp);
	npartp = ram->npart;
	F77read(&(ram->localseed), sizeof(int), IRandNumSize, fp);
	F77read(&(ram->nstar_tot), sizeof(int), 1, fp);
	F77read(&(ram->mstar_tot), sizeof(dptype), 1, fp);
	F77read(&(ram->mstar_lost), sizeof(dptype), 1, fp);
	F77read(&(ram->nsink), sizeof(int), 1, fp);
	dptype *xbuff = (dptype*)Malloc(sizeof(dptype)*npartp,PPTR(xbuff));
	ram->particle = (PmType*)Malloc(sizeof(PmType)*npartp,PPTR(ram->particle));


	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,x);



	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,y);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,z);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,vx);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,vy);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,vz);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,mass);
	idtype *idbuff = (idtype*)Malloc(sizeof(idtype)*npartp, PPTR(idbuff));
	GetPart(idbuff,sizeof(idtype), npartp, fp, ram,particle,id);

	int *ibuff = (int*)Malloc(sizeof(int)*npartp,PPTR(ibuff));
	GetPart(ibuff,sizeof(int), npartp, fp, ram,particle,levelp);
#ifdef OUTPUT_PARTICLE_POTENTIAL
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,potent);
#endif
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,tp);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,zp);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,tpp);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,mass0);
	GetPart(xbuff,sizeof(dptype), npartp, fp, ram,particle,indtab);
	fclose(fp);
	Free(ibuff);
	Free(idbuff);
	Free(xbuff);
	PmType *part = ram->particle;
	int ipartp = 0;
	int nstar = 0;
	for(i=0;i<npartp;i++){
		if(part[i].id >0) {
			part[ipartp] = part[i];
			part[ipartp].mass = part[i].mass* ram->scale_m;
			part[ipartp].x = part[i].x * ram->mpcscale_l;
			part[ipartp].y = part[i].y * ram->mpcscale_l;
			part[ipartp].z = part[i].z * ram->mpcscale_l;
			part[ipartp].vx = part[i].vx * ram->kmscale_v;
			part[ipartp].vy = part[i].vy * ram->kmscale_v;
			part[ipartp].vz = part[i].vz * ram->kmscale_v;
//			part[ipartp].mass0 = part[i].mass0* ram->scale_m;
			if(part[ipartp].tp !=0) nstar ++;
			ipartp ++;
		}
	}
	ram->npart = (npartp= ipartp);
	part = (ram->particle = (PmType*)Realloc(ram->particle, sizeof(PmType)*npartp));

	printf("Total star and dm particles are %d %d from total np=  %d\n", nstar, npartp-nstar, npartp);
	if(0){
		float maxmass = -1.e20;
		for(i=0;i<npartp;i++){
			if(part[i].id <0){
				if(maxmass < part[i].mass) maxmass = part[i].mass;
				chip ++;
			}
		}
		printf("Maximum sink particle mass is %g\n", maxmass*ram->scale_m);
	}
	return npartp;
}
