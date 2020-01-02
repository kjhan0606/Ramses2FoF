/* This is a header file to read RAMSES data in C. */

/* You have to read these lines carefully before running it */
#define NDIM 3
#define NLEVELS 13
#ifdef SINGLE_PRECISION
typedef float dptype;
typedef int idtype;
#else
typedef double dptype;
typedef long long idtype;
#endif
#ifdef QUADHILBERT
typedef long double qdptype;
#else
typedef double qdptype;
#endif
#define MPI_SIZE_T MPI_LONG
/* You have to read these lines carefully before running it */



/*
#define NSPLIT 1024
*/
/*
int NSPLIT= 512;
*/

#define Mpc  (3.0857E24L)
#define kpc  (3.0857E21L)
#define Msun  (2E33L)
#define rhoc  (1.8800000e-29L)
#define mH (1.6600000e-24L)
#define kB (1.3806200e-16L)
#define X (0.76L)


#define IDOFFSET 5000000000


#define IRandNumSize 4

#define YES 1
#define NO 0


typedef struct CommType{
	int ngrid, npart, *igrid, *f;
	dptype *u, *up;
	idtype *fp;
}CommType;

/******************************************************************************/
#define DM 0
#define STAR 1
#define SINK 2
#define HCELL 3
#define GAS 4
typedef struct PmType{
	dptype x,y,z,vx,vy,vz; /* in cMpc/h and km/second */
	dptype mass; /* in unit of Msun/h */
	dptype metal,tp,zp,mass0,tpp,indtab; /* in ramses units */
	idtype id; /* if id < 0, it is a sink particle */
#ifdef OUTPUT_PARTICLE_POTENTIAL
	dptype potent;
#endif
	int levelp; /* Initial level of the particle */
}PmType;

typedef PmType DmType;
typedef PmType StarType;

typedef struct SinkType{
	dptype x,y,z,vx,vy,vz; /* in cMpc/h and km/second */
	dptype mass;  /* in unit of Msun/h */
	dptype tbirth; /* in ramses units */
	dptype Jx,Jy,Jz; /* in [Msun/h km/second kcp] */
	dptype Sx,Sy,Sz; /* in units angle*/
	dptype dMsmbh,dMBH_coarse,dMEd_coarse; /* in unit of Msun/h */
	dptype Esave,Smag,eps; /* in ramses units */
//	dptype stat[2*NDIM+1][NLEVELS];  /* in ramses units */
	int id;
} SinkType;
typedef struct AGNType{ 
	dptype x,y,z,vx,vy,vz; /* in cMpc/h and km/second */ 
	dptype mass;  /* in unit of Msun/h */ 
	dptype tbirth; /* in ramses units */ 
	dptype Jx,Jy,Jz; /* in ramses units */ 
	dptype Sx,Sy,Sz; /* in ramses units */ 
	dptype dMsmbh,dMBH_coarse,dMEd_coarse, Esave,Smag,eps; /* in ramses units */ 
	dptype dMBHovert,Lbol; /* in ramses unit */
	idtype id;
} AGNType;


typedef struct HydroCellType{
	dptype x,y,z,cellsize; /* the center position and size of cell in cMpc/h */
	float vx,vy,vz; /* physical velocity in km/second */
	float den; /* simulation density of the cell. 
				  A factor of scale_d should be multiplied to get the real density in gr/cm^3 */
	float temp; /* gas temperature in K/mu */
	float metallicity;  /* metallicity Z */
	float Fe, H, O; /* in mass fractions */
	int ilevel;
} HydroCellType;
typedef struct GasType{ 
	dptype x,y,z,cellsize; /* the center position and size of cell in cMpc/h */ 
	float vx,vy,vz; /* physical velocity in km/second */ 
	float den; /* simulation density of the cell. 
				  A factor of scale_d should be multiplied to get the real density in gr/cm^3 */ 
	float temp; /* gas temperature in K/mu */ 
	float metallicity;  /* metallicity Z */
	float Fe, H, O; /* in mass fractions */ 
	int ilevel; 
	float mass; /* Msun/h */
	idtype indx;
} GasType;

/******************************************************************************/

typedef struct MeshType{
	int *ilevel,  *headl, *taill,*numbl;
	long *numbtot;
	int *flag1, *flag2, *son, *father;
	int *nbor;
	int headf, tailf, numbf;
	int *headb, *tailb, *numbb;
	dptype *dtold, *dtnew, *xg;
	int *prev, *next, *cpu_map, *cpu_map2;
//	CommType *active;
	qdptype *bound_key;
	dptype *bound_key_restart;
}MeshType;

typedef struct HydroType{
	dptype *uold, *unew;
}HydroType;



typedef struct RamsesType{
	int npart, nsink,nleafcell,ngas;
	int nindsink;
	int ncpu, icpu, ndim;
	int nrestart,nrestart_quad;
	int localseed[IRandNumSize], nstar_tot;
	dptype mstar_tot,mstar_lost;
	int levelmin, nlevelmax,nstep_coarse;
	idtype ngridmax;
	dptype boxlen, time, aexp,amax, H0, omega_m, omega_l, omega_k,omega_b;
	dptype unit_l, unit_d, unit_t;
	int nx,ny,nz;
	int ndomain;
	int overload;
	int ngrid_current,iout, ifout;
	int nstep;
	dptype constant, mass_tot_0, rho_tot;
	dptype omega_h0, aexp_ini, boxlen_ini, hexp, aexp_old, epot_tot_int;
	dptype mass_sph;
	int nbinodes, twotondim;
	int used_mem, used_mem_old, used_mem_tot;
	dptype epot_tot_old,t;
	char ordering[129];
	int cosmo;
	dptype scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2, scale_m,scale_Gyr;
	dptype mpcscale_l, kmscale_v;

	int nvar,nboundary, ncoarse, nener;
	int imetal, ichem, nelt, inener;

	int neq_chem,rt;
	dptype gamma_rad[512];
	dptype gamma,smallr;

	MeshType mesh;
	HydroType hydro;

	dptype xmin,xmax,ymin,ymax,zmin,zmax;



	PmType *particle;
	SinkType *sink;
	HydroCellType *hcell;
	GasType *gas;
	int ramses_sizeof;

}RamsesType;


#define WHERESTR  "[file %s, line %d]: "
#define WHEREARG  __FILE__, __LINE__
#define DEBUGPRINT2(...)       fprintf(stderr, __VA_ARGS__)
#define DEBUGPRINT(_fmt, ...)  DEBUGPRINT2(WHERESTR _fmt, WHEREARG, __VA_ARGS__)
#define DEBUGPRINT0(_fmt) DEBUGPRINT2(WHERESTR _fmt, WHEREARG)



#define F77read(a, size, nmem, fp) do{\
	int _nn,chip;\
	_nn=fread(&chip, sizeof(int), 1,fp);\
	if(chip!=size*nmem) {\
		DEBUGPRINT("Error reading "#a" "#size" "#nmem"  %d :  %d @ %p\n", chip, size*nmem, fp);\
		exit(99);\
	}\
	fread(a, size, nmem,fp);\
	fread(&chip, sizeof(int), 1,fp);\
} while(0)


#define GetPart(xbuff, size, npart, fp, pm, type,mem) do{\
	int iii;\
	F77read(xbuff, size, npart, fp);\
	for(iii=0;iii<npart;iii++) pm->type[iii].mem = xbuff[iii];\
}while(0)


void rd_info(RamsesType *, char *);

int rd_amr(RamsesType *, char *, int);
int rd_hydro(RamsesType *, char *);
int rd_part(RamsesType *, char *);
int rd_sink(RamsesType *, char *);
void rd_info(RamsesType *, char *);
int find_leaf(RamsesType *, int, char *);
int find_leaf_gas(RamsesType *, int, char *);
void GetHydroQ(RamsesType *, int *, int , HydroCellType *);

void cleanup_mesh(RamsesType*, int);
void cleanup_ramses(RamsesType*);

void  SplitDump(RamsesType *, const void *, int , int, int, int, int, int);

int sinksortx(const void*, const void*);
int sinksortid(const void*, const void*);
int agnsortid(const void*, const void*);
int hcellsortx(const void*, const void*);
int gassortx(const void*, const void*);
int starsortx(const void*, const void*);
int dmsortx(const void*, const void*);
void units(RamsesType *);
