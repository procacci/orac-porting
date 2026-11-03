#ifndef AGBNP_PRIV_H
#define AGBNP_PRIV_H

#include "agbnp.h"

#ifdef SINGLE_PREC
typedef float float_a;
#else
typedef double float_a;
#endif

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif 

#define mysqrt(x)   (sqrt(x))
#define mypow(x,y)  (pow(x,y))

#ifndef AGBNP_NEARNEIGHBORS
#define AGBNP_NEARNEIGHBORS (32)
#endif

#ifndef AGBNP_FARNEIGHBORS
#define AGBNP_FARNEIGHBORS (256)
#endif

#ifndef AGBNP_MIN_VOLA
#define AGBNP_MIN_VOLA (0.1)
#endif
#ifndef AGBNP_MIN_VOLB
#define AGBNP_MIN_VOLB (1.0)
#endif

#ifndef AGBNP_MINC_VOLA
#define AGBNP_MINC_VOLA (0.2)
#endif
#ifndef AGBNP_MINC_VOLB
#define AGBNP_MINC_VOLB (2.0)
#endif


#ifndef AGBNP_MAX_OVERLAP_LEVEL
#define AGBNP_MAX_OVERLAP_LEVEL (12)
#endif

#ifndef AGBNP_RADIUS_INCREMENT
#define AGBNP_RADIUS_INCREMENT (0.5)
#endif

#define pi (3.14159265359)

/* conversion factors from spheres to Gaussians */
#define KFC (2.2269859253)
#define PFC (2.5);

#ifdef AGBNP_HB
/* HB types */
#define AGBNP_HB_INACTIVE  (0)
#define AGBNP_HB_POLARH    (1)
#define AGBNP_HB_SPACCPT  (10)
#define AGBNP_HB_SP2ACCPT (20)
#define AGBNP_HB_SP3ACCPT (30)
/* bond length for water sites pseudo-atoms */
#define AGBNP_HB_LENGTH (2.5)
/* radius of water site pseudo-atom */
#define AGBNP_HB_RADIUS (1.4)
#endif

/* a structure that holds 3D Gaussian parameters */
typedef struct gparm_ {
  float_a a; /* Gaussian exponent */
  float_a p; /* Gaussian prefactor */
  float_a c[3]; /* center */
} GParm;

#ifdef AGBNP_HB
/* a pseudo atom for water sites */
typedef struct wat_ {
  float_a pos[3]; /* position */
  float_a r;    /* radius */
  int nparents; /* number of parents */
  int parent[4]; /* parents */
  float_a volume; /* the volume of the ms sphere (4/3) pi r^3 */
  float_a sp;       /* free_volume/volume where the free volume calculated 
		      using only heavy atoms */
  float_a free_volume; /* the volume that doesn't overlap with real atoms */
} WSat;
#endif

/* an atom for the molecular surface calculation */
typedef struct msat_ {
  float_a pos[3]; /* position */
  float_a r;    /* radius */
  int parent1; /* first parent */
  int parent2; /* second parent */
  float_a s;       /* effective density */
} MSat;

typedef struct AGBworkdata_ {
  
  int natoms;
  float_a *vols; /* unscaled atomic volumes */
  float_a *volumep; /* scaled volumes */
  float_a *dera, *deru; /* used in the calculation of derivatives of 
			  AGB energy */
  float_a *q2ab; /* q^2+dera*b */
  float_a *derv; /* used in the calculation of the vdw part of
		   the np energy */
  float_a *abrw;   /* alpha*brw */
  float_a *br1_swf_der; /* derivative of inverse Born radius switching 
			  function */
  int *isheavy;       /* 0 if hydrogen */
  int *nbiat; /* neighbors of atom iat */
  float_a *br1; /* inverse Born radii */
  float_a *br;  /* Born radii */
  float_a *brw; /* 3b^2/(b+rw)^4 for np derivative */
  float_a *alpha; /* ideal alpha + correction alpha */
  float_a *delta; /* ideal delta + correction delta */
  float_a *galpha; /* Gaussian exponents (= kfc/r[i]^2) */
  float_a *gprefac; /* Gaussian prefactors (= pfc) */
  float_a *sp; /* scaling volume factors */
  void **nbdata; /* a buffer to hold pointers to PBC translation vectors */
  float_a *xj, *yj, *zj; /* coordinate buffers to hold PBC-adjusted neighbor positions */

  int *isvfrozen; /* 1 if atom belongs to frozen set AND whose near neighbors 
		     (heavy atoms) are all frozen atoms (constant self-volume).
		     0 otherwise. Frozen hydrogens are logically part of this 
		     set. */
  int *isvvfrozen; /* 1 if all near neighbors are vfrozen */
  int *isbfrozen; /* 1 if atom belongs to frozen set AND whose heavy atom 
		     neighbors (near and far) are all part of the vfrozen 
		     set (constant Born radius). 0 otherwise. */
  float_a *volumep_constant; /* holds the constant part of the self-volume */
  float_a *br1_const; /* holds the constant part of the inverse Born radius */ 
  float_a *gbpair_const; /* holds constant part of GB pair sum */
  float_a *dera_const; /* holds the constant part of Ai's */

  float_a *deru_c; /* holds the constant part of Ui's */
  float_a *derv_c; /* holds the constant part of Vi's */
  

  int nq4cache;    /* i4() memory cache size */
  float  *q4cache; /* a memory cache to store i4() stuff */

  int nv2cache;   /* 2-body memory cache size */
  float  *v2cache; /* a memory cache to store 3-body 
		      overlaps stuff */

  int *nl_indx;   /* index buffer used in neigh. list reordering */ 
  float_a *nl_r2v;/* distance buffer used in neigh. list reordering */ 

  int *nlist;  /* temporary arrays of length 'natoms' for sorting, etc. */
  int *js;
  NeighVector *pbcs;
  void **datas;

  NeighList *near_nl; /* near (d<Ri+Rj) neighbor list for heavy atoms */
  NeighList *far_nl;  /* far Ri+Rj<d<cutoff neigh.list */

  float_a (*dgbdr_h)[3]; /* gradient of GB energy */
  float_a (*dvwdr_h)[3]; /* gradient of vdW energy */

  float_a (*dgbdr_c)[3]; /* constant component of gradient of GB energy */
  float_a (*dvwdr_c)[3]; /* constant component of gradient of vdW energy */

} AGBworkdata;


typedef struct CAVworkdata_ {
  
  int natoms;
  float_a *vols; /* unscaled atomic volumes */
  float_a *volumep; /* scaled volumes */
  float_a *surf_area; /* surface areas */
  int *isheavy;       /* 0 if hydrogen */
  int *nbiat; /* neighbors of atom iat */
  float_a *gamma; /* ideal gamma + correction gamma */
  float_a *gammap; /* gamma corrected by derivative of surface area
			    switching function (for derivative calculation) */
  float_a *galpha; /* Gaussian exponents (= kfc/r[i]^2) */
  float_a *gprefac; /* Gaussian prefactors (= pfc) */

  int *isvfrozen; /* 1 if atom belongs to frozen set AND is surrounded by
		     frozen atoms (constant self-volume). 0 otherwise. */
  int *isvvfrozen; /* 1 if all neighbors are vfrozen */
  float_a *volumep_constant; /* holds the constant part of the self-volume */
  float_a *surfarea_constant;/* holds the constant part of the surface area */
  int nv2cache;   /* 2-body memory cache size */
  float  *v2cache; /* a memory cache to store 3-body 
		      overlaps stuff */

  void **nbdata; /* a buffer to hold pointers to PBC translation vectors */
  float_a *xj, *yj, *zj; /* coordinate buffers to hold PBC-adjusted neighbor 
			   positions */

  NeighList *near_nl; /* near (d<Ri+Rj) neighbor list for heavy atoms */

  int *nl_indx;   /* index buffer used in neigh. list reordering */ 
  float_a *nl_r2v;/* distance buffer used in neigh. list reordering */ 

  int *nlist;  /* temporary arrays of length 'natoms' for sorting, etc. */
  int *js;
  NeighVector *pbcs;
  void **datas;

  float_a (*decav_h)[3];
  float_a (*decav_const)[3];


} CAVworkdata;


typedef struct AGBNPdata_ {
  int in_use; /* TRUE if in use */

  int natoms; /* Numb. of atoms. Atoms are numbered from 0 to natoms-1 */

  float_a *x, *y, *z; /* coordinates */
  float_a *r;  /* vdw radii (define solute volume for reaction field
		 calculation) */

  float_a *rcav; /* cavity radii (define solute cavity for surface area
		   calculation) */

  float_a *charge; /* the partial charge of each atom */

  float_a *igamma, *sgamma; /* np parameters */
  float_a *ialpha, *salpha;
  float_a *idelta, *sdelta;

  float_a *ab; /* screening parameters */

#ifdef AGBNP_HB
  int *hbtype;          /* HB atom type 0=inactive, 1=hydrogen donor,
			   10=sp acceptor, 20=sp2 acceptor, 30=sp3 acceptor */
  float_a *hbcorr;       /* HB to solvent correction strength */
#endif

  int nheavyat, *iheavyat; /* number of and list of heavy atoms */
  int nhydrogen, *ihydrogen; /* number of and list of hydrogen atoms */
  int ndummy, *idummy; /* number of and list of dummy atoms */

  int do_frozen; /* 1 if it should turn on optimizations relating to
		    frozen atoms */
  int *isfrozen; /* 1 if a frozen atom, 0 otherwise */

  float_a dielectric_in, dielectric_out; /* dielectric constants */

  NeighList *neigh_list; /* Verlet neighbor list (NULL if not using) */
  NeighList *excl_neigh_list; /* Verlet neighbor list of excluded neighbors
				 (NULL if not using) */

  int dopbc;            /* > 0 if doing PBC's */
  int nsym;             /* number of space group operations */
  int docryst;          /* > 0 if doing crystal PBC */
  int ssize;            /* size of coordinate buffers xs, ys, zs */
  float_i *xs, *ys, *zs; /* coordinates of symmetric images (for crystal
			   PBC's) */
  float_a (*rot)[3][3];  /* cell rotation matrixes */

  NeighList *conntbl; /* connection table, NULL if not using dummy atoms
			 for molecular surface */
#ifdef AGBNP_HB
  int do_w;             /* not zero if doing water sites */
  int nwsat;            /* number of water sites pseudo atoms */
  int wsat_size;        /* size of ws atom list (wsat) */
  WSat *wsat;            /* list of ws atoms in the system */
  NeighList *wsat_list; /* list of ms atoms for each real atom */
#endif

  int do_ms;            /* not zero if computing molecular surface */
  int nmsat;            /* number of molecular surface pseudo atoms */
  int msat_size;        /* size of ms atom list (msat) */
  MSat *msat;           /* list of ms atoms in the system */
  NeighList *msat_list; /* list of ms atoms for each real atom */

  AGBworkdata *agbw; /* AGB work data */
  CAVworkdata *cavw; /* Cavity work data */
  
  int nprocs;  /* number of OPENMP threads (>1 if multithreading) */
  int maxprocs;  /* max number of OPENMP threads */
  AGBworkdata **agbw_p; /* store pointers to thread local copies of */
  CAVworkdata **cavw_p; /* of work data                             */   
  
  int doalt;
  float_a *occupancy;   /* occupancy of alternate conformation */
  int *alt_site;       /* id of alternate conf. site (starting from 1)*/
  int *alt_id;         /* id of alternate conf. in site (starting from 1)*/

} AGBNPdata;


/* 1D and 2D lookup tables */
typedef struct c1table_  {
  int n;
  float_a dx;
  float_a dxinv;
  float_a yinf;
  float_a (*data)[4]; /* for each node: a3, a2, a1, a0 */ 
} C1Table;
typedef struct c1table2d_  {
  int ny;           /* number of nodes */
  float_a dy;       /* spacing */
  float_a dyinv;    /* inverse of spacing */
  C1Table **table;  /* array of tables along x */
} C1Table2D;


/*                                            *
 * macros                                     *
 *                                            */
#ifdef NOTNOW
#define I4_FAR( rij, Ri, Rj) ( \
  _i4rij2 = rij*rij, \
  _i4u1 = rij+Rj,  \
  _i4u2 = rij-Rj,  \
  _i4u3 = _i4u1*_i4u2,  \
  _i4u4 = 0.5*log(_i4u1/_i4u2),  \
  _i4dr4 = _i4twopi*( (Rj/(rij*_i4u3))*(1. - 2.*_i4rij2/_i4u3 ) + _i4u4/_i4rij2 ), \
  _i4twopi*(Rj/_i4u3 - _i4u4/rij) )
#endif
float_a i4p(float_a rij, float_a Ri, float_a Rj, 
	    float_a *dr);
#define I4_FAR( rij, Ri, Rj) ( i4p((rij),(Ri),(Rj),&_i4dr4))


/*                                            *
 * function prototypes of local functions     *
 *                                            */

int agbnp_reset(AGBNPdata *data);
int agbnp_reset_agbworkdata(AGBworkdata *agbw);
int agbnp_allocate_agbworkdata(int natoms, AGBNPdata *agb, AGBworkdata *agbw);
int agbnp_delete_agbworkdata(AGBworkdata *agbw);
int agbnp_init_agbworkdata(AGBNPdata *agbdata, AGBworkdata *agbw);
int agbnp_reset_cavworkdata(CAVworkdata *cavw);
 int agbnp_allocate_cavworkdata(int natoms, AGBNPdata *agb, CAVworkdata *cavw);
 int agbnp_delete_cavworkdata(CAVworkdata *cavw);
 int agbnp_init_cavworkdata(AGBNPdata *agbdata, CAVworkdata *cavw);
 int agbnp_tag_ok(int tag);


 float_a i4(float_a rij, float_a Ri, float_a Rj, float_a *dr);
#ifdef VOLUME_OVERLAP
 float_a volume_overlap(float_a (*x)[3], float_a *r,
			     float_a (*dr)[3], float_a *dR, 
			     float_a (*d2rR)[AGBNP_MAX_OVERLAP_LEVEL][3]);
#endif
 float_a ogauss_2body_smpl(float_a d2, float_a p1, float_a p2, 
				float_a c1, float_a c2);
 float_a swf_area(float_a x, float_a *fp);
 float_a swf_vol3(float_a x, float_a *fp, float_a *fpp, float_a a, float_a b);
 float_a ogauss_incremental(int n, float_a (*x)[3], float_a *a, float_a *p,
				 float_a *r,
				 float_a anm1, float_a pnm1, float_a cnm1[3],
				 float_a *an, float_a *pn, float_a cn[3],
				 float_a volmina, float_a volminb,
				 float_a (*dr)[3], float_a *dR, 
				 float_a (*d2rR)[AGBNP_MAX_OVERLAP_LEVEL][3]);
 float_a swf_invbr(float_a beta, float_a *fp);
 int agbnp_getatomneighbors(AGBNPdata *agb, int iat, int *nbnum, 
			    int bsize, int *nbiat);
 int agbnp_getatomneighbors_pbc(AGBNPdata *agb, int iat, int *nbnum, 
				  int bsize, int *nbiat,
				  float_a *xj, float_a *yj, float_a *zj,
				  void **nbdata);
 int agbnp_neighbor_lists(AGBNPdata *agb, AGBworkdata *agbw,
				float_a *x, float_a *y, float_a *z);
 int agbnp_neighbor_lists_pbc(AGBNPdata *agb, AGBworkdata *agbw,
				float_a *x, float_a *y, float_a *z);
 int agbnp_self_volumes(AGBNPdata *agb, AGBworkdata *agbw,
			      float_a *x, float_a *y, float_a *z);
 int agbnp_self_volumes_pbc(AGBNPdata *agb, AGBworkdata *agbw,
			      float_a *x, float_a *y, float_a *z);
 int agbnp_scaling_factors(AGBNPdata *agb, AGBworkdata *agbw_h);
int agbnp_inverse_born_radii_pbc(AGBNPdata *agb, AGBworkdata *agbw_h,
				 float_a *x, float_a *y, float_a *z,
				 int init_frozen);
int agbnp_inverse_born_radii(AGBNPdata *agb, AGBworkdata *agbw_h,
			     float_a *x, float_a *y, float_a *z,
			     int init_frozen);
 int agbnp_born_radii(AGBNPdata *agb, AGBworkdata *agbw_h);
 int agbnp_reset_derivatives(AGBNPdata *agb, AGBworkdata *agbw_h);
 int agbnp_gb_energy_pbc(AGBNPdata *agb, AGBworkdata *agbw_h,
			   float_a *x, float_a *y, float_a *z,
			   float_a *egb_selfs, float_a *egb_pairs,
			   int init_frozen);
 int agbnp_gb_energy(AGBNPdata *agb, AGBworkdata *agbw_h,
			   float_a *x, float_a *y, float_a *z,
			   float_a *egb_selfs, float_a *egb_pairs,
			   int init_frozen);
int agbnp_gb_energy_inner_pbc(AGBNPdata *agb, int iat, int natoms, 
		int nb, int *nbiat, float_a *x, float_a *y, float_a *z,
		float_a qiat, float_a biat,
		float_a *charge, float_a *br, float_a *ab, float_a *dera,
		float_a (*dgbdr)[3], float_a *egb_pair, 
		float_a dielectric_factor,
		int initfrozen, float_a (*dgbdr_c)[3],
		float_a *gbpair_const, float_a *dera_const,
		int *isbfrozen,
		float_a *xj, float_a *yj, float_a *zj,
		float_a (*rot)[3][3]);
int agbnp_gb_energy_inner(AGBNPdata *agb, int iat, int natoms, 
		int nb, int *nbiat, float_a *x, float_a *y, float_a *z,
		float_a qiat, float_a biat,
		float_a *charge, float_a *br, float_a *ab, float_a *dera,
		float_a (*dgbdr)[3], float_a *egb_pair, 
		float_a dielectric_factor,
		int initfrozen, float_a (*dgbdr_c)[3],
		float_a *gbpair_const, float_a *dera_const,
		int *isbfrozen);
 int agbnp_gb_ders_constvp(AGBNPdata *agb, AGBworkdata *agbw_h,
				 float_a *x, float_a *y, float_a *z,
				 int init_frozen);
 int agbnp_gb_ders_constvp_pbc(AGBNPdata *agb, AGBworkdata *agbw_h,
				 float_a *x, float_a *y, float_a *z,
				 int init_frozen);
 int agbnp_der_vp(AGBNPdata *agb, AGBworkdata *agbw_h,
			float_a *x, float_a *y, float_a *z,
			int init_frozen);
 int agbnp_der_vp_pbc(AGBNPdata *agb, AGBworkdata *agbw_h,
			float_a *x, float_a *y, float_a *z,
			int init_frozen);
 int agbnp_gb_deruv(AGBNPdata *agb, AGBworkdata *agbw_h, 
			  int init_frozen);
 int agbnp_gb_deruv_pbc(AGBNPdata *agb, AGBworkdata *agbw_h, 
			  int init_frozen);
 int agb_vdw_energy(AGBNPdata *agb, float_i *x, float_i *y, float_i *z,
			  float_i *sp, float_i *br, float_i *egb, 
			  float_i (*dgbdr)[3],
			  float_i *evdw, float_i (*dvwdr)[3],
			  float_i *ecorr);
 int agbnp_neighbor_list_cavity(AGBNPdata *agb, CAVworkdata *cavw,
				      float_a *x, float_a *y, float_a *z);
 int agbnp_neighbor_list_cavity_pbc(AGBNPdata *agb, CAVworkdata *cavw,
				      float_a *x, float_a *y, float_a *z);
 int agbnp_surface_areas(AGBNPdata *agb, CAVworkdata *cavw,
			       float_a *x, float_a *y, float_a *z);
 int agbnp_surface_areas_pbc(AGBNPdata *agb, CAVworkdata *cavw,
			       float_a *x, float_a *y, float_a *z);
 int agbnp_cavity_ders(AGBNPdata *agb, CAVworkdata *cavw,
			     float_a *x, float_a *y, float_a *z);
 int agbnp_cavity_ders_pbc(AGBNPdata *agb, CAVworkdata *cavw,
			     float_a *x, float_a *y, float_a *z);
 int agb_cavity_energy(AGBNPdata *agb, float_i *x, float_i *y, float_i *z, 
		       float_i *mol_volume, float_i *surf_area,
		       float_i *ecav, float_i *ecorr, float_i (*decav)[3]);
 int agbnp_init_frozen_cavity(AGBNPdata *agb, 
			      float_i *x, float_i *y, float_i *z);
 int agbnp_set_vfrozen(AGBNPdata *agb, AGBworkdata *agbw, 
			     CAVworkdata *cavw);
 int agbnp_set_bfrozen(AGBNPdata *agb, AGBworkdata *agbw); 
 int agbnp_init_frozen_agb(AGBNPdata *agb, 
				 float_i *x, float_i *y, float_i *z);

 int agbnp_create_msatoms(AGBNPdata *agb);
 float_a agbnp_pol_switchfunc(float_a x, float_a xa, float_a xb,
			      float_a *fp, float_a *fpp);

 float_a altweight_pair(AGBNPdata *agb, int iat, int jat);
 float_a altweight_incremental(AGBNPdata *agb, int nat, 
			       int *iat, float_a iat_occupancy, int jat);
 int agbnp_get_neigh_coords(AGBNPdata *agb, 
				  int iat, NeighList *neigh_list, 
				  float_a *xj, float_a *yj, float_a *zj);
int mymax(int a, int b);
void rtvec(float_a y[3], float_a rot[3][3], float_a x[3]);


int agbnp_create_wsatoms(AGBNPdata *agb);
int agbnp_create_ws_ofatom(AGBNPdata *agb, int iat, int iws, int *nws);
int agbnp_create_ws_atoms_ph(AGBNPdata *agb, int iat, int iws, int *nws);
void agbnp_place_wat_hydrogen(float_a xd, float_a yd, float_a zd,
			      float_a xh, float_a yh, float_a zh, 
			      float_a d, 
			      float_a *xw, float_a *yw, float_a *zw);
int agbnp_ws_free_volumes(AGBNPdata *agb);
void agbnp_place_wat_lp_sp2(float_a xa, float_a ya, float_a za, 
			    float_a xr, float_a yr, float_a zr, 
			    float_a xr1, float_a yr1, float_a zr1, 
			    float_a xr2, float_a yr2, float_a zr2,
			    float_a d,
			    float_a *xw1, float_a *yw1, float_a *zw1,
			    float_a *xw2, float_a *yw2, float_a *zw2);
void agbnp_place_wat_lp_sp3_o(float_a xa, float_a ya, float_a za, 
			      float_a xr1, float_a yr1, float_a zr1, 
			      float_a xr2, float_a yr2, float_a zr2,
			      float_a d,
			      float_a *xw1, float_a *yw1, float_a *zw1,
			      float_a *xw2, float_a *yw2, float_a *zw2);
int agbnp_ws_correction(AGBNPdata *agb, float_a *ce);
#ifdef NOTNOW
 int agbnp_find_wparents(AGBNPdata *agb, int iat, int *jat, int *kat);
 int agbnp_calc_waxes(AGBNPdata *agb, int iat, int jat, int kat,
			   float_a origin[3], 
			   float_a x[3], float_a y[3], float_a z[3]);
 int agbnp_init_watoms(AGBNPdata *agb);
 int agbnp_position_watoms(AGBNPdata *agb);
#endif
void agbnp_cpyd(float_a *local, float_i *caller, int n);


#if defined(USE_SSE) && !defined(USE_PBC)
__m128 vexpf(__m128 vy);
void init_vexpf();
float pow2lt1(float x);
float expbf(float y);
#endif

/* functions to create and evaluate lookup tables for i4() */
#define F4LOOKUP_MAXA (20.0)
#define F4LOOKUP_SMOOTHA (10.0)
#define F4LOOKUP_NA (128)
#define F4LOOKUP_MAXB (8.0)
#define F4LOOKUP_NB (64)
int agbnp_fill_ctable(int n, float_a *x, float_a *y, float_a *yp,
		      C1Table *c1table);
int agbnp_create_ctablef4(int n, float_a amax, float_a b, 
			  C1Table **c1table);
int agbnp_create_ctablef42d(int na, float_a amax, 
			    int nb, float_a bmax, 
			    C1Table2D **table2d);
int agbnp_interpolate_ctable(C1Table *c1table, float_a x, 
			     float_a *f, float_a *fp);
int agbnp_interpolate_ctablef42d(C1Table2D *table2d, float_a x, float_a y,
				 float_a *f, float_a *fp);
int agbnp_init_i4p(void);

int plop_agb_cavity_energy(AGBNPdata *agb, int natom_asu, float_i *x, float_i *y, float_i *z,
                           float_i *mol_volume, float_i *surf_area,
                           float_i *ecav, float_i *ecorr);

int plop_agb_vdw_energy(AGBNPdata *agb, int natom_asu, float_i *x, float_i *y, float_i *z,
                        float_i *sp, float_i *br, float_i *evdw, float_i *ecorr);

int agbnp_fcompare( const void *val1, const void *val2 );
void agbnp_fsortindx( int pnval, float_a *val, int *indx );
int agbnp_nblist_reorder(AGBworkdata *agbw, NeighList *nl, int iat, int *indx);
int agbnp_nblist_reorder_cavity( CAVworkdata *cavw, NeighList *nl, int iat, int *indx);

#endif
