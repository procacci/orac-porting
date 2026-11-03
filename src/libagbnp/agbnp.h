#ifndef AGBNP_H
#define AGBNP_H

#include "libnblist.h"

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif

#ifndef FLOAT_I_H
#define FLOAT_I_H
#ifdef SINGLE_PREC_INTERFACE
typedef float float_i;
#else
typedef double float_i;
#endif
#endif

#define AGBNP_OK (2)
#define AGBNP_ERR (-1)

/* Initializes libagbnp library.*/
int agbnp_initialize( void );

/* Terminate libagbnp library. */
void agbnp_terminate( void );

/* creates a new public instance of an agbnp structure */
int agbnp_new(int *tag, int natoms, 
	      float_i *x, float_i *y, float_i *z, float_i *r, 
	      float_i *charge, float_i dielectric_in, float_i dielectric_out,
	      float_i *igamma, float_i *sgamma,
	      float_i *ialpha, float_i *salpha,
	      float_i *idelta, float_i *sdelta,
	      float_i *ab, 
#ifdef AGBNP_HB
	      int *hbtype, float_i *hbcorr,
#endif
	      int nhydrogen, int *ihydrogen, 
	      int ndummy, int *idummy,
	      int *isfrozen,
	      NeighList *neigh_list, NeighList *excl_neigh_list,
	      int dopbc, int nsym, int ssize,
              float_i *xs, float_i *ys, float_i *zs,
	      float_i (*rot)[3][3], NeighList *conntbl);

/* deletes a AGBNP object */
int agbnp_delete(int tag);

/* generalized born and van der Waals energies */
int agbnp_agb_energy(int tag, float_i *x, float_i *y, float_i *z,
		     float_i *sp, float_i *br, float_i *egb, float_i (*dgbdr)[3],
		     float_i *evdw, float_i (*dvwdr)[3], float_i *ecorr);

/* cavity energy */
int agbnp_cavity_energy(int tag, float_i *x, float_i *y, float_i *z,
			float_i *mol_volume, float_i *surf_area, 
			float_i *ecav, float_i *ecorr,
			float_i (*decav)[3]);

/* return born radii and scaled radii only (no energies, no derivatives) */
int agbnp_bornr(int tag, float_i *x, float_i *y, float_i *z,
		float_i *sp, float_i *br);


/* (re)-initialize constant terms from frozen atoms */
int agbnp_init_frozen(int tag, float_i *x, float_i *y, float_i *z, int *isfrozen);

/* set alternate conformation information */
int agbnp_set_altconfs(int tag, float_i *occupancy, int *alt_site, int *alt_id);

int plop_agbnp_cavity_energy(int tag, int natom_asu, float_i *x, float_i *y, float_i *z,
                             float_i *mol_volume, float_i *surf_area,
                             float_i *ecav, float_i *ecorr);

int plop_agbnp_agb_energy(int tag, int natom_asu, float_i *x, float_i *y, float_i *z,
                          float_i *sp, float_i *br, float_i *evdw, float_i *ecorr);

#endif
