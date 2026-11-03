#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

#include "libnblist.h"
#include "agbnp.h"

#ifndef AGBNP_MAXSTRING
#define AGBNP_MAXSTRING 1024
#endif

int main(void){

  int natoms = 0;
  NeighList *neigh_list = NULL, *excl_neigh_list = NULL; /* ptrs to Verlet neighbor lists */
  int nheavyat = 0, *iheavyat = NULL; /* number of and list of heavy atoms */
  int nhydrogen = 0, *ihydrogen = NULL; /* number of and list of hydrogen atoms */
  int ndummy = 0, *idummy = NULL; /* number of and list of dummy atoms */
  int *isfrozen = NULL; /* 1 if atom is frozen, 0 otherwise */
  double *x = NULL, *y = NULL, *z = NULL; /* atom coordinates */
  double *r = NULL, *rcav = NULL;         /* atom radii */
  double *charge = NULL;    /* atom charges */
  double *sp = NULL;        /* scaled volume factors returned by agbnp() */
  double *br = NULL;        /* born radii returned by agbnp() */
  double mol_volume;
  double *surf_area = NULL; /* atoms vdw surface areas returned by agbnp */
  double (*dgbdr)[3] = NULL; /* positional derivatives of GB energy */
  double (*dvwdr)[3] = NULL; /* positional derivatives of NP energy */
  double (*decav)[3] = NULL; /* positional derivatives of cavity energy */
  double *gamma = NULL;      /* gamma np parameters */
  double *alpha = NULL;      /* alpha np parameters */
  double *delta = NULL;      /* delta np parameters */
  double *sgamma = NULL;      /* gamma np correction parameters */
  double *salpha = NULL;      /* alpha np correction parameters */
  double *sdelta = NULL;      /* delta np correction parameters */
  double *ab = NULL;          /* GB pair distance correction factor */
  int agbnp_tag = -1;
  int agbnpfrz_tag = -1;
  double dielectric_in = 1.0;
  double dielectric_out = 80.0;
  double ecav, ecorr_cav;
  double egb, evdw, ecorr_vdw;
  double agb, np, ecorr, egnp;
  int dopbc = 0;
  NeighList *conntbl = NULL;
  int i, iat;

  /* read number of atoms */
  if(scanf("%d",&natoms) != 1){
    fprintf(stderr,"agbnp_test(): error reading number of atoms from stdin.\n");
    return 1;
  }

  /* allocate parameter arrays */
  x = (double *)realloc(x, natoms*sizeof(double));
  y = (double *)realloc(y, natoms*sizeof(double));
  z = (double *)realloc(z, natoms*sizeof(double));
  r = (double *)realloc(r, natoms*sizeof(double));
  rcav = (double *)realloc(rcav, natoms*sizeof(double));
  charge = (double *)realloc(charge, natoms*sizeof(double));
  iheavyat = (int *)realloc(iheavyat, natoms*sizeof(int));
  ihydrogen = (int *)realloc(ihydrogen, natoms*sizeof(int));
  idummy = (int *)realloc(idummy, natoms*sizeof(int));
  sp = (double *)realloc(sp, natoms*sizeof(double));
  br = (double *)realloc(br, natoms*sizeof(double));
  surf_area = (double *)realloc(surf_area, natoms*sizeof(double));
  dgbdr = (double (*)[3])realloc(dgbdr, natoms*sizeof(double [3]));
  gamma = (double *)realloc(gamma, natoms*sizeof(double));
  alpha = (double *)realloc(alpha, natoms*sizeof(double));
  delta = (double *)realloc(delta, natoms*sizeof(double));
  sgamma = (double *)realloc(sgamma, natoms*sizeof(double));
  salpha = (double *)realloc(salpha, natoms*sizeof(double));
  sdelta = (double *)realloc(sdelta, natoms*sizeof(double));
  ab = (double *)realloc(ab, natoms*sizeof(double));
  dvwdr = (double (*)[3])realloc(dvwdr, natoms*sizeof(double [3]));
  decav = (double (*)[3])realloc(decav, natoms*sizeof(double [3]));
  if(!(x && y && z && r && charge && iheavyat &&
       ihydrogen && idummy && sp && br && surf_area &&
       dgbdr && gamma && alpha && delta && sgamma &&
       salpha && sdelta && ab && dvwdr && decav)){
    fprintf(stderr,"agbnp_test(): unable to allocate memory for one or more parameter arrays.\n");
    return 1;
  }

  for(iat=0;iat<natoms;iat++){
    if(scanf("%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &i, 
	     &(x[iat]), &(y[iat]), &(z[iat]),
	     &(r[iat]), &(rcav[iat]),
	     &(charge[iat]),
	     &(gamma[iat]), &(sgamma[iat]),
	     &(alpha[iat]), &(salpha[iat]),
	     &(delta[iat]), &(sdelta[iat]),
	     &(ab[iat])) != 14){
      fprintf(stderr,"agbnp_test(): error reading parameters from stdin.\n");
      return 1;
    }
  }
  if(scanf("%d",&nhydrogen) != 1){
      fprintf(stderr,"agbnp_test(): error reading nhydrogen from stdin.\n");
      return 1;
  }
  for(i=0;i<nhydrogen;i++){
    if(scanf("%d",&(ihydrogen[i])) != 1){
      fprintf(stderr,"agbnp_test(): error reading ihydrogen from stdin.\n");
      return 1;
    }
  }

  /* initialize libagbnp */
  if(agbnp_initialize() != AGBNP_OK){
    fprintf(stderr, "agbnp_test: error in agbnp_initialize().\n");
    return 1;
  }

  /* create a libagbnp object */
  if(agbnp_new(&(agbnp_tag), natoms, 
	       x, y, z, r, 
	       charge, dielectric_in, dielectric_out,
	       gamma, sgamma,
	       alpha, salpha,
	       delta, sdelta,
	       ab,
	       nhydrogen, ihydrogen, 
	       ndummy, idummy, isfrozen,
	       neigh_list, excl_neigh_list,
	       dopbc, 0, 0, NULL, NULL, NULL, NULL, NULL) != AGBNP_OK){
    fprintf(stderr, "agbnp_test: error in agbnp_new().\n");
    return 1;
  }

  /* initialize frozen atoms constant terms */
  if(agbnp_init_frozen(agbnp_tag, x, y, z, isfrozen) 
     != AGBNP_OK){
    fprintf(stderr, "agbnp_test: error in agbnp_init_frozen().\n");
    return 1;
  }

  /* compute cavity energy */
  if(agbnp_cavity_energy(agbnp_tag, 
			 x, y, z,
			 &(mol_volume), surf_area, 
			 &ecav, &ecorr_cav, 
			 decav ) != AGBNP_OK){
    fprintf(stderr, "agbnp_test: error in agbnp_cavity_energy().\n");
    return 1;
  }
#ifdef NOTNOW
  ecav = 0.0;
  ecorr_cav = 0.0;
  memset(decav,0,3*agbnp_w->natoms*sizeof(double));
#endif

  /* calls agbnp() */
  if(agbnp_agb_energy(agbnp_tag, x, y, z,
		      sp,  br, &egb, dgbdr,
		      &evdw, dvwdr, &ecorr_vdw) != AGBNP_OK){
    fprintf(stderr, "agbnp_test: error in agbnp_agb_energy().\n");
    return 1;
  }
#ifdef NOTNOW
  egb = 0.0;
  evdw = 0.0;
  ecorr_vdw = 0.0;
  memset(dgbdr,0,3*natoms*sizeof(double));
  memset(dvwdr,0,3*natoms*sizeof(double));
#endif

  printf("AGB: %lf\n", egb);
  printf("CAV: %lf\n", ecav);
  printf("VDW: %lf\n", evdw);
  printf("NP : %lf\n", ecav+evdw);
  printf("CAVcorr: %lf\n", ecorr_cav);
  printf("VDWcorr: %lf\n", ecorr_vdw);
  printf("CORR: %lf\n", ecorr_cav + ecorr_vdw);
  printf("AGBNP+CORR: %lf\n", egb+ecav+evdw+ecorr_cav+ecorr_vdw);


  /* now do the same thing with frozen atoms */
  isfrozen = (int *)calloc(natoms,sizeof(int));
  if(!isfrozen){
    fprintf(stderr, "agbnp_test: unable to allocate isfrozen array.\n");
    return 1;
  }
  for(iat=111;iat<191;iat++){ /* freezes residues 48 thru 53 */ 
    isfrozen[iat] = 1;
  }
  if(agbnp_new(&(agbnpfrz_tag), natoms, 
	       x, y, z, r, 
	       charge, dielectric_in, dielectric_out,
	       gamma, sgamma,
	       alpha, salpha,
	       delta, sdelta,
	       ab,
	       nhydrogen, ihydrogen, 
	       ndummy, idummy, isfrozen,
	       neigh_list, excl_neigh_list,
	       dopbc, 0, 0, NULL, NULL, NULL, NULL, NULL) != AGBNP_OK){
    fprintf(stderr, "agbnp_test: error in agbnp_new().\n");
    return 1;
  }
  

  /* initialize frozen atoms constant terms */
  if(agbnp_init_frozen(agbnp_tag, x, y, z, isfrozen) 
     != AGBNP_OK){
    fprintf(stderr, "agbnp_test: error in agbnp_init_frozen().\n");
    return 1;
  }

  /* compute cavity energy */
  if(agbnp_cavity_energy(agbnp_tag, 
			 x, y, z,
			 &(mol_volume), surf_area, 
			 &ecav, &ecorr_cav, 
			 decav ) != AGBNP_OK){
    fprintf(stderr, "agbnp_test: error in agbnp_cavity_energy().\n");
    return 1;
  }
#ifdef NOTNOW
  ecav = 0.0;
  ecorr_cav = 0.0;
  memset(decav,0,3*agbnp_w->natoms*sizeof(double));
#endif

  /* calls agbnp() */
  if(agbnp_agb_energy(agbnp_tag, x, y, z,
		      sp,  br, &egb, dgbdr,
		      &evdw, dvwdr, &ecorr_vdw) != AGBNP_OK){
    fprintf(stderr, "agbnp_test: error in agbnp_agb_energy().\n");
    return 1;
  }
#ifdef NOTNOW
  egb = 0.0;
  evdw = 0.0;
  ecorr_vdw = 0.0;
  memset(dgbdr,0,3*natoms*sizeof(double));
  memset(dvwdr,0,3*natoms*sizeof(double));
#endif

  printf("\nResults with Frozen Atoms (should be the same as before)\n");

  printf("AGB: %lf\n", egb);
  printf("CAV: %lf\n", ecav);
  printf("VDW: %lf\n", evdw);
  printf("NP : %lf\n", ecav+evdw);
  printf("CAVcorr: %lf\n", ecorr_cav);
  printf("VDWcorr: %lf\n", ecorr_vdw);
  printf("CORR: %lf\n", ecorr_cav + ecorr_vdw);
  printf("AGBNP+CORR: %lf\n", egb+ecav+evdw+ecorr_cav+ecorr_vdw);


  if(agbnp_delete(agbnp_tag) != AGBNP_OK){
    fprintf(stderr, "agbnp_test: error in agbnp_delete().\n");
    return 1;
  }

  if(agbnp_delete(agbnpfrz_tag) != AGBNP_OK){
    fprintf(stderr, "agbnp_test: error in agbnp_delete().\n");
    return 1;
  }

  agbnp_terminate();

  return 0;
}


