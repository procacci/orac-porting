#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_SSE
#include <xmmintrin.h>
#endif

#include "agbnp.h"
#include "libnblist.h"
#include "agbnp_private.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef MMTIMER
#include "mmtimer.h"
static int cavity_timer = -1;
static int agbvdw_timer = -1;
#endif


#ifndef USE_PBC

#ifdef USE_SSE
static const float ln2 = 0.693147181;
static const float ln2inv = 1.442695041;
static const float twoexpa[8] = {
  1.0000000000,
  0.6931471810,
  0.2402265073,
  0.0555041088,
  0.0096181291,
  0.0013333558,
  0.0001540353,
  0.0000152527
};
static __m128 vln2inv;
static __m128 vtwoexpa[8];
void init_vexpf(){
  int i;
  vln2inv = _mm_set_ps(ln2inv,ln2inv,ln2inv,ln2inv);
  for(i=0;i<7;i++){
    vtwoexpa[i] = _mm_set_ps(twoexpa[i],twoexpa[i],twoexpa[i],twoexpa[i]);
  }
}
#endif


/*                                                     *
 * Data structures that holds list of AGBNP structures *
 *                                                     */
static int agbnp_initialized = FALSE;
/* a pointer to all the agbnp strctures which have been allocated */
static AGBNPdata *agbdata_list = NULL;
/* the number of allocated agbnp structures */
int agbdata_allocated = 0;
/* the number of atom lists currently in use */
int agbdata_used = 0;
/* the initial number of agbnp structures allocated */
static const int AGBDATA_INITIAL_NUMBER = 1;
/* the increment of agbnp structures slots for reallocation */
static const int AGBDATA_INCREMENT = 1;


/* lookup table for i4/u function */
static C1Table2D *agbnp_f4c1table2d = NULL;

/* Initializes libagbnp library.*/
int agbnp_initialize( void ){
  int i;

  if(agbnp_initialized){
    return AGBNP_OK;
  }

  /* creates list of pointers of agbnp structures managed by library */
  agbdata_list = (AGBNPdata * ) 
    calloc(AGBDATA_INITIAL_NUMBER, sizeof(AGBNPdata) );
  if(!agbdata_list) {
    fprintf(stderr,"agbnp_initialize(): error allocating memory for %d agbnp objects.\n", AGBDATA_INITIAL_NUMBER);
     return AGBNP_ERR;
  }
  agbdata_allocated = AGBDATA_INITIAL_NUMBER;
  agbdata_used = 0; 

  /* reset lists values */
  for(i=0; i<agbdata_allocated; i++){
    /* set in_use=0 among other things */
    agbnp_reset(&(agbdata_list[i]));
  }

  agbnp_initialized = TRUE;
  
  return AGBNP_OK;
}

/* Terminate libagbnp library. */
void agbnp_terminate( void ){
  int i;

  if(!agbnp_initialized){
    fprintf(stderr,"agbnp_terminate(): agbnp library is not initialized.\n");
    return;
  }

  /* deletes allocated structures */
  for(i = 0; i<agbdata_allocated ; i++){
    if(agbdata_list[i].in_use){
      if(agbnp_delete(i) != AGBNP_OK){
	fprintf(stderr, "agbnp_terminate(): error in deleting agbnp structure %d.\n", i);
      }
    }
  }

  /* free list of allocated lists */
  free(agbdata_list);
  agbdata_list = NULL;
  agbdata_allocated = 0;
  agbdata_used = 0;
}

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
	      float_i (*rot)[3][3], NeighList *conntbl){

  int do_ms = 0;

  int slot, il, i, iat, indx;
  int *iswhat;
  AGBNPdata *agbdata;
#ifdef _OPENMP
  int iproc;
#endif

  if(!agbnp_initialized){
    fprintf(stderr,"agbnp_new(): agbnp library is not initialized.\n");
    return AGBNP_ERR;
  }

  /* find an allocated structure not in use */
  slot = 0;
  while(agbdata_list[slot].in_use){

    if(slot+1 >= agbdata_allocated){
      /* can't find it, expand list */
      agbdata_list = (AGBNPdata * ) 
	realloc(agbdata_list, 
	   (agbdata_allocated+AGBDATA_INCREMENT)*sizeof(AGBNPdata) );
      if(!agbdata_list){
	fprintf(stderr,"agbnp_new(): error reallocating memory for %d agbnp objects.\n", agbdata_allocated+AGBDATA_INCREMENT);
	return AGBNP_ERR;
      }
      /* reset new lists */
      for(il=0; il < AGBDATA_INCREMENT ; il++){
	agbnp_reset(&(agbdata_list[agbdata_allocated+il]));
      }
      /* increase number of allocated lists */
      agbdata_allocated += AGBDATA_INCREMENT;
    }

    /* next slot */
    slot += 1;
  }

  /* this is the index of the available structure */
  *tag = slot;
  agbdata = &(agbdata_list[slot]);

  /* reset new structure */
  agbnp_reset(agbdata);

  /* set natoms */
  agbdata->natoms = natoms;

  /* allocates and set coordinates */
  agbdata->x = (float_a *)calloc(natoms,sizeof(float_a));
  agbdata->y = (float_a *)calloc(natoms,sizeof(float_a));
  agbdata->z = (float_a *)calloc(natoms,sizeof(float_a));
  if(!(agbdata->x && agbdata->y && agbdata->z)){
    fprintf(stderr,"agbnp_new(): error allocating memory for coordinates (%d doubles)\n",3*natoms);
    return AGBNP_ERR;
  }
  for(iat = 0; iat < natoms; iat++){
    agbdata->x[iat] = x[iat];
    agbdata->y[iat] = y[iat];
    agbdata->z[iat] = z[iat];
  }

  /* allocates and set radii */
  agbdata->r = (float_a *)calloc(natoms,sizeof(float_a));
  agbdata->rcav = (float_a *)calloc(natoms,sizeof(float_a));
  if(!(agbdata->r && agbdata->rcav)){
     fprintf(stderr,"agbnp_new(): error allocating memory for atomic radii (%d doubles)\n",2*natoms);
    return AGBNP_ERR;
  }
  for(iat = 0; iat < natoms; iat++){
    agbdata->r[iat] = r[iat];
    agbdata->rcav[iat] = r[iat] + AGBNP_RADIUS_INCREMENT;
  }

  /* allocates and set charges */
  agbdata->charge = (float_a *)calloc(natoms,sizeof(float_a));
  if(!agbdata->charge){
     fprintf(stderr,"agbnp_new(): error allocating memory for atomic partial charges (%d doubles)\n",natoms);
    return AGBNP_ERR;
  }
  for(iat = 0; iat < natoms; iat++){
    agbdata->charge[iat] = charge[iat];
  }

  /* allocates and set np parameters */
  agbdata->igamma = (float_a *)calloc(natoms,sizeof(float_a));
  agbdata->sgamma = (float_a *)calloc(natoms,sizeof(float_a));
  agbdata->ialpha = (float_a *)calloc(natoms,sizeof(float_a));
  agbdata->salpha = (float_a *)calloc(natoms,sizeof(float_a));
  agbdata->idelta = (float_a *)calloc(natoms,sizeof(float_a));
  agbdata->sdelta = (float_a *)calloc(natoms,sizeof(float_a));
  agbdata->ab = (float_a *)calloc(natoms,sizeof(float_a));
  if(!(agbdata->igamma && agbdata->sgamma && agbdata->ialpha &&
       agbdata->salpha && agbdata->idelta && agbdata->sdelta &&
       agbdata->ab)){
    fprintf(stderr,"agbnp_new(): error allocating memory for non-polar parameters (%d doubles)\n",6*natoms);
    return AGBNP_ERR;
  }
#ifdef AGBNP_HB
  agbdata->hbtype = (int *)calloc(natoms,sizeof(int));
  agbdata->hbcorr = (float_a *)calloc(natoms,sizeof(float_a));
  if(!(agbdata->hbtype && agbdata->hbcorr)){
    fprintf(stderr,"agbnp_new(): error allocating memory for non-polar parameters\n");
    return AGBNP_ERR;
  }
#endif
  for(iat=0;iat<natoms;iat++){
    agbdata->igamma[iat] = igamma[iat];
    agbdata->sgamma[iat] = sgamma[iat];
    agbdata->ialpha[iat] = ialpha[iat];
    agbdata->salpha[iat] = salpha[iat];
    agbdata->idelta[iat] = idelta[iat];
    agbdata->sdelta[iat] = sdelta[iat];
    agbdata->ab[iat]     = ab[iat];
#ifdef AGBNP_HB
    agbdata->hbtype[iat] = hbtype[iat];
    agbdata->hbcorr[iat] = hbcorr[iat];
#endif
  }

  /* allocates and set list of hydrogens, heavy atoms and dummy atoms */
  iswhat = (int *)calloc(natoms,sizeof(int)); /* 1 = heavy atom, 
						 2 = hydrogen
						 3 = dummy */
  agbdata->iheavyat = (int *)calloc(natoms,sizeof(int));
  agbdata->ihydrogen = (int *)calloc(natoms,sizeof(int));
  agbdata->idummy = (int *)calloc(natoms,sizeof(int));
  if(!(iswhat && agbdata->iheavyat && agbdata->ihydrogen && agbdata->idummy)){
    fprintf(stderr,"agbnp_new(): error allocating memory for lists of hydrogens, heavy atoms and dummy atoms (%d integers)\n",4*natoms);
    return AGBNP_ERR;
  }

  for(i=0;i<nhydrogen;i++){
    iat = ihydrogen[i];
    iswhat[iat] = 2;
    agbdata->ihydrogen[i] = iat;
  }
  agbdata->nhydrogen = nhydrogen;

  for(i=0;i<ndummy;i++){
    iat = idummy[i];
    if(iswhat[iat]){
      fprintf(stderr,"agbnp_new(): hydrogen atom %d cannot be also set as a dummy atom.\n",iat);
      return AGBNP_ERR;
    }
    iswhat[iat] = 3;
    agbdata->idummy[i] = iat;
  }
  agbdata->ndummy = ndummy;

  /* an atom that is not an hydrogen atom or a dummy atom is a heavy atom */
  agbdata->nheavyat = 0;
  for(iat=0;iat<natoms;iat++){
    if(!iswhat[iat]){
      agbdata->iheavyat[agbdata->nheavyat++] = iat;
    }
  }
  free(iswhat);

  /* set dielectric contstants */
  agbdata->dielectric_in = dielectric_in;
  agbdata->dielectric_out = dielectric_out;

  /* set frozen atom flag */
  if(isfrozen){
    agbdata->do_frozen = 1;
    agbdata->isfrozen = (int *)calloc(natoms,sizeof(int));
    if(!agbdata->isfrozen){
      fprintf(stderr,"agbnp_new(): unable to allocate isfrozen array (%d ints)\n",natoms);
      return AGBNP_ERR;
    }
    for(iat=0;iat<natoms;iat++){
      agbdata->isfrozen[iat] = isfrozen[iat];
    }
  }

  /* allocate and set default occupancies */
  agbdata->occupancy = (float_a *)calloc(natoms,sizeof(float_a));
  if(!agbdata->occupancy){
    fprintf(stderr,"agbnp_new(): unable to allocate occupancy array (%d doubles)\n",natoms);
      return AGBNP_ERR;
  }
  for(iat=0;iat<natoms;iat++){
    agbdata->occupancy[iat] = 1.0;
  }

  /* set pointers to Verlet neighbor lists if given */
  if(neigh_list){
    agbdata->neigh_list = neigh_list;
    if(excl_neigh_list){
      agbdata->excl_neigh_list = excl_neigh_list;
    }
  }

  /* PBC stuff */
  agbdata->nsym = 1;
  if(dopbc){
    int i,j,isym;
    if(!neigh_list){
      fprintf(stderr,"agbnp_new(): A neighbor list is needed with periodic boundary conditions.\n");
      return AGBNP_ERR;
    }
    agbdata->dopbc = 1;
    agbdata->nsym = nsym;
    if(nsym>1){
      /*
      if(!neigh_list->idx_remap){
	fprintf(stderr,"agbnp_new(): The neighbor list index remapping array is needed with crystal periodic boundary conditions.\n");
	return AGBNP_ERR;
      }
      */
      agbdata->docryst = 1;
      /* assign crystal cell coordinate buffer */ 
      if(ssize < nsym*agbdata->natoms){
	fprintf(stderr,"agbnp_new(): spcified coordinate buffer size (%d) is smaller than number of atoms in unit cell (%d).\n",ssize, nsym*agbdata->natoms);
	return AGBNP_ERR;
      }
      if(!(xs && ys && zs)){
	fprintf(stderr,"agbnp_new(): one of more coordinate buffers pointers needed for crystal PBC's are NULL.\n");
	return AGBNP_ERR;
      }
      agbdata->xs = xs;
      agbdata->ys = ys;
      agbdata->zs = zs;

      /* allocate space and assign rotation matrixes */ 
      agbdata->rot = (float_a (*)[3][3])calloc(nsym,sizeof(float_a [3][3]));
      if(!agbdata->rot){
	fprintf(stderr,"agbnp_new(): unable to allocate array of cell rotation matrixes.\n");
        return AGBNP_ERR;
      }
      for(isym=0;isym<nsym;isym++){
	for(i=0;i<3;i++){
	  for(j=0;j<3;j++){
	    agbdata->rot[isym][i][j] = rot[isym][i][j];
	  }
	}
      }

    }
  }


  /* set for molecular surface */
  agbdata->do_ms = do_ms;

  /* copy connection table to libagbnp internal structure */
  if(conntbl){
#ifdef AGBNP_HB
    agbdata->do_w = 1;
#endif
    agbdata->conntbl = (NeighList *)calloc(1,sizeof(NeighList));
    nblist_reset_neighbor_list(agbdata->conntbl);
    if(nblist_reallocate_neighbor_list(agbdata->conntbl,conntbl->natoms,
				       conntbl->neighl_size) != NBLIST_OK){
      fprintf(stderr, "agbnp_new(): unable to allocate connection table (size = %d ints).\n",conntbl->neighl_size );
      return AGBNP_ERR;
    }
    /* copy neighbor list */
    memcpy(agbdata->conntbl->nne, conntbl->nne, (conntbl->natoms)*sizeof(int));
    memcpy(agbdata->conntbl->neighl1, conntbl->neighl1, 
	   (conntbl->neighl_size)*sizeof(int));
    /* rebuilds array of pointers */
    indx = 0;
    for(iat=0;iat<natoms;iat++){
      agbdata->conntbl->neighl[iat] = &(agbdata->conntbl->neighl1[indx]);
      indx += agbdata->conntbl->nne[iat];
    }
  }

#ifdef AGBNP_WRITE_DATA
  printf("%d\n",natoms);
  /*  printf("i x y z  r rcav igamma sgamma ialpha salpha idelta sdelta ab\n"); */
  for(iat=0;iat<natoms;iat++){
#ifdef AGBNP_HB
    printf("%d %f %f %f %f %f %f %f %f %f %f %f %f %f %d %f\n",iat,
	   agbdata->x[iat],agbdata->y[iat],agbdata->z[iat],
	   agbdata->r[iat],agbdata->rcav[iat],
	   agbdata->charge[iat], 
	   agbdata->igamma[iat],agbdata->sgamma[iat],
	   agbdata->ialpha[iat], agbdata->salpha[iat], 
	   agbdata->idelta[iat] , agbdata->sdelta[iat],  
	   agbdata->ab[iat], agbdata->hbtype[iat], agbdata->hbcorr[iat]);
#else
    printf("%d %f %f %f %f %f %f %f %f %f %f %f %f %f\n",iat,
	   agbdata->x[iat],agbdata->y[iat],agbdata->z[iat],
	   agbdata->r[iat],agbdata->rcav[iat],
	   agbdata->charge[iat], 
	   agbdata->igamma[iat],agbdata->sgamma[iat],
	   agbdata->ialpha[iat], agbdata->salpha[iat], 
	   agbdata->idelta[iat] , agbdata->sdelta[iat],  
	   agbdata->ab[iat]);
#endif
  }
  printf("%d\n",nhydrogen);
  for(i=0;i<nhydrogen;i++){
    printf("%d\n",ihydrogen[i]);
  }
#endif

  /* allocates work arrays */
  agbdata->agbw = (AGBworkdata *)calloc(1,sizeof(AGBworkdata));
  if(!agbdata->agbw){
    fprintf(stderr,"agbnp_new(): error allocating memory for AGB work data structure.\n");
    return AGBNP_ERR;
  }
  agbnp_reset_agbworkdata(agbdata->agbw);
  if(agbnp_allocate_agbworkdata(natoms, agbdata, agbdata->agbw) != AGBNP_OK){
    fprintf(stderr,"agbnp_new(): error in agbnp_allocate_agbworkdata()\n");
    return AGBNP_ERR;
  }
  agbnp_init_agbworkdata(agbdata,agbdata->agbw);

  agbdata->cavw = (CAVworkdata *)calloc(1,sizeof(CAVworkdata));
  if(!agbdata->cavw){
    fprintf(stderr,"agbnp_new(): error allocating memory for Cavity work data structure.\n");
    return AGBNP_ERR;
  }
  agbnp_reset_cavworkdata(agbdata->cavw);
  if(agbnp_allocate_cavworkdata(natoms, agbdata, agbdata->cavw) != AGBNP_OK){
    fprintf(stderr,"agbnp_new(): error in agbnp_allocate_cavworkdata()\n");
    return AGBNP_ERR;
  }
  agbnp_init_cavworkdata(agbdata,agbdata->cavw);

#ifdef _OPENMP
  /* allocates work arrays for each thread */

#pragma omp parallel
  {
#pragma omp single
    {
      agbdata->nprocs = omp_get_num_threads();
    }
  }

  fprintf(stderr,"\n agbnp_new(): info: using %5d OpenMP thread(s).\n\n",
	 agbdata->nprocs);
  agbdata->maxprocs = omp_get_max_threads();
  agbdata->agbw_p = (AGBworkdata **)calloc(agbdata->maxprocs,sizeof(AGBworkdata *));
  agbdata->cavw_p = (CAVworkdata **)calloc(agbdata->maxprocs,sizeof(CAVworkdata *));
  if(!(agbdata->agbw_p && agbdata->cavw_p)){
    fprintf(stderr,"agbnp_new(): error allocating memory for array of pointers to thread work data structures.\n");
    return AGBNP_ERR;
  }
  for(iproc=0;iproc<agbdata->maxprocs;iproc++){
    agbdata->agbw_p[iproc] = (AGBworkdata *)calloc(1,sizeof(AGBworkdata));
    if(!agbdata->agbw_p[iproc]){
      fprintf(stderr,"agbnp_new(): error allocating memory for AGB work data structure.\n");
      return AGBNP_ERR;
    }
    agbnp_reset_agbworkdata(agbdata->agbw_p[iproc]);
    if(agbnp_allocate_agbworkdata(natoms,agbdata,agbdata->agbw_p[iproc]) != AGBNP_OK){
      fprintf(stderr,"agbnp_new(): error in agbnp_allocate_agbworkdata()\n");
      return AGBNP_ERR;
    }
    agbnp_init_agbworkdata(agbdata,agbdata->agbw_p[iproc]);
    
    agbdata->cavw_p[iproc] = (CAVworkdata *)calloc(1,sizeof(CAVworkdata));
    if(!agbdata->cavw_p[iproc]){
      fprintf(stderr,"agbnp_new(): error allocating memory for Cavity work data structure.\n");
      return AGBNP_ERR;
    }
    agbnp_reset_cavworkdata(agbdata->cavw_p[iproc]);
    if(agbnp_allocate_cavworkdata(natoms,agbdata,agbdata->cavw_p[iproc]) != AGBNP_OK){
      fprintf(stderr,"agbnp_new(): error in agbnp_allocate_cavworkdata()\n");
      return AGBNP_ERR;
    }
    agbnp_init_cavworkdata(agbdata,agbdata->cavw_p[iproc]);
  }

#endif /* _OPENMP */

#ifdef AGBNP_HB
  /* create water sites */
  if(agbdata->do_w){
    if(agbnp_create_wsatoms(agbdata) != AGBNP_OK){
      fprintf(stderr,"agbnp_new(): error in agbnp_create_wsatoms()\n");
      return AGBNP_ERR;
    }
  }
#endif

  /* create molecular surface atoms */
  if(agbdata->do_ms){
    if(agbnp_create_msatoms(agbdata) != AGBNP_OK){
      fprintf(stderr,"agbnp_new(): error in agbnp_create_msatoms()\n");
      return AGBNP_ERR;
    }
  }

#ifdef USE_SSE
  /* initializes lookup table routine for vectorized exp() */
  init_vexpf();
#endif

  /* initializes lookup table version of i4 */
  if(agbnp_init_i4p() != AGBNP_OK){
    fprintf(stderr,"agbnp_new(): error in agbnp_init_i4p()\n");
    return AGBNP_ERR;
  }

  /* set in_use=TRUE */
  agbdata->in_use = TRUE;

  /* update number of active structures */
  agbdata_used += 1;

  return AGBNP_OK;
}

/* deletes a public instance of an agbnp structure */
int agbnp_delete(int tag){
    AGBNPdata *agb;

  if(!agbnp_initialized){
    fprintf(stderr,"agbnp_delete(): agbnp library is not initialized.\n");
    return AGBNP_ERR;
  }
  if(!agbnp_tag_ok(tag)){
    fprintf(stderr,"agbnp_delete(): invalid tag %d.\n",tag);
    return AGBNP_ERR;
  }
  /* pointer to agb data structure */
  agb = &(agbdata_list[tag]);

  if(agb->x){ free(agb->x); agb->x = NULL;}
  if(agb->y){ free(agb->y); agb->y = NULL;}
  if(agb->z){ free(agb->z); agb->z = NULL;}
  if(agb->r){ free(agb->r); agb->r = NULL;}
  if(agb->rcav){ free(agb->rcav); agb->rcav = NULL;}
  if(agb->charge){ free(agb->charge); agb->charge = NULL;}
  if(agb->igamma){ free(agb->igamma); agb->igamma = NULL;}
  if(agb->sgamma){ free(agb->sgamma); agb->sgamma = NULL;}
  if(agb->ialpha){ free(agb->ialpha); agb->ialpha = NULL;}
  if(agb->salpha){ free(agb->salpha); agb->salpha = NULL;}
  if(agb->idelta){ free(agb->idelta); agb->idelta = NULL;}
  if(agb->sdelta){ free(agb->sdelta); agb->sdelta = NULL;}
  if(agb->ab){ free(agb->ab); agb->ab = NULL;}
#ifdef AGBNP_HB
  if(agb->hbtype){ free(agb->hbtype); agb->hbtype = NULL;}
  if(agb->hbcorr){ free(agb->hbcorr); agb->hbcorr = NULL;}
#endif
  if(agb->iheavyat){ free(agb->iheavyat); agb->iheavyat = NULL;}
  if(agb->ihydrogen){ free(agb->ihydrogen); agb->ihydrogen = NULL;}
  if(agb->idummy){ free(agb->idummy); agb->idummy = NULL;}
  if(agb->isfrozen){ free(agb->isfrozen) ; agb->isfrozen = NULL;}
  if(agb->rot){ free(agb->rot) ; agb->rot = NULL; }

  if(agb->agbw){
    agbnp_delete_agbworkdata(agb->agbw); free(agb->agbw) ; agb->agbw = NULL;}
  if(agb->cavw){
        agbnp_delete_cavworkdata(agb->cavw); free(agb->cavw) ; agb->cavw = NULL;}
  if(agb->conntbl){
    nblist_delete_neighbor_list(agb->conntbl); free(agb->conntbl) ; agb->conntbl = NULL; }

#ifdef AGBNP_HB
  if(agb->wsat){ free(agb->wsat); agb->wsat = NULL;}
  if(agb->wsat_list){
    nblist_delete_neighbor_list(agb->wsat_list); free(agb->wsat_list) ; agb->wsat_list = NULL; }
#endif

  if(agb->msat){ free(agb->msat); agb->msat = NULL;}
  if(agb->msat_list){
    nblist_delete_neighbor_list(agb->msat_list); free(agb->msat_list) ; agb->msat_list = NULL; }

#ifdef _OPENMP
  if(agb->maxprocs > 0){
    int iproc;
    if(agb->agbw_p){
      for(iproc=0;iproc<agb->maxprocs;iproc++){
	if(agb->agbw_p[iproc]){
	  agbnp_delete_agbworkdata(agb->agbw_p[iproc]);
	}
	agb->agbw_p[iproc] = NULL;
      }
      free(agb->agbw_p); agb->agbw_p = NULL;
    }
    if(agb->cavw_p){
      for(iproc=0;iproc<agb->maxprocs;iproc++){
	if(agb->cavw_p[iproc]){
	  agbnp_delete_cavworkdata(agb->cavw_p[iproc]);
	}
	agb->cavw_p[iproc] = NULL;
      }
      free(agb->cavw_p); agb->cavw_p = NULL;
    }
  }
#endif

  if(agb->occupancy) { free(agb->occupancy); agb->occupancy = NULL; }
  if(agb->alt_site) { free(agb->alt_site); agb->alt_site = NULL; }
  if(agb->alt_id) { free(agb->alt_id); agb->alt_id = NULL; }

  agbnp_reset(agb);
  return AGBNP_OK;
}

/* return agb+vdw energies and derivatives */
int agbnp_agb_energy(int tag, float_i *x, float_i *y, float_i *z,
		     float_i *sp, float_i *br, float_i *egb, float_i (*dgbdr)[3],
		     float_i *evdw, float_i (*dvwdr)[3], float_i *ecorr){
  AGBNPdata *agb;

  if(!agbnp_initialized){
    fprintf(stderr,"agbnp_agb_energy(): agbnp library is not initialized.\n");
    return AGBNP_ERR;
  }
  if(!agbnp_tag_ok(tag)){
    fprintf(stderr,"agbnp_agb_energy(): invalid tag %d.\n",tag);
    return AGBNP_ERR;
  }

  /* pointer to agb data structure */
  agb = &(agbdata_list[tag]);

  if(agb_vdw_energy(agb, x, y, z, sp, br, egb, dgbdr,
		    evdw, dvwdr, ecorr) != AGBNP_OK){
    fprintf(stderr,"agbnp_agb_energy(): error in agb_vdw_energy().\n");
    return AGBNP_ERR;
  }

  return AGBNP_OK;
}

int plop_agbnp_agb_energy(int tag, int natom_asu, float_i *x, float_i *y, float_i *z,
                          float_i *sp, float_i *br, float_i *evdw, float_i *ecorr){
  AGBNPdata *agb;

  if(!agbnp_initialized){
    fprintf(stderr,"plop_agbnp_agb_energy(): agbnp library is not initialized.\n");
    return AGBNP_ERR;
  }
  if(!agbnp_tag_ok(tag)){
    fprintf(stderr,"plop_agbnp_agb_energy(): invalid tag %d.\n",tag);
    return AGBNP_ERR;
  }

  /* pointer to agb data structure */
  agb = &(agbdata_list[tag]);

  if(plop_agb_vdw_energy(agb, natom_asu, x, y, z, sp, br, evdw, ecorr) != AGBNP_OK){
    fprintf(stderr,"plop_agbnp_agb_energy(): error in agb_vdw_energy().\n");
    return AGBNP_ERR;
  }

  return AGBNP_OK;
}


/* returns cavity energy and derivatives */
int agbnp_cavity_energy(int tag, float_i *x, float_i *y, float_i *z,
			float_i *mol_volume, float_i *surf_area, 
			float_i *ecav, float_i *ecorr,
			float_i (*decav)[3]){
  AGBNPdata *agb;

  if(!agbnp_initialized){
    fprintf(stderr,"agbnp_cavity_energy(): agbnp library is not initialized.\n");
    return AGBNP_ERR;
  }
  if(!agbnp_tag_ok(tag)){
    fprintf(stderr,"agbnp_cavity_energy(): invalid tag %d.\n",tag);
    return AGBNP_ERR;
  }

  /* pointer to agb data structure */
  agb = &(agbdata_list[tag]);

  if(agb_cavity_energy(agb, x, y, z, mol_volume, surf_area, 
		       ecav, ecorr, decav) != AGBNP_OK){
    fprintf(stderr,"agbnp_cavity_energy(): error in agb_cavity_energy().\n");
    return AGBNP_ERR;
  }

  return AGBNP_OK;  
}

int plop_agbnp_cavity_energy(int tag, int natom_asu, float_i *x, float_i *y, float_i *z,
                             float_i *mol_volume, float_i *surf_area, 
                             float_i *ecav, float_i *ecorr) {
  AGBNPdata *agb;

  if(!agbnp_initialized){
    fprintf(stderr,"plop_agbnp_cavity_energy(): agbnp library is not initialized.\n");
    return AGBNP_ERR;
  }
  if(!agbnp_tag_ok(tag)){
    fprintf(stderr,"plop_agbnp_cavity_energy(): invalid tag %d.\n",tag);
    return AGBNP_ERR;
  }

  /* pointer to agb data structure */
  agb = &(agbdata_list[tag]);

  if(plop_agb_cavity_energy(agb, natom_asu, x, y, z, mol_volume, surf_area, ecav, ecorr) != AGBNP_OK){
    fprintf(stderr,"plop_agbnp_cavity_energy(): error in agb_cavity_energy().\n");
    return AGBNP_ERR;
  }

  return AGBNP_OK;  
}


/* (re)-initialize constant terms from frozen atoms */
int agbnp_init_frozen(int tag, float_i *x, float_i *y, float_i *z, int *isfrozen){

  AGBNPdata *agb;

  if(!agbnp_initialized){
    fprintf(stderr,"agbnp_init_frozen(): agbnp library is not initialized.\n");
    return AGBNP_ERR;
  }
  if(!agbnp_tag_ok(tag)){
    fprintf(stderr,"agbnp_init_frozen(): invalid tag %d.\n",tag);
    return AGBNP_ERR;
  }

  /* pointer to agb data structure */
  agb = &(agbdata_list[tag]);

  /* if isfrozen is undefined just turn off frozen atoms optimizations */
  if(!isfrozen){
    agb->do_frozen = 0;
    return AGBNP_OK;
  }

  /* reset frozen atom list */
  agb->do_frozen = 1;
  agb->isfrozen = (int *)calloc(agb->natoms,sizeof(int));
  if(!agb->isfrozen){
    fprintf(stderr,"agbnp_init_frozen(): unable to allocate isfrozen array (%d ints)\n",agb->natoms);
    return AGBNP_ERR;
  }
  memcpy(agb->isfrozen, isfrozen, agb->natoms*sizeof(int));

  /* initialize cavity free energy constant terms */
  if(agbnp_init_frozen_cavity(agb, x, y, z) != AGBNP_OK){
    fprintf(stderr,"agbnp_init_frozen(): error in agb_init_frozen_cavity().\n");
    return AGBNP_ERR;
  }

  /* initialize AGB/wdW free energy constant terms */
  if(agbnp_init_frozen_agb(agb, x, y, z) != AGBNP_OK){
    fprintf(stderr,"agbnp_init_frozen(): error in agb_init_frozen_agb().\n");
    return AGBNP_ERR;
  }

  return AGBNP_OK;
}

/* set alternate conformation information */
int agbnp_set_altconfs(int tag, float_i *occupancy, int *alt_site, int *alt_id){
  AGBNPdata *agb;
  int natoms, iat;

  if(!agbnp_initialized){
    fprintf(stderr,"agbnp_set_altconfs(): agbnp library is not initialized.\n");
    return AGBNP_ERR;
  }
  if(!agbnp_tag_ok(tag)){
    fprintf(stderr,"agbnp_set_altconfs(): invalid tag %d.\n",tag);
    return AGBNP_ERR;
  }
  /* pointer to agb data structure */
  agb = &(agbdata_list[tag]);

  /* allow alternate conformations only if PBC is on, otherwise will
     need optimized routines able to do alternate conformations
     without coordinates from neighbor list */
  if(!agb->dopbc){
    fprintf(stderr,"agbnp_set_altconfs(): sorry AGBNP supports alternate conformations only with periodic boundary conditions. If you really need alternate conformations but not PBCs try setting PBC's with a very large simulation box.\n");
    return AGBNP_ERR;  
  }
  
  natoms = agb->natoms;
  agb->occupancy = (float_a *)realloc(agb->occupancy, natoms*sizeof(float_a));
  agb->alt_site = (int *)realloc(agb->alt_site, natoms*sizeof(int));
  agb->alt_id = (int *)realloc(agb->alt_id, natoms*sizeof(int));
  if(!(agb->occupancy && agb->alt_site && agb->alt_id)){
    fprintf(stderr, "agbnp_set_altconfs(): unable to allocate alternate conformations buffers (%d doubles and %d ints).\n", natoms, 2*natoms);
    return AGBNP_ERR;
  }
  agbnp_cpyd(agb->occupancy, occupancy, natoms);
  memcpy(agb->alt_site, alt_site, natoms*sizeof(int));
  memcpy(agb->alt_id, alt_id, natoms*sizeof(int));

  /* re-initializes work arrays vols ... to take into account occupancies */
  agbnp_init_agbworkdata(agb,agb->agbw);
  agbnp_init_cavworkdata(agb,agb->cavw);
#ifdef _OPENMP
 {
   int iproc;
   /* init work arrays for each thread */
   agb->nprocs = omp_get_num_threads();
   for(iproc=0;iproc<agb->maxprocs;iproc++){
     agbnp_init_agbworkdata(agb,agb->agbw_p[iproc]);
     agbnp_init_cavworkdata(agb,agb->cavw_p[iproc]);
   }
 }
#endif

  agb->doalt = 1;
  return AGBNP_OK;
}

/* check if it is a valid tag */
 int agbnp_tag_ok(int tag){
  /* check ranges */
  if(tag < 0) return FALSE;
  if(tag >= agbdata_allocated) return FALSE;
  
  /* check if list is in use */
  if(!agbdata_list[tag].in_use) return FALSE;

  return TRUE;
}


/* reset an agbnp structure */
 int agbnp_reset(AGBNPdata *data){
  data->in_use = FALSE;
  data->natoms = 0;
  data->x = data->y = data->z = NULL;
  data->r = NULL;
  data->rcav = NULL;
  data->charge = NULL;
  data->igamma = data->sgamma = NULL;
  data->ialpha = data->salpha = NULL;
  data->idelta = data->sdelta = NULL;
  data->ab = NULL;
#ifdef AGBNP_HB
  data->hbtype = NULL;
  data->hbcorr = NULL;
#endif
  data->nheavyat = 0;
  data->iheavyat = NULL;
  data->nhydrogen = 0;
  data->ihydrogen = NULL;
  data->ndummy = 0;
  data->idummy = NULL;
  data->isfrozen = NULL;
  data->dielectric_in = data->dielectric_out = -1.0;
  data->neigh_list = data->excl_neigh_list = NULL;
  data->dopbc = 0;
  data->nsym = 0;
  data->docryst = 0;
  data->ssize = 0;
  data->xs = data->ys = data->zs = NULL;
  data->rot = NULL;
  data->conntbl = NULL;
#ifdef AGBNP_HB
  data->do_w = 0;
  data->nwsat = 0;
  data->wsat_size = 0;
  data->wsat = NULL;
  data->wsat_list = NULL;
#endif
  data->do_ms = 0;
  data->nmsat = 0;
  data->msat_size = 0;
  data->msat = NULL;
  data->msat_list = NULL;
  data->agbw = NULL;
  data->cavw = NULL;
  data->nprocs = data->maxprocs = 0;
  data->agbw_p = NULL;
  data->cavw_p = NULL;
  data->doalt = 0;
  data->occupancy = NULL;
  data->alt_site = NULL;
  data->alt_id = NULL;
  return AGBNP_OK;
}

 int agbnp_reset_agbworkdata(AGBworkdata *agbw){
  agbw->natoms = 0;
  agbw->vols = agbw->volumep = NULL;
  agbw->dera = agbw->deru = agbw->derv = NULL;
  agbw->q2ab = agbw->abrw = NULL;
  agbw->br1_swf_der = NULL;
  agbw->isheavy = NULL;
  agbw->nbiat = NULL;
  agbw->nbdata = NULL;
  agbw->xj = agbw->yj = agbw->zj = NULL;
  agbw->br = agbw->br1 = agbw->brw = NULL;
  agbw->alpha = agbw->delta = NULL;
  agbw->galpha = agbw->gprefac = NULL;
  agbw->sp = NULL;
  agbw->isvfrozen = NULL;
  agbw->isvvfrozen = NULL;
  agbw->isbfrozen = NULL;
  agbw->volumep_constant = NULL;
  agbw->br1_const = NULL;
  agbw->gbpair_const = NULL;
  agbw->dera_const = agbw->deru_c = agbw->derv_c = NULL;
  agbw->nq4cache = 0;
  agbw->q4cache = NULL;
  agbw->nv2cache = 0;
  agbw->v2cache = NULL;
  agbw->near_nl = NULL;
  agbw->far_nl = NULL;
  agbw->dgbdr_h = NULL;
  agbw->dvwdr_h = NULL;
  agbw->dgbdr_c = NULL;
  agbw->dvwdr_c = NULL;
  agbw->nl_r2v = NULL;
  agbw->nl_indx = NULL;
  agbw->nlist = NULL;
  agbw->js = NULL;
  agbw->pbcs = NULL;
  agbw->datas = NULL;
  return AGBNP_OK;
}

 int agbnp_reset_cavworkdata(CAVworkdata *cavw){
  cavw->natoms = 0;
  cavw->vols = cavw->volumep = NULL;
  cavw->isheavy = NULL;
  cavw->nbiat = NULL;
  cavw->nbdata = NULL;
  cavw->xj = cavw->yj = cavw->zj = NULL;
  cavw->gamma = cavw->gammap = NULL;
  cavw->galpha = cavw->gprefac = NULL;
  cavw->isvfrozen = NULL;
  cavw->isvvfrozen = NULL;
  cavw->volumep_constant = NULL;
  cavw->surfarea_constant = NULL;
  cavw->nv2cache = 0;
  cavw->v2cache = NULL;
  cavw->near_nl = NULL;
  cavw->decav_h = NULL;
  cavw->decav_const = NULL;
  cavw->nl_r2v = NULL;
  cavw->nl_indx = NULL;
  cavw->nlist = NULL;
  cavw->js = NULL;
  cavw->pbcs = NULL;
  cavw->datas = NULL;
  return AGBNP_OK;
}

 int agbnp_allocate_agbworkdata(int natoms, AGBNPdata *agb, 
				      AGBworkdata *agbw){
  int nsym = agb->nsym;
  agbw->natoms = natoms;
  agbw->vols = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->volumep = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->dera = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->deru = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->derv = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->q2ab = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->abrw = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->br1_swf_der = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->isheavy = (int *)calloc(natoms,sizeof(int));
  agbw->nbiat = (int *)calloc(nsym*natoms,sizeof(int));
  agbw->nbdata = (void **)calloc(nsym*natoms,sizeof(void *));
  agbw->xj = (float_a *)calloc(nsym*natoms,sizeof(float_a));
  agbw->yj = (float_a *)calloc(nsym*natoms,sizeof(float_a));
  agbw->zj = (float_a *)calloc(nsym*natoms,sizeof(float_a));
  agbw->br = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->br1 = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->brw = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->alpha = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->delta = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->galpha = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->gprefac = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->sp = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->isvfrozen = (int *)calloc(natoms,sizeof(int));
  agbw->isvvfrozen = (int *)calloc(natoms,sizeof(int));
  agbw->isbfrozen = (int *)calloc(natoms,sizeof(int));
  agbw->volumep_constant = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->br1_const = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->gbpair_const = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->dera_const = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->deru_c = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->derv_c = (float_a *)calloc(natoms,sizeof(float_a));
  agbw->dgbdr_h = (float_a (*)[3])calloc(natoms,sizeof(float_a [3]));
  agbw->dvwdr_h = (float_a (*)[3])calloc(natoms,sizeof(float_a [3]));
  agbw->dgbdr_c = (float_a (*)[3])calloc(natoms,sizeof(float_a [3]));
  agbw->dvwdr_c = (float_a (*)[3])calloc(natoms,sizeof(float_a [3]));
  agbw->nl_r2v = (float_a *)calloc(nsym*natoms, sizeof(float_a));
  agbw->nl_indx = (int *)calloc(nsym*natoms, sizeof(int));
  agbw->nlist = (int *)calloc(nsym*natoms,sizeof(int));
  agbw->js = (int *)calloc(nsym*natoms,sizeof(int));
  agbw->pbcs = (NeighVector *)calloc(nsym*natoms,sizeof(NeighVector));
  agbw->datas = (void **)calloc(nsym*natoms,sizeof(void *));
  if(!(agbw->vols && agbw->volumep && agbw->dera && agbw->deru && 
       agbw->derv && agbw->q2ab && agbw->abrw && agbw->br1_swf_der && 
       agbw->isheavy && agbw->nbiat && agbw->nbdata && 
       agbw->xj && agbw->yj && agbw->zj &&
       agbw->br && agbw->br1 && agbw->brw  && 
       agbw->alpha && agbw->delta && agbw->galpha && agbw->gprefac && 
       agbw->sp && agbw->isvfrozen && agbw->isvvfrozen && agbw->isbfrozen && 
       agbw->volumep_constant && 
       agbw->br1_const && agbw->gbpair_const &&  agbw->dera_const &&
       agbw->deru_c && agbw->derv_c &&
       agbw->dgbdr_h && agbw->dvwdr_h && agbw->dgbdr_c && agbw->dvwdr_c &&
       agbw->nl_r2v && agbw->nl_indx &&
       agbw->nlist && agbw->js && agbw->pbcs && agbw->datas)){
    fprintf(stderr, "agbnp_allocate_agbworkdata(): unable to allocate memory for AGB work data structure.\n");
    return AGBNP_ERR;
  }

  agbw->near_nl = (NeighList *)malloc(1*sizeof(NeighList));
  if(!agbw->near_nl){
    fprintf(stderr, "agbnp_allocate_agbworkdata(): unable to allocate memory for near_nl (%d NeighList).\n", 1);
    return AGBNP_ERR;
  }
  nblist_reset_neighbor_list(agbw->near_nl);
  if(agb->dopbc){
    /* to store pointers to PBC translation vectors */
    agbw->near_nl->data = 1;
  }
  if(nblist_reallocate_neighbor_list(agbw->near_nl, natoms,
				     natoms*AGBNP_NEARNEIGHBORS) != NBLIST_OK){
      fprintf(stderr,"agbnp_allocate_agbworkdata(): unable to allocate near_nl neighbor list (natoms=%d, size=%d)\n", natoms, natoms*AGBNP_NEARNEIGHBORS);
      return AGBNP_ERR;
  }

  agbw->far_nl = (NeighList *)malloc(1*sizeof(NeighList));
  if(!agbw->far_nl){
    fprintf(stderr, "agbnp_allocate_agbworkdata(): unable to allocate memory for far_nl (%d NeighList).\n", 1);
    return AGBNP_ERR;
  }
  nblist_reset_neighbor_list(agbw->far_nl);
  if(agb->dopbc){
    /* to store pointers to PBC translation vectors */
    agbw->far_nl->data = 1;
  }
  if(nblist_reallocate_neighbor_list(agbw->far_nl, natoms,
				     natoms*AGBNP_FARNEIGHBORS) != NBLIST_OK){
      fprintf(stderr,"agbnp_allocate_agbworkdata(): unable to allocate far_nl neighbor list (natoms=%d, size=%d)\n", natoms, natoms*AGBNP_FARNEIGHBORS);
      return AGBNP_ERR;
  }
  
  return AGBNP_OK;
}

 int agbnp_delete_agbworkdata(AGBworkdata *agbw){
  if(agbw->vols){ free(agbw->vols); agbw->vols = NULL;}
  if(agbw->volumep){free(agbw->volumep); agbw->volumep = NULL;}
  if(agbw->dera){ free(agbw->dera); agbw->dera = NULL;}
  if(agbw->deru){ free(agbw->deru); agbw->deru = NULL;}
  if(agbw->derv){ free(agbw->derv); agbw->derv = NULL;}
  if(agbw->q2ab){ free(agbw->q2ab); agbw->q2ab = NULL;}
  if(agbw->abrw){ free(agbw->abrw); agbw->abrw = NULL;}
  if(agbw->br1_swf_der){free(agbw->br1_swf_der); agbw->br1_swf_der = NULL;}
  if(agbw->isheavy){free(agbw->isheavy); agbw->isheavy = NULL;}
  if(agbw->nbiat){free(agbw->nbiat); agbw->nbiat = NULL;}
  if(agbw->nbdata){free(agbw->nbdata); agbw->nbdata = NULL;}
  if(agbw->xj){free(agbw->xj); agbw->xj = NULL;}
  if(agbw->yj){free(agbw->yj); agbw->yj = NULL;}
  if(agbw->zj){free(agbw->zj); agbw->zj = NULL;}
  if(agbw->br){free(agbw->br); agbw->br = NULL;}
  if(agbw->br1){free(agbw->br1); agbw->br1 = NULL;}
  if(agbw->brw){free(agbw->brw); agbw->brw = NULL;}
  if(agbw->alpha){free(agbw->alpha); agbw->alpha = NULL;}
  if(agbw->delta){free(agbw->delta); agbw->delta = NULL;}
  if(agbw->galpha){free(agbw->galpha); agbw->galpha = NULL;}
  if(agbw->gprefac){free(agbw->gprefac); agbw->gprefac = NULL;}
  if(agbw->sp){free(agbw->sp); agbw->sp = NULL;}
  if(agbw->isvfrozen){free(agbw->isvfrozen); agbw->isvfrozen = NULL;}
  if(agbw->isvvfrozen){free(agbw->isvvfrozen); agbw->isvvfrozen = NULL;}
  if(agbw->isbfrozen){free(agbw->isbfrozen); agbw->isbfrozen = NULL;}
  if(agbw->volumep_constant){
    free(agbw->volumep_constant); agbw->volumep_constant = NULL;}
  if(agbw->br1_const){free(agbw->br1_const); agbw->br1_const = NULL;}
  if(agbw->gbpair_const){free(agbw->gbpair_const); agbw->gbpair_const = NULL;}
  if(agbw->dera_const){free(agbw->dera_const); agbw->dera_const = NULL;}
  if(agbw->deru_c){free(agbw->deru_c); agbw->deru_c = NULL;}
  if(agbw->derv_c){free(agbw->derv_c); agbw->derv_c = NULL;}
  if(agbw->q4cache){free(agbw->q4cache); agbw->q4cache = NULL;}
  if(agbw->v2cache){free(agbw->v2cache); agbw->v2cache = NULL;}
  if(agbw->dgbdr_h){free(agbw->dgbdr_h); agbw->dgbdr_h = NULL;}
  if(agbw->dvwdr_h){free(agbw->dvwdr_h); agbw->dvwdr_h = NULL;}
  if(agbw->dgbdr_c){free(agbw->dgbdr_c); agbw->dgbdr_c = NULL;}
  if(agbw->dvwdr_c){free(agbw->dvwdr_c); agbw->dvwdr_c = NULL;}

  if(agbw->nl_r2v)  {free(agbw->nl_r2v); agbw->nl_r2v = NULL;}
  if(agbw->nl_indx) {free(agbw->nl_indx); agbw->nl_indx = NULL;}
  if(agbw->nlist)  {free(agbw->nlist); agbw->nlist  = NULL;}
  if(agbw->js)  {free(agbw->js); agbw->js  = NULL;}
  if(agbw->pbcs)  {free(agbw->pbcs); agbw->pbcs  = NULL;}
  if(agbw->datas)  {free(agbw->datas); agbw->datas  = NULL;}  

  if(agbw->near_nl){
    nblist_delete_neighbor_list(agbw->near_nl);
    agbw->near_nl = NULL;
  }
  if(agbw->far_nl){
    nblist_delete_neighbor_list(agbw->far_nl);
    agbw->far_nl = NULL;
  }
  return AGBNP_OK;
}

 int agbnp_allocate_cavworkdata(int natoms, AGBNPdata *agb,
				      CAVworkdata *cavw){
  int nsym = agb->nsym;
  cavw->natoms = natoms;
  cavw->vols = (float_a *)calloc(natoms,sizeof(float_a));
  cavw->volumep = (float_a *)calloc(natoms,sizeof(float_a));
  cavw->surf_area = (float_a *)calloc(natoms,sizeof(float_a));
  cavw->isheavy = (int *)calloc(natoms,sizeof(int));
  cavw->nbiat = (int *)calloc(nsym*natoms,sizeof(int));
  cavw->nbdata = (void **)calloc(nsym*natoms,sizeof(void *));
  cavw->xj = (float_a *)calloc(nsym*natoms,sizeof(float_a));
  cavw->yj = (float_a *)calloc(nsym*natoms,sizeof(float_a));
  cavw->zj = (float_a *)calloc(nsym*natoms,sizeof(float_a));
  cavw->gamma = (float_a *)calloc(natoms,sizeof(float_a));
  cavw->gammap = (float_a *)calloc(natoms,sizeof(float_a));
  cavw->galpha = (float_a *)calloc(natoms,sizeof(float_a));
  cavw->gprefac = (float_a *)calloc(natoms,sizeof(float_a));
  cavw->isvfrozen = (int *)calloc(natoms,sizeof(int));
  cavw->isvvfrozen = (int *)calloc(natoms,sizeof(int));
  cavw->volumep_constant = (float_a *)calloc(natoms,sizeof(float_a));
  cavw->surfarea_constant = (float_a *)calloc(natoms,sizeof(float_a));
  cavw->decav_h = (float_a (*)[3])calloc(natoms,sizeof(float_a [3]));
  cavw->decav_const = (float_a (*)[3])calloc(natoms,sizeof(float_a [3]));
  cavw->nl_r2v = (float_a *)calloc(nsym*natoms, sizeof(float_a));
  cavw->nl_indx = (int *)calloc(nsym*natoms, sizeof(int));
  cavw->nlist = (int *)calloc(nsym*natoms,sizeof(int));
  cavw->js = (int *)calloc(nsym*natoms,sizeof(int));
  cavw->pbcs = (NeighVector *)calloc(nsym*natoms,sizeof(NeighVector));
  cavw->datas = (void **)calloc(nsym*natoms,sizeof(void *));
  if(!(cavw->vols && cavw->volumep &&  
       cavw->isheavy && cavw->nbiat && cavw->nbdata &&
       cavw->xj && cavw->yj && cavw->zj && 
       cavw->gamma && cavw->gammap && cavw->galpha && cavw->gprefac && 
       cavw->isvfrozen && cavw->isvvfrozen && 
       cavw->volumep_constant && cavw->surfarea_constant &&
       cavw->decav_h && cavw->decav_const &&
       cavw->nl_r2v && cavw->nl_indx &&
       cavw->nlist && cavw->js && cavw->pbcs && cavw->datas)){
    fprintf(stderr, "agbnp_allocate_cavworkdata(): unable to allocate memory for Cavity work data structure.\n");
    return AGBNP_ERR;
  }

  cavw->near_nl = (NeighList *)malloc(1*sizeof(NeighList));
  if(!cavw->near_nl){
    fprintf(stderr, "agbnp_allocate_cavworkdata(): unable to allocate memory for near_nl (%d NeighList).\n", 1);
    return AGBNP_ERR;
  }
  nblist_reset_neighbor_list(cavw->near_nl);
  if(agb->dopbc){
    /* to store pointers to PBC translation vectors */
    cavw->near_nl->data = 1;
  }
  if(nblist_reallocate_neighbor_list(cavw->near_nl, natoms,
				     natoms*AGBNP_NEARNEIGHBORS) != NBLIST_OK){
      fprintf(stderr,"agbnp_allocate_cavworkdata(): unable to allocate near_nl neighbor list (natoms=%d, size=%d)\n", natoms, natoms*AGBNP_NEARNEIGHBORS);
      return AGBNP_ERR;
  }
  
  return AGBNP_OK;
}

 int agbnp_delete_cavworkdata(CAVworkdata *cavw){
  if(cavw->vols){ free(cavw->vols); cavw->vols = NULL;}
  if(cavw->volumep){free(cavw->volumep); cavw->volumep = NULL;}
  if(cavw->surf_area){free(cavw->surf_area); cavw->surf_area = NULL;}
  if(cavw->isheavy){free(cavw->isheavy); cavw->isheavy = NULL;}
  if(cavw->nbiat){free(cavw->nbiat); cavw->nbiat = NULL;}
  if(cavw->nbdata){free(cavw->nbdata); cavw->nbdata = NULL;}
  if(cavw->xj){free(cavw->xj); cavw->xj = NULL;}
  if(cavw->yj){free(cavw->yj); cavw->yj = NULL;}
  if(cavw->zj){free(cavw->zj); cavw->zj = NULL;}
  if(cavw->gamma){free(cavw->gamma); cavw->gamma = NULL;}
  if(cavw->gammap){free(cavw->gammap); cavw->gammap = NULL;}
  if(cavw->galpha){free(cavw->galpha); cavw->galpha = NULL;}
  if(cavw->gprefac){free(cavw->gprefac); cavw->gprefac = NULL;}
  if(cavw->isvfrozen){free(cavw->isvfrozen); cavw->isvfrozen = NULL;}
  if(cavw->isvvfrozen){free(cavw->isvvfrozen); cavw->isvvfrozen = NULL;}
  if(cavw->volumep_constant){
    free(cavw->volumep_constant); cavw->volumep_constant = NULL;}
  if(cavw->surfarea_constant){
    free(cavw->surfarea_constant); cavw->surfarea_constant = NULL;}
  if(cavw->decav_h){free(cavw->decav_h); cavw->decav_h = NULL;}
  if(cavw->decav_const){free(cavw->decav_const); cavw->decav_const = NULL;}

  if(cavw->nl_r2v)  {free(cavw->nl_r2v); cavw->nl_r2v = NULL;}
  if(cavw->nl_indx) {free(cavw->nl_indx); cavw->nl_indx = NULL;}
  if(cavw->nlist)  {free(cavw->nlist); cavw->nlist  = NULL;}
  if(cavw->js)  {free(cavw->js); cavw->js  = NULL;}
  if(cavw->pbcs)  {free(cavw->pbcs); cavw->pbcs  = NULL;}
  if(cavw->datas)  {free(cavw->datas); cavw->datas  = NULL;} 

  if(cavw->near_nl){
    nblist_delete_neighbor_list(cavw->near_nl);
    cavw->near_nl = NULL;
  }
  return AGBNP_OK;
}

 int agbnp_init_agbworkdata(AGBNPdata *agb, AGBworkdata *agbw){
  int i, iat;
  float_a c4;

  /* set isheavy array */
  memset(agbw->isheavy,0,agbw->natoms*sizeof(int));
  for(i=0;i<agb->nheavyat;i++){
    agbw->isheavy[agb->iheavyat[i]] = 1;
  }

  /* set np parameters */
  for(iat=0;iat<agbw->natoms;iat++){
    agbw->alpha[iat] = agb->ialpha[iat] + agb->salpha[iat];
    agbw->delta[iat] = agb->idelta[iat] + agb->sdelta[iat];
  }

  /* initialize atomic volumes */
  c4 = 4.0*pi/3.0;
  for(iat=0;iat<agbw->natoms;iat++){
    agbw->vols[iat] = c4*agb->occupancy[iat]*pow(agb->r[iat],3);
  }
  
  /* Gaussian parameters initialization */
  for(i=0;i<agb->nheavyat;i++){
    iat = agb->iheavyat[i];
    agbw->galpha[iat] = KFC/(agb->r[iat]*agb->r[iat]);
    agbw->gprefac[iat] = PFC;
  }

  return AGBNP_OK;
}

 int agbnp_init_cavworkdata(AGBNPdata *agb, CAVworkdata *cavw){
  int i, iat;
  float_a c4;

  /* set isheavy array */
  memset(cavw->isheavy,0,cavw->natoms*sizeof(int));
  for(i=0;i<agb->nheavyat;i++){
    cavw->isheavy[agb->iheavyat[i]] = 1;
  }

  /* set np parameters */
  for(iat=0;iat<cavw->natoms;iat++){
    cavw->gamma[iat] = agb->igamma[iat] + agb->sgamma[iat];
  }

  /* initialize atomic volumes */
  c4 = 4.0*pi/3.0;
  for(iat=0;iat<cavw->natoms;iat++){
    cavw->vols[iat] = c4*agb->occupancy[iat]*pow(agb->rcav[iat],3);
  }

  /* Gaussian parameters initialization */
  for(i=0;i<agb->nheavyat;i++){
    iat = agb->iheavyat[i];
    cavw->galpha[iat] = KFC/(agb->rcav[iat]*agb->rcav[iat]);
    cavw->gprefac[iat] = PFC;
  }

  return AGBNP_OK;
}

/*
 float_a i4_far(float_a rij, float_a Ri, float_a Rj, 
		     float_a *dr, float_a *dRj){
  float_a u1, u2, u3,u4;
  float_a rij2 = rij*rij;
  float_a q;
  static const float_a twopi = 2.0*pi;

  u1 = rij+Rj;
  u2 = rij-Rj;
  u3 = u1*u2;
  u4 = 0.5*log(u1/u2);
  q = twopi*(Rj/u3 - u4/rij);
  *dr =  twopi*( (Rj/(rij*u3))*(1. - 2.*rij2/u3 ) + u4/rij2 );
  *dRj = twopi*( 2.0*Rj*Rj/(u3*u3) );

  return q;
}
*/


/* gives the alternate conformation weight for each interaction */
 float_a altweight_pair(AGBNPdata *agb, int iat, int jat){
  if(!agb->doalt){
    return 1.0;
  }
  if(agb->alt_site[iat]==agb->alt_site[jat]){
    if(agb->alt_id[iat]==agb->alt_id[jat]){
      return agb->occupancy[iat];
    }else{
      return 0.0;
    }
  }
  return agb->occupancy[iat]*agb->occupancy[jat];
}
 float_a altweight_incremental(AGBNPdata *agb, int nat, 
				    int *iat, float_a iat_occupancy, int jat){
  int i;
  int *alt_site, *alt_id;
  int samesite = 0;
  if(!agb->doalt){
    return 1.0;
  }
  alt_site=agb->alt_site;
  alt_id=agb->alt_id;
  for(i=0;i<nat;i++){
    if(alt_site[iat[i]]==alt_site[jat]){
      if(alt_id[iat[i]]!=alt_id[jat]){
	return 0.0;
      }else{
	samesite=1;
      }
    }
  }
  if(samesite){
    return iat_occupancy;
  }
  return iat_occupancy*agb->occupancy[jat];
}

#ifdef VOLUME_OVERLAP
 float_a volume_overlap(float_a (*x)[3], float_a *r,
			     float_a (*dr)[3], float_a *dR, 
			     float_a (*d2rR)[AGBNP_MAX_OVERLAP_LEVEL][3]){
  float_a dji[3], rji2, rji,Rij;
  float_a Ri, Rj;
  float_a cossji, cossij;
  float_a hi,hj,hi2,hj2, Rhi, Rhj, vol;
  float_a dwri,dwrj;
  float_a dwriRi, dwriRj, dwrjRi, dwrjRj;
  int i;
  float_a pi3 =  pi/3.0;
  
 
  dji[0] = x[1][0]-x[0][0];
  dji[1] = x[1][1]-x[0][1];
  dji[2] = x[1][2]-x[0][2];
  
  rji2 = dji[0]*dji[0]+dji[1]*dji[1]+dji[2]*dji[2];

  dr[0][0] = dr[0][1] = dr[0][2] = 0.0;
  dr[1][0] = dr[1][1] = dr[1][2] = 0.0;
  dR[0] = dR[1] = 0.0;
  memset(d2rR,0,12*sizeof(float_a));

  /* distance test for separated spheres */
  Ri = r[0];
  Rj = r[1];
  Rij = Ri+Rj;
  if( rji2 > Rij*Rij ) return 0.0;
  /* distance test for totally included spheres */
  Rij = Ri-Rj;
  if( rji2 <= Rij*Rij ){
    if(Rij < 0.0){           /* i is smaller */
      vol = 4.*pi3*Ri*Ri*Ri;
      dR[0] = 4.*pi*Ri*Ri;
    }else{                   /* j is smaller */
      vol = 4.*pi3*Rj*Rj*Rj;
      dR[1] = 4.*pi*Rj*Rj;
    }
    return vol;
  }

  rji = sqrt(rji2);

  cossji = (Ri*Ri - Rj*Rj + rji2)/(2.*rji*Ri);
  cossij = (Rj*Rj - Ri*Ri + rji2)/(2.*rji*Rj);

  hi = Ri*(1.-cossji); hi2 = hi*hi; Rhi = (2.*Ri-hi);
  hj = Rj*(1.-cossij); hj2 = hj*hj; Rhj = (2.*Rj-hj);

  vol = pi3*( hi2*(3.*Ri-hi) + hj2*(3.*Rj-hj) );

  /* positional derivatives */
  dwri = pi*hi*Rhi/rji;
  dwrj = -dwri;

  dr[0][0] = dji[0]*dwri;
  dr[0][1] = dji[1]*dwri;
  dr[0][2] = dji[2]*dwri;

  dr[1][0] = dji[0]*dwrj;
  dr[1][1] = dji[1]*dwrj;
  dr[1][2] = dji[2]*dwrj;
  
  /* radius derivatives */
  dR[0] = 2. * pi * Ri * hi;
  dR[1] = 2. * pi * Rj * hj;

  /* pos-radius mixed derivatives */
  dwriRi =  2. * pi * Ri * ( Rj - hj)/rji2;
  dwriRj =  2. * pi * Rj * ( Ri - hi)/rji2;
  dwrjRi = -dwriRi;
  dwrjRj = -dwriRj;

  for(i=0;i<3;i++){
    d2rR[0][0][i] = dwriRi*dji[i];
    d2rR[0][1][i] = dwriRj*dji[i];
    d2rR[1][0][i] = dwrjRi*dji[i];
    d2rR[1][1][i] = dwrjRj*dji[i];
  }

  return vol;
}
#endif

 float_a swf_area(float_a x, float_a *fp){
  static const float_a a2 = 5.*5.;
  float_a t, f, x2;
 
  if(x<0.0){
    *fp = 0.0;
    return (float_a)0.0;
  }
  x2 = x*x;
  t = x/(a2 + x2);
  f = x*t;
  *fp  = (2.*t)*(1. - f);
  return f;
}

 float_a agbnp_pol_switchfunc(float_a x, float_a xa, float_a xb,
				   float_a *fp, float_a *fpp){
  float_a u,d,u2,u3,f;
  if(x > xb) {
    if(fp) *fp = 0.0;
    if(fpp) *fpp = 0.0;
    return 1.0;
  }
  if(x < xa) {
    if(fp) *fp = 0.0;
    if(fpp) *fpp = 0.0;
    return 0.0;
  }
  d = 1./(xb - xa);
  u = (x - xa)*d;
  u2 = u*u;
  u3 = u*u2;
  f = u3*(10.-15.*u+6*u2);
  if(fp){   /* first derivative */
    *fp = d*30.*u2*(1. - 2.*u + u2);
  }
  if(fpp){ /* second derivative */
    *fpp = d*d*60.*u*(1. - 3.*u + 2.*u2);
  }
  return f;
}


 float_a swf_vol3(float_a x, float_a *fp, float_a *fpp, 
		       float_a a, float_a b){
  float_a f, s, sp, spp;

  if(x>b){
    *fp = 1.0;
    *fpp = 0.0;
    return x;
  }
  if(x<a){
    *fp = 0.0;
    *fpp = 0.0;
    return (float_a)0.0;
  }
  s = agbnp_pol_switchfunc(x, a, b, &sp, &spp);
  f = s*x;
  *fp = s + x*sp;
  *fpp = 2.*sp + x*spp;

  return f;
}

 float_a ogauss_2body_smpl(float_a d2, float_a p1, float_a p2, 
			   float_a c1, float_a c2){
  float_a deltai = 1./(c1+c2);
  float_a p = p1*p2;
  float_a kappa, gvol;

#ifdef USE_SSE
  kappa = expbf(-c1*c2*d2*deltai);
#else
  kappa = exp(-c1*c2*d2*deltai);
#endif
  gvol = p*kappa*pow(pi*deltai,1.5);

  return gvol;
}

/* a switching function for the inverse born radius (beta)
   so that if beta is negative -> beta' = minbeta
   and otherwise beta' = beta^3/(beta^2+a^2) + minbeta
*/ 
 float_a swf_invbr(float_a beta, float_a *fp){
  /* the maximum born radius is 50.0 Ang. */
  static const float_a a  = 0.02;
  static const float_a a2 = 0.02*0.02;
  float_a t;

  if(beta<0.0){
    *fp = 0.0;
    return a;
  }
  t = sqrt(a2 + beta*beta);

  *fp  = beta/t;
  return t;
}

/* a macro to calculate 3 b^2/(b+r)^4 */
#define AGBNP_BRW(b,r) ( _agbnp_brw1 = (b) + (r) , _agbnp_brw2 = _agbnp_brw1 * _agbnp_brw1 , _agbnp_brw3 = _agbnp_brw2 * _agbnp_brw2 ,  3.*((b)*(b))/_agbnp_brw3 )

/* utility to copy caller arrays into local arrays */
void agbnp_cpyd(float_a *local, float_i *caller, int n){
  if(sizeof(float_i)==sizeof(float_a)){
    memcpy(local,caller,n*sizeof(float_a));
  }else{
    int i;
    for(i=0;i<n;i++){
      local[i] = caller[i];
    }
  }
}

float_a agbnp_dist(float_a x1, float_a y1, float_a z1, 
		   float_a x2, float_a y2, float_a z2){
  return sqrt((x2-x1)*(x2-x1)+
	      (y2-y1)*(y2-y1)+
	      (z2-z1)*(z2-z1));
}

/* returns an array with the coordinates of the neighbors */
 int agbnp_get_neigh_coords(AGBNPdata *agb, 
				  int iat, NeighList *neigh_list, 
				  float_a *xj, float_a *yj, float_a *zj){
  float_a *xa = agb->x;
  float_a *ya = agb->y;
  float_a *za = agb->z;
  float_i *xs = agb->xs;
  float_i *ys = agb->ys;
  float_i *zs = agb->zs;
  int nneigh, *neigh;
  NeighList *ext_neigh_list = agb->neigh_list;
  NeighVector **data_index;
  NeighVector *pbc_trans;
  int *int2ext;
  int i, jat, jatom;

  nneigh = neigh_list->nne[iat];
  neigh = neigh_list->neighl[iat];

  if(!agb->dopbc){ 

    /* not doing PBC's, just copy coordinates */
    for(i=0;i<nneigh;i++){
      jat = neigh[i];
      xj[i] = xa[jat];
      yj[i] = ya[jat];
      zj[i] = za[jat];
    }

  }else{

    if(agb->nsym == 1){

      

      /* apply PBC's but use standard coordinates */
      data_index = (NeighVector **)neigh_list->data_index[iat];
      for(i=0;i<nneigh;i++){
	pbc_trans = data_index[i];
	jat = neigh[i];
	xj[i] = xa[jat] - pbc_trans->x;
	yj[i] = ya[jat] - pbc_trans->y;
	zj[i] = za[jat] - pbc_trans->z;
      }

    }else{

      /* crystal PBC's, use external coordinate buffers with external 
	 atom indexes*/
      int2ext = ext_neigh_list->int2ext;
      data_index = (NeighVector **)neigh_list->data_index[iat];
      for(i=0;i<nneigh;i++){
        pbc_trans = data_index[i];
        jat = neigh[i];
	jatom = int2ext ? int2ext[jat] : jat;
        xj[i] = xs[jatom] - pbc_trans->x;
        yj[i] = ys[jatom] - pbc_trans->y;
        zj[i] = zs[jatom] - pbc_trans->z;
      }

    }
  }

  return AGBNP_OK;
}


int mymax(int a, int b){
  return a > b ? a : b;
}

/* transpose x vector matrix multiplication. To rotate gradients */
 void rtvec(float_a y[3], float_a rot[3][3], float_a x[3]){
  int i,j;
  for(i=0;i<3;i++){
    y[i] = 0.0;
    for(j=0;j<3;j++){
      y[i] += x[j]*rot[j][i];
    }
  }
}



 int agbnp_scaling_factors(AGBNPdata *agb, AGBworkdata *agbw_h){
  int i,iat;
  int natoms = agb->natoms;
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  AGBworkdata *agbw = agb->agbw;
  float_a *volumep = agbw->volumep;
  float_a *vols = agbw->vols;
  float_a *sp = agbw->sp;
#ifdef _OPENMP
  float_a *volumep_h = agbw_h->volumep;
  float_a *sp_h = agbw_h->sp;
#endif
  int do_frozen = agb->do_frozen;
  int *isvfrozen = agb->agbw->isvfrozen;
  float_a *volumep_constant = agb->agbw->volumep_constant;

#ifdef _OPENMP

#pragma omp single
  {
    memset(volumep,0,natoms*sizeof(float_a));
  }
#pragma omp barrier
#pragma omp critical
  {
    /* reduce volumep from thread accumulators */ 
    for(i=0;i<nheavyat;i++){
      iat = iheavyat[i];
      volumep[iat] += volumep_h[iat];
    }
  }
#endif

#pragma omp barrier
#pragma omp single
  {
    /* add volume because volumep accumulator was started from
       zero (see agbnp_self_volumes) */
    for(iat=0;iat<natoms;iat++){
      volumep[iat] += vols[iat];
    }
    if(do_frozen){ /* set constant self-volumes */ 
      for(i = 0; i < nheavyat; i++){
	iat = iheavyat[i];
	if(isvfrozen[iat]){
	  volumep[iat] = volumep_constant[iat];
	}
      }
    }
    /* compute scaled volume factors */
    for(i=0;i<nheavyat;i++){
      iat = iheavyat[i];
      sp[iat] = volumep[iat]/vols[iat];
    }
  }
#pragma omp barrier

#ifdef _OPENMP
  /* copy scaling factors to threads
     (may not be necessary) */
  memcpy(sp_h,sp,natoms*sizeof(float_a));
#endif

  return AGBNP_OK;
}

 int agbnp_born_radii(AGBNPdata *agb, AGBworkdata *agbw_h){
  int natoms = agb->natoms;
  int iat;
  float_a fp;
  AGBworkdata *agbw = agb->agbw;
  float_a *br = agbw->br;
  float_a *br1 = agbw->br1;
  float_a *r = agb->r;
  float_a *br1_swf_der = agbw->br1_swf_der;
  float_a *brw = agbw->brw;
  float_a *br1_const = agb->agbw->br1_const;
  int do_frozen =agb->do_frozen;
#ifdef _OPENMP
  float_a *brw_h = agbw_h->brw;
  float_a *br1_swf_der_h = agbw_h->br1_swf_der;
  float_a *br1_h = agbw_h->br1;
  float_a *br_h = agbw_h->br;
#endif
  static const float_a rw = 1.4;  /* water radius offset for np 
				    energy function */
  float_a _agbnp_brw1, _agbnp_brw2, _agbnp_brw3; 

#ifdef _OPENMP
  /* if multiple threads reduce inverse born radii. br1 is master copy. */
#pragma omp single
  {
    for(iat=0;iat<natoms;iat++){
      br1[iat] = 1./r[iat];
    }
  }
#pragma omp barrier
#pragma omp critical
  {
    for(iat=0;iat<natoms;iat++){
      br1[iat] += br1_h[iat];
    }
  }
#pragma omp barrier

#else /* _OPENMP */

  /* if no multithreading, master copy has been calculated. Add inverse
     of vdw radius because br1 accumulator started from zero. */
  for(iat=0;iat<natoms;iat++){
    br1[iat] += 1./r[iat];
  }

#endif

#pragma omp single
  if(do_frozen) {/* add constant part to inverse Born radius */
    for(iat=0;iat<natoms;iat++){
      br1[iat] += br1_const[iat]; 
    }
  }


#pragma omp barrier
#pragma omp single
  {
    float_a biat;
    /* br1[], br1_swf_der[], br[], brw[] are all master copies */
    for(iat = 0; iat < natoms ; iat++){
      /* filters out very large or, worse, negative born radii */
      br1[iat] = swf_invbr(br1[iat], &fp);
      /* save derivative of filter function for later use */
      br1_swf_der[iat] = fp;
      /* calculates born radius from inverse born radius */
      br[iat] = 1./br1[iat];
      biat = br[iat];
      brw[iat] = AGBNP_BRW(biat,rw); /* 3*b^2/(b+rw)^4 for np derivative */
    }
  }
#pragma omp barrier

#ifdef _OPENMP
  /* copy to threads local arrays */
  memcpy(br1_h,br1,natoms*sizeof(float_a));
  memcpy(br_h,br,natoms*sizeof(float_a));
  memcpy(br1_swf_der_h,br1_swf_der,natoms*sizeof(float_a));
  memcpy(brw_h,brw,natoms*sizeof(float_a));
#endif

  return AGBNP_OK;
}


 int agbnp_reset_derivatives(AGBNPdata *agb, AGBworkdata *agbw_h){
  int natoms = agb->natoms;

  memset(agbw_h->dgbdr_h,0,3*natoms*sizeof(float_a));
  memset(agbw_h->dvwdr_h,0,3*natoms*sizeof(float_a));
  memset(agbw_h->dera,0,natoms*sizeof(float_a));
  memset(agbw_h->deru,0,natoms*sizeof(float_a));
  memset(agbw_h->derv,0,natoms*sizeof(float_a));
#ifdef _OPENMP
#pragma omp single
  {
    memset(agb->agbw->dgbdr_h,0,3*natoms*sizeof(float_a));
    memset(agb->agbw->dvwdr_h,0,3*natoms*sizeof(float_a));
    memset(agb->agbw->dera,0,natoms*sizeof(float_a));
    memset(agb->agbw->deru,0,natoms*sizeof(float_a));
    memset(agb->agbw->derv,0,natoms*sizeof(float_a));

    /* initialize constant terms */
    if(agb->agbw->dgbdr_c) 
      memset(agb->agbw->dgbdr_c,0,3*natoms*sizeof(float_a));
    if(agb->agbw->dvwdr_c)
      memset(agb->agbw->dgbdr_c,0,3*natoms*sizeof(float_a));
    if(agb->agbw->dera_const)
      memset(agb->agbw->dera_const,0,natoms*sizeof(float_a));
    if(agb->agbw->deru_c)
      memset(agb->agbw->deru_c,0,natoms*sizeof(float_a));
    if(agb->agbw->derv_c)
      memset(agb->agbw->derv_c,0,natoms*sizeof(float_a));

  }
#pragma omp barrier
#endif

  return AGBNP_OK;
}




/* Evaluates Ui's and Vi's */
 int agbnp_gb_deruv(AGBNPdata *agb, AGBworkdata *agbw_h, 
			  int init_frozen){
  int iq4cache=0;
  int i, iat, j, jat;
  float_a q;

  int natoms = agb->natoms;
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  int *isheavy = agbw_h->isheavy;
  int nhydrogen = agb->nhydrogen;
  int *ihydrogen = agb->ihydrogen;
  NeighList *near_nl = agbw_h->near_nl;
  NeighList *far_nl = agbw_h->far_nl;
  float_a *q2ab = agbw_h->q2ab;
  float_a *abrw = agbw_h->abrw;
  float_a *deru = agbw_h->deru;
  float_a *derv = agbw_h->derv;
  float *q4cache = agbw_h->q4cache;
  float_a *deru_m = agb->agbw->deru;
  float_a *derv_m = agb->agbw->derv;
  float_a *vols = agbw_h->vols;
  int do_frozen = agb->do_frozen;
  int *isvfrozen = agb->agbw->isvfrozen;
  float_a *deru_c = agb->agbw->deru_c;
  float_a *derv_c = agb->agbw->derv_c;

  /* Loop over near heavy atom neighbors */
  for(i=0;i<nheavyat;i++){
    iat = iheavyat[i];
    for(j=0;j<near_nl->nne[iat];j++){
      jat = near_nl->neighl[iat][j] % natoms;
      if(do_frozen){
	/* skip vfrozen pairs */
	if(isvfrozen[iat] && isvfrozen[jat]) continue;
      }
      /* get from cache */
      q = q4cache[iq4cache++];
      iq4cache += 1;
      deru[iat] += q2ab[jat]*q;
      derv[iat] += abrw[jat]*q;
      if(init_frozen){
	if(isvfrozen[iat] && isvfrozen[jat]){
#pragma omp critical
	  {
	    deru_c[iat] += q2ab[jat]*q;
	    derv_c[iat] += abrw[jat]*q;
	  }
	}
      }
      /* get from cache */
      q = q4cache[iq4cache++];
      iq4cache += 1;
      deru[jat] += q2ab[iat]*q;
      derv[jat] += abrw[iat]*q;
      if(init_frozen){
	if(isvfrozen[iat] && isvfrozen[jat]){
#pragma omp critical
	  {
	     deru_c[jat] += q2ab[iat]*q;
	     derv_c[jat] += abrw[iat]*q;
	  }
	}
      }

    }
  }
  /* loop over heavy atoms far neighbors */
  for(i=0;i<nheavyat;i++){
    iat = iheavyat[i];
    for(j=0;j<far_nl->nne[iat];j++){
      jat = far_nl->neighl[iat][j] % natoms;
      if(do_frozen){
	if(isvfrozen[iat] && isvfrozen[jat]) continue;
      }
      /* get from cache */
      if(isheavy[jat]){
	q = q4cache[iq4cache++];
	iq4cache += 1;
	deru[iat] += q2ab[jat]*q;
	derv[iat] += abrw[jat]*q;
	if(init_frozen){
	  if(isvfrozen[iat] && isvfrozen[jat]){
#pragma omp critical
	    {
	      deru_c[iat] += q2ab[jat]*q;
	      derv_c[iat] += abrw[jat]*q;
	    }
	  }
	}
	/* get from cache */
	q = q4cache[iq4cache++];
	iq4cache += 1;
	deru[jat] += q2ab[iat]*q;
	derv[jat] += abrw[iat]*q;
	if(init_frozen){
	  if(isvfrozen[iat] && isvfrozen[jat]){
#pragma omp critical
	    {
	      deru_c[jat] += q2ab[iat]*q;
	      derv_c[jat] += abrw[iat]*q;
	    }
	  }
	}
      }else{
	q = q4cache[iq4cache++];
	iq4cache += 1;
	deru[iat] += q2ab[jat]*q;
	derv[iat] += abrw[jat]*q;
	if(init_frozen){
	  if(isvfrozen[iat] && isvfrozen[jat]){
#pragma omp critical
	    {
	      deru_c[jat] += q2ab[iat]*q;
	      derv_c[jat] += abrw[iat]*q;
	    }
	  }
	}
      }
    }
  }  
  for(i=0;i<nhydrogen;i++){ /* loop for hydrogen-heavy interactions */
    iat = ihydrogen[i];
    for(j=0;j< far_nl->nne[iat];j++){
      jat = far_nl->neighl[iat][j] % natoms;
      if(do_frozen){
	if(isvfrozen[jat] && isvfrozen[iat]) continue;
      }
      /* get from cache */
      if(isheavy[jat]){
	q = q4cache[iq4cache++];
	iq4cache += 1;
	deru[jat] += q2ab[iat]*q;
	derv[jat] += abrw[iat]*q;
	if(init_frozen){
	  if(isvfrozen[iat] && isvfrozen[jat]){
#pragma omp critical
	    {
	      deru_c[jat] += q2ab[iat]*q;
	      derv_c[jat] += abrw[iat]*q;
	    }
	  }
	}
      }
    }
  }

#ifdef _OPENMP
#pragma omp critical
  /* reduction of derv and deru*/
  {
    for(iat=0;iat<natoms;iat++){
      deru_m[iat] += deru[iat];
      derv_m[iat] += derv[iat];
    }
  }
#endif

  if(do_frozen){
    /* add constant portions */
#pragma omp single
    {
      for(iat=0;iat<natoms;iat++){
	deru_m[iat] += deru_c[iat];
	derv_m[iat] += derv_c[iat];
      }
    }
  }

#pragma omp barrier
#pragma omp single
  for(iat=0;iat<natoms;iat++){
    deru_m[iat] /= vols[iat];
    derv_m[iat] /= vols[iat];
  }
#pragma omp barrier

#ifdef _OPENMP
  /* copy reduced quantities to threads */
  memcpy(deru,deru_m,natoms*sizeof(float_a));
  memcpy(derv,derv_m,natoms*sizeof(float_a));
#endif

  return AGBNP_OK;
}

/* Set vfrozen set (the set of frozen heavy atoms with only frozen neighbors)*/
 int agbnp_set_vfrozen(AGBNPdata *agb, AGBworkdata *agbw, 
			     CAVworkdata *cavw){
  int i,j,iat,jat;
  int natoms = agb->natoms;
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  int *isfrozen = agb->isfrozen;
  int *isvfrozen;
  int *isvvfrozen;
  NeighList *near_nl;

  /* work with agb or cavity data depending on which one is defined */
  if(agbw){
    isvfrozen = agb->agbw->isvfrozen;
    isvvfrozen = agb->agbw->isvvfrozen;
    near_nl = agbw->near_nl;
  }else if(cavw){
    isvfrozen = agb->cavw->isvfrozen;
    isvvfrozen = agb->cavw->isvvfrozen;
    near_nl = cavw->near_nl;
  }else{
    fprintf(stderr, "agbnp_set_vfrozen(): error: undefined arguments.\n");
    return AGBNP_ERR;
  }

#pragma omp single
  {
    /* initialize isvfrozen. 
       Start with all frozen atoms in vfrozen set. (Frozen hydrogens 
       are logically part of vfrozen set) */
    memcpy(isvfrozen,isfrozen,natoms*sizeof(int));
  }
#pragma omp barrier

  /* (i,j) i<j is a neighbor pair:
     if i is free j is not vfrozen
     if j is free i is not vfrozen */
  for(i=0;i<nheavyat;i++){
    iat = iheavyat[i];
    if(near_nl->nne[iat]<=0) continue;
    for(j=0;j<near_nl->nne[iat];j++){

      /* near-nl includes only heavy atoms */
      jat = near_nl->neighl[iat][j] % natoms;

      /* if iat is free jat cannot be vfrozen */
      if(!isfrozen[iat]){
#ifdef _OPENMP
#pragma omp critical
#endif
	isvfrozen[jat] = 0;
      }

      /* if jat is free iat cannot be vfrozen */
      if(!isfrozen[jat]){
#ifdef _OPENMP
#pragma omp critical
#endif
	isvfrozen[iat] = 0;
      }

    }
  }

  /* set vvfrozen holds those vfrozen atoms whose near neighbors are all 
     vfrozen */
#pragma omp barrier
#pragma omp single
  {
    memcpy(isvvfrozen,isvfrozen,natoms*sizeof(int));
  }
#pragma omp barrier
  for(i=0;i<nheavyat;i++){
    iat = iheavyat[i];
    if(near_nl->nne[iat]<=0) continue;
    for(j=0;j<near_nl->nne[iat];j++){
      jat = near_nl->neighl[iat][j] % natoms;
      if(!isvfrozen[jat]){
#ifdef _OPENMP
#pragma omp critical
#endif
	isvvfrozen[iat] = 0;
      }
      if(!isvfrozen[iat]){
#ifdef _OPENMP
#pragma omp critical
#endif
	isvvfrozen[jat] = 0;
      }
    }
  }
#pragma omp barrier


  return AGBNP_OK;
}

/* Set bfrozen set (the set of frozen atoms with constan Born radius) */
 int agbnp_set_bfrozen(AGBNPdata *agb, AGBworkdata *agbw){ 
  int i,j,iat,jat;
  int natoms = agb->natoms;
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  int *isbfrozen = agb->agbw->isbfrozen;
  /* assumes that vfrozen set has been already constructed */
  int *isvfrozen = agb->agbw->isvfrozen;
  NeighList *near_nl = agbw->near_nl;
  NeighList *far_nl = agbw->far_nl;
  int *isheavy = agbw->isheavy;

#pragma omp single
  {
    /* initialize isbfrozen. 
       Start with all vfrozen atoms in bfrozen set. (Frozen hydrogens 
       are logically part of vfrozen set) */
    memcpy(isbfrozen,isvfrozen,natoms*sizeof(int));
  }
#pragma omp barrier

  /* (i,j) i<j is a neighbor pair:
     if heavy atom i is not vfrozen j is not bfrozen
     if heavy atom j is not vfrozen i is not bfrozen */

  /* near neighbors first (these are all heavy atoms) */
  for(i=0;i<nheavyat;i++){
    iat = iheavyat[i];
    if(near_nl->nne[iat]<=0) continue;
    for(j=0;j<near_nl->nne[iat];j++){

      /* near-nl includes only heavy atoms */
      jat = near_nl->neighl[iat][j] % natoms;

      /* if iat is not vfrozen jat cannot be bfrozen */
      if(!isvfrozen[iat]){
#ifdef _OPENMP
#pragma omp critical
#endif
	isbfrozen[jat] = 0;
      }

      /* if jat is not vfrozen iat cannot be bfrozen */
      if(!isvfrozen[jat]){
#ifdef _OPENMP
#pragma omp critical
#endif
	isbfrozen[iat] = 0;
      }

    }
  }

  /* Now far neighbors */ 
  for(iat=0;iat<natoms;iat++){
    if(far_nl->nne[iat]<=0) continue;
    for(j=0;j<far_nl->nne[iat];j++){

      jat = far_nl->neighl[iat][j] % natoms;

      /* if iat is a heavy atoms and is not vfrozen jat cannot be bfrozen */
      if(isheavy[iat] && !isvfrozen[iat]){
#ifdef _OPENMP
#pragma omp critical
#endif
	isbfrozen[jat] = 0;
      }

      /* if jat is a heavy atom and is not vfrozen iat cannot be bfrozen */
      if(isheavy[jat] && !isvfrozen[jat]){
#ifdef _OPENMP
#pragma omp critical
#endif
	isbfrozen[iat] = 0;
      }

    }
  }
#pragma omp barrier

  return AGBNP_OK;
}



/* initialize constant terms for generalized Born and vdW calculation */
 int agbnp_init_frozen_agb(AGBNPdata *agb, 
				 float_i *x, float_i *y, float_i *z){

  int i, iat; /* counters */
  int init_frozen = 1;
  float_a egb_self, egb_pair;

  AGBworkdata *agbw = agb->agbw; /* shared work space */
  AGBworkdata *agbw_h;           /* work space for this thread */
#ifdef _OPENMP
  int iproc;
#endif
  int res, error = 0, nop = 0;
  int natoms = agb->natoms;
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  int *isvfrozen = agbw->isvfrozen;
  float_a *volumep = agbw->volumep;
  float_a *volumep_constant = agbw->volumep_constant;

  if(!agb->do_frozen) {
    return AGBNP_OK;
  }

  /* copy coordinates to internal buffers */
  agbnp_cpyd(agb->x,x,agb->natoms);
  agbnp_cpyd(agb->y,y,agb->natoms);
  agbnp_cpyd(agb->z,z,agb->natoms);

  /* now act as if no atoms are frozen */
  agb->do_frozen = 0;

#ifdef _OPENMP
#pragma omp parallel private(agbw_h, iproc, error)
#endif
  {

#ifdef _OPENMP
    /* get local work structure */
    iproc = omp_get_thread_num();
    agbw_h = agb->agbw_p[iproc];
#else
    agbw_h = agbw;
#endif /* _OPENMP */

    error = 0;

    /* constructs neighbor lists */
    if(agb->dopbc || agb->doalt){
      res = agbnp_neighbor_lists_pbc(agb, agbw_h, agb->x, agb->y, agb->z);
    }else{
      res = agbnp_neighbor_lists(agb, agbw_h, agb->x, agb->y, agb->z);
    }
    if(res != AGBNP_OK){
      fprintf(stderr,"agbnp_init_frozen_agb(): error in agbnp_neighbor_lists()\n");
      error = 1;
      goto ERROR;
    }

    /* set vfrozen set (frozen atoms which have only frozen near neighbors) */
    if(agbnp_set_vfrozen(agb, agbw_h, NULL) != AGBNP_OK){
      fprintf(stderr,"agbnp_init_frozen_agb(): error in agbnp_set_vfrozen()\n");
      error = 1;
      goto ERROR;
    }

    /* set bfrozen set (vfrozen atoms with constant Born radius, with only
       vfrozen neighbors) */
    if(agbnp_set_bfrozen(agb, agbw_h) != AGBNP_OK){
      fprintf(stderr,"agbnp_init_frozen_agb(): error in agbnp_set_bfrozen()\n");
      error = 1;
      goto ERROR;
    }

    /* computes self volumes  as  if no atoms are frozen */
    if(agb->dopbc || agb->doalt){
      res = agbnp_self_volumes_pbc(agb, agbw_h, agb->x, agb->y, agb->z);
    }else{
      res = agbnp_self_volumes(agb, agbw_h, agb->x, agb->y, agb->z);
    }
    if(res != AGBNP_OK){
      fprintf(stderr, "agbnp_init_frozen_agb(): error in agbnp_self_volumes()\n");
      error = 1;
      goto ERROR;
    }

    /* calculate volume scaling factors */
    if(agbnp_scaling_factors(agb, agbw_h) != AGBNP_OK){
      fprintf(stderr, "agbnp_init_frozen_agb(): error in agbnp_scaling_factors()\n");
      error = 1;
      goto ERROR;
    }

    /* assign constant self_volumes */
#pragma omp barrier
#pragma omp single
    {
      /* set constant self volumes */
      memset(volumep_constant,0,natoms*sizeof(float_a));
      for(i = 0; i < nheavyat; i++){
	iat = iheavyat[i];
	if(isvfrozen[iat]){
	  volumep_constant[iat] = volumep[iat];
	}
      }
    }
#pragma omp barrier

    /* calculates inverse Born radii, saves constant portion (init_frozen) */
    if(agb->dopbc || agb->doalt){
      res = agbnp_inverse_born_radii_pbc(agb, agbw_h, agb->x, agb->y, agb->z, init_frozen);
    }else{
      res = agbnp_inverse_born_radii(agb, agbw_h, agb->x, agb->y, agb->z, init_frozen);
    }
    if(res != AGBNP_OK){
      fprintf(stderr, "agbnp_init_frozen_agb(): error in agbnp_inverse_born_radii()\n");
      error = 1;
      goto ERROR;
    }

    /* calculates Born radii and related quantities (brw, etc.) */
    if(agbnp_born_radii(agb, agbw_h) != AGBNP_OK){
      fprintf(stderr, "agbnp_init_frozen_agb(): error in agbnp_born_radii()\n");
      error = 1;
      goto ERROR;
    }

    if(agbnp_reset_derivatives(agb,agbw_h) != AGBNP_OK){
      fprintf(stderr, "agbnp_init_frozen_agb(): error in agbnp_reset_derivatives\n");
      error = 1;
      goto ERROR;
    }

    /* calculates GB energy, saves constant portions (init_frozen) */
    if(agb->dopbc || agb->doalt){
      res = agbnp_gb_energy_pbc(agb, agbw_h, agb->x, agb->y, agb->z, 
				&egb_self, &egb_pair, init_frozen);
    }else{
      res = agbnp_gb_energy(agb, agbw_h, agb->x, agb->y, agb->z, 
			    &egb_self, &egb_pair, init_frozen); 
    }
    if(res != AGBNP_OK){
      fprintf(stderr, "agbnp_init_frozen_agb(): error in agbnp_gb_energy()\n");
      error = 1;
      goto ERROR;
    }

    /* GB derivatives contribution at constant self volumes */
    if(agb->dopbc || agb->doalt){
      res = agbnp_gb_ders_constvp_pbc(agb, agbw_h, agb->x, agb->y, agb->z, 
				      init_frozen);
    }else{
      res = agbnp_gb_ders_constvp(agb, agbw_h, agb->x, agb->y, agb->z, 
				  init_frozen);
    }
    if(res != AGBNP_OK){
      fprintf(stderr, "agbnp_init_frozen_agb(): error in agbnp_gb_ders_constvp()\n");
      error = 1;
      goto ERROR;
    }

    /* evaluation of Ui's and Vi's*/
    if(agbnp_gb_deruv(agb, agbw_h, init_frozen) != AGBNP_OK){
      fprintf(stderr, "agbnp_init_frozen_agb(): error in agbnp_gb_deruvp()\n");
      error = 1;
      goto ERROR;
    }

    /* derivatives due to changes in self volumes */
    if(agb->dopbc || agb->doalt){
      res = agbnp_der_vp_pbc(agb, agbw_h, agb->x, agb->y, agb->z, init_frozen);
    }else{
      res = agbnp_der_vp(agb, agbw_h, agb->x, agb->y, agb->z, init_frozen);
    }
    if(res != AGBNP_OK){
      fprintf(stderr, "agbnp_init_frozen_agb(): error in agbnp_der_vp()\n");
      error = 1;
      goto ERROR;
    }

  ERROR:
    nop;

  }/* #pragma omp parallel */

  if(error){
    return AGBNP_ERR;
  }


  /* puts frozen flag back up */
  agb->do_frozen = 1;

  return AGBNP_OK;
}


 int agb_vdw_energy(AGBNPdata *agb, float_i *x, float_i *y, float_i *z,
		    float_i *sp, float_i *br, float_i *egb, float_i (*dgbdr)[3],
		    float_i *evdw, float_i (*dvwdr)[3],
		    float_i *ecorr){

  static const float_a tokcalmol = 332.0; /* conversion to kcal/mol */
  int i, ki; /* counters */
  float_a a;
  int iat; /* atomic counter */
  static const float_a rw = 1.4;  /* water radius offset for np 
				    energy function */
  int no_init_frozen = 0;
  float_a dielectric_factor = 
    -0.5*(1./agb->dielectric_in - 1./agb->dielectric_out);

  int natoms = agb->natoms;
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  float_a *ialpha = agb->ialpha;
  float_a *salpha = agb->salpha;
  float_a *idelta = agb->idelta;
  float_a *sdelta = agb->sdelta;
  float_a egb_self, egb_pair;
  
  AGBworkdata *agbw = agb->agbw; /* shared work space */
  AGBworkdata *agbw_h;           /* work space for this thread */
#ifdef _OPENMP
  int iproc;
#endif
  int res, error = 0, nop = 0;

#ifdef MMTIMER
  static MMtimer timer1 = -1, timer2 = -1, timer3 = -1;
  static MMtimer timer4 = -1, timer5 = -1, timer6 = -1, timer7 = -1;
  double user_time, system_time, real_time;
#endif

#ifdef MMTIMER
  if(agbvdw_timer<0) mmtimer_new( &agbvdw_timer );
  if(agbvdw_timer<0) {
    fprintf(stderr,"agb_vdw_energy(): error initializing agbvdw_timer\n");
    error = 1;
    goto ERROR;
  }
  mmtimer_start( agbvdw_timer );
  if(timer1<0) mmtimer_new( &timer1);
  if(timer2<0) mmtimer_new( &timer2);
  if(timer3<0) mmtimer_new( &timer3);
  if(timer4<0) mmtimer_new( &timer4);
  if(timer5<0) mmtimer_new( &timer5);
  if(timer6<0) mmtimer_new( &timer6);
  if(timer7<0) mmtimer_new( &timer7);
  if(timer1<0 || timer2<0 || timer3<0 || timer4<0 ||
     timer5<0 || timer6<0 || timer7<0){
    fprintf(stderr,"agb_vdw_energy(): error initializing timers\n");
    error = 1;
    goto ERROR;
  }
#endif

  /* copy coordinates in local arrays */
  agbnp_cpyd(agb->x,x,agb->natoms);
  agbnp_cpyd(agb->y,y,agb->natoms);
  agbnp_cpyd(agb->z,z,agb->natoms);

#ifdef _OPENMP
#pragma omp parallel private(agbw_h, iproc, error)
#endif
  {

#ifdef _OPENMP
  /* get local work structure */
  iproc = omp_get_thread_num();
  agbw_h = agb->agbw_p[iproc];
#else
  agbw_h = agbw;
#endif /* _OPENMP */

  error = 0;


#ifdef MMTIMER
  mmtimer_start( timer1 );
#endif
  /* constructs neighbor lists */
  if(agb->dopbc || agb->doalt){
    res = agbnp_neighbor_lists_pbc(agb, agbw_h, agb->x, agb->y, agb->z);
  }else{
    res = agbnp_neighbor_lists(agb, agbw_h, agb->x, agb->y, agb->z);
  }
  if(res != AGBNP_OK){
    fprintf(stderr,"agb_vdw_energy(): error in agbnp_neighbor_lists()\n");
    error = 1;
    goto ERROR;
  }
#ifdef MMTIMER
  mmtimer_stop(timer1);
#endif
#ifdef MMTIMER
  mmtimer_start( timer2 );
#endif
  /* computes self volumes  */
  if(agb->dopbc || agb->doalt){
    res = agbnp_self_volumes_pbc(agb, agbw_h, agb->x, agb->y, agb->z);
  }else{
    res = agbnp_self_volumes(agb, agbw_h, agb->x, agb->y, agb->z);
  }
  if(res != AGBNP_OK){
    fprintf(stderr, "agb_vdw_energy(): error in agbnp_self_volumes()\n");
    error = 1;
    goto ERROR;
  }

  /* calculate volume scaling factors */
  if(agbnp_scaling_factors(agb, agbw_h) != AGBNP_OK){
    fprintf(stderr, "agb_vdw_energy(): error in agbnp_scaling_factors()\n");
    error = 1;
    goto ERROR;
  }
#ifdef MMTIMER
  mmtimer_stop(timer2);
#endif
#ifdef MMTIMER
  mmtimer_start( timer3 );
#endif
  /* calculates inverse Born radii */
  if(agb->dopbc || agb->doalt){
    res = agbnp_inverse_born_radii_pbc(agb, agbw_h, agb->x, agb->y, agb->z,
				       no_init_frozen); 
  }else{
    res = agbnp_inverse_born_radii(agb, agbw_h, agb->x, agb->y, agb->z, 
				   no_init_frozen);
  }
  if(res != AGBNP_OK){
    fprintf(stderr, "agb_vdw_energy(): error in agbnp_inverse_born_radii()\n");
    error = 1;
    goto ERROR;
  }

  /* calculates Born radii and related quantities (brw, etc.) */
  if(agbnp_born_radii(agb, agbw_h) != AGBNP_OK){
    fprintf(stderr, "agb_vdw_energy(): error in agbnp_born_radii()\n");
    error = 1;
    goto ERROR;
  }
#ifdef MMTIMER
 mmtimer_stop( timer3 );
#endif

#pragma omp barrier
#pragma omp single
  /* calculates van der waals energy */
  {
    for(iat=0;iat<natoms;iat++){
      br[iat] = agb->agbw->br[iat];
    }
    *evdw = 0.0;
    *ecorr = 0.0;
    for(iat=0;iat<natoms;iat++){
      a = 1.0/(br[iat]+rw);
      a = pow(a,3);
      *evdw += agb->occupancy[iat]*(ialpha[iat]*a + idelta[iat]);
      *ecorr += agb->occupancy[iat]*(salpha[iat]*a + sdelta[iat]);
    }
  }

  if(agbnp_reset_derivatives(agb,agbw_h) != AGBNP_OK){
    fprintf(stderr, "agb_vdw_energy(): error in agbnp_reset_derivatives\n");
    error = 1;
    goto ERROR;
  }
#ifdef MMTIMER
  mmtimer_start( timer4 );
#endif
  /* Evaluates solvation energy, Ai's and derivatives of GB energy 
     at constant Born radii */
  if(agb->dopbc || agb->doalt){ 
    res = agbnp_gb_energy_pbc(agb, agbw_h, agb->x, agb->y, agb->z, 
			      &egb_self, &egb_pair, no_init_frozen);
  }else{
    res = agbnp_gb_energy(agb, agbw_h, agb->x, agb->y, agb->z, 
			  &egb_self, &egb_pair, no_init_frozen);
  }
  if(res != AGBNP_OK){
    fprintf(stderr, "agb_vdw_energy(): error in agbnp_gb_energy()\n");
    error = 1;
    goto ERROR;
  }

#pragma omp barrier
#pragma omp single
  {
    /* calculates GB energy */
    egb_self *=  tokcalmol*dielectric_factor;
    egb_pair *= tokcalmol*dielectric_factor;
    *egb = egb_self + egb_pair;
  }
#ifdef MMTIMER
 mmtimer_stop( timer4 );
#endif

#ifdef MMTIMER
  mmtimer_start( timer5 );
#endif
  /* GB derivatives contribution at constant self volumes */
  if(agb->dopbc || agb->doalt){ 
    res = agbnp_gb_ders_constvp_pbc(agb, agbw_h, agb->x, agb->y, agb->z, 
				    no_init_frozen);
  }else{
    res = agbnp_gb_ders_constvp(agb, agbw_h, agb->x, agb->y, agb->z, 
				no_init_frozen);
  }
  if(res != AGBNP_OK){
    fprintf(stderr, "agb_vdw_energy(): error in agbnp_gb_ders_constvp()\n");
    error = 1;
    goto ERROR;
  }
#ifdef MMTIMER
 mmtimer_stop( timer5 );
#endif

#ifdef MMTIMER
  mmtimer_start( timer6 );
#endif
  /* evaluation of Ui's and Vi's*/
  if(agbnp_gb_deruv(agb, agbw_h, no_init_frozen) != AGBNP_OK){
    fprintf(stderr, "agb_vdw_energy(): error in agbnp_gb_deruvp()\n");
    error = 1;
    goto ERROR;
  }
#ifdef MMTIMER
 mmtimer_stop( timer6 );
#endif
#ifdef MMTIMER
  mmtimer_start( timer7 );
#endif

  /* derivatives due to changes in self volumes */
  if(agb->dopbc || agb->doalt){ 
    res = agbnp_der_vp_pbc(agb, agbw_h, agb->x, agb->y, agb->z, 
			   no_init_frozen);
  }else{
    res = agbnp_der_vp(agb, agbw_h, agb->x, agb->y, agb->z, 
		       no_init_frozen);
  }
  if(res != AGBNP_OK){
    fprintf(stderr, "agb_vdw_energy(): error in agbnp_der_vp()\n");
    error = 1;
    goto ERROR;
  } 


#ifdef MMTIMER
 mmtimer_stop( timer7 );
#endif

#ifdef _OPENMP
#pragma omp critical
  /* reduction of derivatives */
  {
    for(iat=0;iat<natoms;iat++){
      for(i=0;i<3;i++){
	agbw->dgbdr_h[iat][i] += agbw_h->dgbdr_h[iat][i];
	agbw->dvwdr_h[iat][i] += agbw_h->dvwdr_h[iat][i];
      }
    }
  }
#endif
  
  if(agb->do_frozen){
    /* add constant terms to derivatives */
#pragma omp single
    {
      float_a (*dgbdr_c)[3] = agb->agbw->dgbdr_c;
      float_a (*dvwdr_c)[3] = agb->agbw->dvwdr_c;
      for(iat=0;iat<natoms;iat++){
	for(i=0;i<3;i++){
	  agbw->dgbdr_h[iat][i] += dgbdr_c[iat][i];
	  agbw->dvwdr_h[iat][i] += dvwdr_c[iat][i];
	}
      }
    }
  }

  ERROR:
  nop;

  }/* #pragma omp parallel */

  if(error){
    return AGBNP_ERR;
  }

  /* return derivatives */
  for(iat=0;iat<natoms;iat++){
    for(ki=0;ki<3;ki++){
      dgbdr[iat][ki] = tokcalmol*agbw->dgbdr_h[iat][ki];
      dvwdr[iat][ki] = agbw->dvwdr_h[iat][ki];
    }
  }

  /* returns scaled volume factors */
  for(iat=0;iat<natoms;iat++){
    sp[iat] = 1.0;
  }
  for(i=0;i<nheavyat;i++){
    iat = iheavyat[i];
    sp[iat] = agbw->sp[iat];
  }

  /* return born radii */
  for(iat=0;iat<natoms;iat++){
    br[iat] = agbw->br[iat];
  }


#ifdef MMTIMER
 mmtimer_stop( agbvdw_timer );
 printf("\n");
 mmtimer_report( timer1, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
	"AGB NB", user_time, system_time, real_time);
 mmtimer_report( timer2, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
        "AGB Self. Vol.", user_time, system_time, real_time);
 mmtimer_report( timer3, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
        "AGB Born radii", user_time, system_time, real_time);
 mmtimer_report( timer4, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
	"AGB energy", user_time, system_time, real_time);
 mmtimer_report( timer5, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
        "AGB der1", user_time, system_time, real_time);
 mmtimer_report( timer6, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
        "AGB U/V", user_time, system_time, real_time);
 mmtimer_report( timer7, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
        "AGB der2", user_time, system_time, real_time);
 mmtimer_report( agbvdw_timer, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
        "AGB+VDW tot timer:", user_time, system_time, real_time);
#endif

  return AGBNP_OK;
}


int plop_agb_vdw_energy(AGBNPdata *agb, int natom_asu, float_i *x, float_i *y, float_i *z,
                        float_i *sp, float_i *br, float_i *evdw, float_i *ecorr){
  
  static const float_a tokcalmol = 332.0; /* conversion to kcal/mol */
  int i, ki; /* counters */
  float_a a;
  int iat; /* atomic counter */
  static const float_a rw = 1.4;  /* water radius offset for np 
				    energy function */
  int no_init_frozen = 0;
  float_a dielectric_factor = 
    -0.5*(1./agb->dielectric_in - 1./agb->dielectric_out);

  int natoms = agb->natoms;
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  float_a *ialpha = agb->ialpha;
  float_a *salpha = agb->salpha;
  float_a *idelta = agb->idelta;
  float_a *sdelta = agb->sdelta;
  float_a egb_self, egb_pair;
  
  AGBworkdata *agbw = agb->agbw; /* shared work space */
  AGBworkdata *agbw_h;           /* work space for this thread */
#ifdef _OPENMP
  int iproc;
#endif
  int res, error = 0, nop = 0;

#ifdef MMTIMER
  static MMtimer timer1 = -1, timer2 = -1, timer3 = -1;
  double user_time, system_time, real_time;
#endif

#ifdef MMTIMER
  if(agbvdw_timer<0) mmtimer_new( &agbvdw_timer );
  if(agbvdw_timer<0) {
    fprintf(stderr,"plop_agb_vdw_energy(): error initializing agbvdw_timer\n");
    error = 1;
    goto ERROR;
  }
  mmtimer_start( agbvdw_timer );
  if(timer1<0) mmtimer_new( &timer1);
  if(timer2<0) mmtimer_new( &timer2);
  if(timer3<0) mmtimer_new( &timer3);


  if(timer1<0 || timer2<0 || timer3<0) {
    fprintf(stderr,"plop_agb_vdw_energy(): error initializing timers\n");
    error = 1;
    goto ERROR;
  }
#endif

  /* copy coordinates in local arrays */
  agbnp_cpyd(agb->x,x,agb->natoms);
  agbnp_cpyd(agb->y,y,agb->natoms);
  agbnp_cpyd(agb->z,z,agb->natoms);

#ifdef _OPENMP
#pragma omp parallel private(agbw_h, iproc, error)
#endif
  {

#ifdef _OPENMP
  /* get local work structure */
  iproc = omp_get_thread_num();
  agbw_h = agb->agbw_p[iproc];
#else
  agbw_h = agbw;
#endif /* _OPENMP */

  error = 0;


#ifdef MMTIMER
  mmtimer_start( timer1 );
#endif
  /* constructs neighbor lists */
  if(agb->dopbc || agb->doalt){
    res = agbnp_neighbor_lists_pbc(agb, agbw_h, agb->x, agb->y, agb->z);
  }else{
    res = agbnp_neighbor_lists(agb, agbw_h, agb->x, agb->y, agb->z);
  }
  if(res != AGBNP_OK){
    fprintf(stderr,"agb_vdw_energy(): error in agbnp_neighbor_lists()\n");
    error = 1;
    goto ERROR;
  }
#ifdef MMTIMER
  mmtimer_stop(timer1);
#endif
#ifdef MMTIMER
  mmtimer_start( timer2 );
#endif
  /* computes self volumes  */
  if(agb->dopbc || agb->doalt){
    res = agbnp_self_volumes_pbc(agb, agbw_h, agb->x, agb->y, agb->z);
  }else{
    res = agbnp_self_volumes(agb, agbw_h, agb->x, agb->y, agb->z);
  }
  if(res != AGBNP_OK){
    fprintf(stderr, "agb_vdw_energy(): error in agbnp_self_volumes()\n");
    error = 1;
    goto ERROR;
  }

  /* calculate volume scaling factors */
  if(agbnp_scaling_factors(agb, agbw_h) != AGBNP_OK){
    fprintf(stderr, "agb_vdw_energy(): error in agbnp_scaling_factors()\n");
    error = 1;
    goto ERROR;
  }
#ifdef MMTIMER
  mmtimer_stop(timer2);
#endif
#ifdef MMTIMER
  mmtimer_start( timer3 );
#endif
  /* calculates inverse Born radii */
  if(agb->dopbc || agb->doalt){
    res = agbnp_inverse_born_radii_pbc(agb, agbw_h, agb->x, agb->y, agb->z,
				       no_init_frozen); 
  }else{
    res = agbnp_inverse_born_radii(agb, agbw_h, agb->x, agb->y, agb->z, 
				   no_init_frozen);
  }
  if(res != AGBNP_OK){
    fprintf(stderr, "agb_vdw_energy(): error in agbnp_inverse_born_radii()\n");
    error = 1;
    goto ERROR;
  }

  /* calculates Born radii and related quantities (brw, etc.) */
  if(agbnp_born_radii(agb, agbw_h) != AGBNP_OK){
    fprintf(stderr, "agb_vdw_energy(): error in agbnp_born_radii()\n");
    error = 1;
    goto ERROR;
  }
#ifdef MMTIMER
 mmtimer_stop( timer3 );
#endif

#pragma omp barrier
#pragma omp single
  /* calculates van der waals energy */
  {
    for(iat=0;iat<natoms;iat++){
      br[iat] = agb->agbw->br[iat];
    }
    *evdw = 0.0;
    *ecorr = 0.0;
    for(iat=0;iat<natoms;iat++){
      if(iat >= natom_asu) continue;
      a = 1.0/(br[iat]+rw);
      a = pow(a,3);
      *evdw += agb->occupancy[iat]*(ialpha[iat]*a + idelta[iat]);
      *ecorr += agb->occupancy[iat]*(salpha[iat]*a + sdelta[iat]);
    }
  }

  ERROR:
  nop;

  }/* #pragma omp parallel */

  if(error){
    return AGBNP_ERR;
  }

  /* returns scaled volume factors */
  for(iat=0;iat<natoms;iat++){
    sp[iat] = 1.0;
  }
  for(i=0;i<nheavyat;i++){
    iat = iheavyat[i];
    sp[iat] = agbw->sp[iat];
  }

  /* return born radii */
  for(iat=0;iat<natoms;iat++){
    br[iat] = agbw->br[iat];
  }

#ifdef MMTIMER
 mmtimer_stop( agbvdw_timer );
 printf("\n");
 mmtimer_report( timer1, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
	"AGB NB", user_time, system_time, real_time);
 mmtimer_report( timer2, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
        "AGB Self. Vol.", user_time, system_time, real_time);
 mmtimer_report( timer3, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
        "AGB Born radii", user_time, system_time, real_time);
#endif

  return AGBNP_OK;
}


int agbnp_init_frozen_cavity(AGBNPdata *agb, float_i *x, float_i *y, float_i *z){
  int natoms = agb->natoms;
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  float_a *gamma = agb->cavw->gamma;

  CAVworkdata *cavw = agb->cavw;
  CAVworkdata *cavw_h;
  float_a *volumep = cavw->volumep;
  float_a *volumep_constant = cavw->volumep_constant;
  float_a *surfarea_constant = cavw->surfarea_constant;
  int *isvfrozen = agb->cavw->isvfrozen;
  float_a *gammap = agb->cavw->gammap;
  float_a (*decav)[3] = agb->cavw->decav_h;
  float_a (*decav_const)[3] = agb->cavw->decav_const;

#ifdef _OPENMP
  int iproc;
#endif
  int error = 0, nop = 0;
  int i, iat;
  float_a s, a, f, fp;
  int res;

  if(!agb->do_frozen) {
    return AGBNP_OK;
  }

  /* copy coordinates to internal buffers */
  agbnp_cpyd(agb->x,x,agb->natoms);
  agbnp_cpyd(agb->y,y,agb->natoms);
  agbnp_cpyd(agb->z,z,agb->natoms);

  /* now act as if no atoms are frozen */
  agb->do_frozen = 0;

#ifdef _OPENMP
#pragma omp parallel private(cavw_h,iproc,error)
#endif
  {

#ifdef _OPENMP
    /* get local work structure */
    iproc = omp_get_thread_num();
    cavw_h = agb->cavw_p[iproc];
#else
    cavw_h = cavw;
#endif /* _OPENMP */

    error = 0;

    /* constructs neighbor list */
    if(agb->dopbc || agb->doalt){
      res = agbnp_neighbor_list_cavity_pbc(agb, cavw_h, agb->x, agb->y, agb->z);
    }else{
      res = agbnp_neighbor_list_cavity(agb, cavw_h, agb->x, agb->y, agb->z);
    }
    if(res != AGBNP_OK){
      fprintf(stderr,"agbnp_init_frozen_cavity(): error in agbnp_neighbor_list_cavity()\n");
      error = 1;
      goto ERROR;
    }

    /* set vfrozen set (frozen atoms which have only frozen neighbors) */
    if(agbnp_set_vfrozen(agb, NULL, cavw_h) != AGBNP_OK){
      fprintf(stderr,"agbnp_init_frozen_cavity(): error in agbnp_set_vfrozen()\n");
      error = 1;
      goto ERROR;
    }

    /* calculate self-volumes and surface areas as if no atoms are frozen */
    if(agb->dopbc || agb->doalt){
      res = agbnp_surface_areas_pbc(agb, cavw_h, agb->x, agb->y, agb->z);
    }else{
      res = agbnp_surface_areas(agb, cavw_h, agb->x, agb->y, agb->z);
    }
    if(res != AGBNP_OK){
      fprintf(stderr,"agbnp_init_frozen_cavity(): error in agbnp_surface_areas()\n");
      error = 1;
      goto ERROR;
    }
    
    /* assign constant surface areas and self_volumes */
#pragma omp barrier
#pragma omp single
    {
      /* constant self volumes */
      memset(volumep_constant,0,natoms*sizeof(float_a));
      for(i = 0; i < nheavyat; i++){
	iat = iheavyat[i];
	if(isvfrozen[iat]){
	  volumep_constant[iat] = volumep[iat];
	}
      }
      /* constant surface areas */
      memset(surfarea_constant,0,natoms*sizeof(float_a));
      for(i=0;i<nheavyat;i++){
	iat = iheavyat[i];
	if(isvfrozen[iat]){
	  surfarea_constant[iat] = cavw->surf_area[iat];
	}
	a = cavw->surf_area[iat]/agb->occupancy[iat];
	f = swf_area(a, &fp);
	s = cavw->surf_area[iat]*f;
	gammap[iat] = gamma[iat]*(f+a*fp);
      }
    }
#pragma omp barrier

    /* calculate derivatives */
    if(agb->dopbc || agb->doalt){
      res = agbnp_cavity_ders_pbc(agb, cavw_h, agb->x, agb->y, agb->z);
    }else{
      res = agbnp_cavity_ders(agb, cavw_h, agb->x, agb->y, agb->z);
    }
    if(res != AGBNP_OK){
      fprintf(stderr,"agbnp_cavity_energy(): error in agbnp_cavity_ders()\n");
      error = 1;
      goto ERROR;
    }

#pragma omp barrier
#pragma omp single
    {
      /* constant cavity energy gradients */
      memset(decav_const,0,3*natoms*sizeof(float_a));
      for(i=0;i<nheavyat;i++){
	iat = iheavyat[i];
	if(isvfrozen[iat]){
	  decav_const[iat][0] = decav[iat][0];
	  decav_const[iat][1] = decav[iat][1];
	  decav_const[iat][2] = decav[iat][2];
	}
      }
    }

  ERROR:
    nop;

  }

  if(error){
    return AGBNP_ERR;
  }

  /* puts frozen flag back up */
  agb->do_frozen = 1;
 
 return AGBNP_OK;
}


 int agb_cavity_energy(AGBNPdata *agb, float_i *x, float_i *y, float_i *z, 
		       float_i *mol_volume, float_i *surf_area,
		       float_i *ecav, float_i *ecorr, float_i (*decav)[3]){

  int i, ki; /* counter */
  int iat; /* atomic counter */
  int res;
  int natoms = agb->natoms;
  float_a *igamma = agb->igamma; /* ideal gamma */
  float_a *gamma = agb->cavw->gamma; /* ideal gamma + corrction gamma */
  float_a *sgamma = agb->sgamma; /* correction gamma */
  float_a *gammap = agb->cavw->gammap; /* gamma corrected by derivative of 
					 surface area switching function 
					 (for derivative calculation) */
  float_a *volumep = agb->cavw->volumep; /* self-volumes */
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  float_a a, f, fp;
  float_a (*decav_h)[3] = agb->cavw->decav_h;
  float_a (*decav_const)[3] = agb->cavw->decav_const;
  CAVworkdata *cavw = agb->cavw;
  CAVworkdata *cavw_h;
#ifdef _OPENMP
  int iproc;
#endif
  int error = 0, nop = 0;
  int *isvfrozen = cavw->isvfrozen;
  float_a *volumep_constant = cavw->volumep_constant;
  float_a *surfarea_constant = cavw->surfarea_constant;

#ifdef MMTIMER
  static MMtimer timer1 = -1, timer2 = -1, timer3 = -1;
  double user_time, system_time, real_time;
#endif

#ifdef MMTIMER
  if(cavity_timer < 0) mmtimer_new( &cavity_timer );
  if(timer1<0) mmtimer_new( &timer1);
  if(timer2<0) mmtimer_new( &timer2);
  if(timer3<0) mmtimer_new( &timer3);
  mmtimer_start( cavity_timer );
#endif

  /* copy coordinates to internal buffers */
  agbnp_cpyd(agb->x,x,agb->natoms);
  agbnp_cpyd(agb->y,y,agb->natoms);
  agbnp_cpyd(agb->z,z,agb->natoms);

#ifdef _OPENMP
#pragma omp parallel private(cavw_h,iproc,error)
#endif
{

#ifdef _OPENMP
  /* get local work structure */
  iproc = omp_get_thread_num();
  cavw_h = agb->cavw_p[iproc];
#else
  cavw_h = cavw;
#endif /* _OPENMP */

#ifdef MMTIMER
  mmtimer_start( timer1 );
#endif

  error = 0;

  /* constructs neighbor list */
  if(agb->dopbc || agb->doalt){
    res = agbnp_neighbor_list_cavity_pbc(agb, cavw_h, agb->x, agb->y, agb->z);
  }else{
    res = agbnp_neighbor_list_cavity(agb, cavw_h, agb->x, agb->y, agb->z);
  }
  if(res != AGBNP_OK){
    fprintf(stderr,"agb_cavity_energy(): error in agbnp_neighbor_list_cavity()\n");
    error = 1;
    goto ERROR;
  }

#ifdef MMTIMER
  mmtimer_stop( timer1 );
#endif

#ifdef MMTIMER
  mmtimer_start( timer2 );
#endif

  /* calculates self volumes and surface areas*/
  if(agb->dopbc || agb->doalt){
    res = agbnp_surface_areas_pbc(agb, cavw_h, agb->x, agb->y, agb->z);
  }else{
    res = agbnp_surface_areas(agb, cavw_h, agb->x, agb->y, agb->z);
  }
  if(res != AGBNP_OK){
    fprintf(stderr,"agb_cavity_energy(): error in agbnp_surface_areas()\n");
    error = 1;
    goto ERROR;
  }

#ifdef MMTIMER
  mmtimer_stop( timer2 );
#endif


#pragma omp barrier
#pragma omp single
  {
    /* volume of molecule */
    *mol_volume = 0.0;
    if(!agb->do_frozen){
      for(i = 0; i < nheavyat; i++){
	iat = iheavyat[i];
	*mol_volume += volumep[iat];
      }
    }else{
      for(i = 0; i < nheavyat; i++){
	iat = iheavyat[i];
	if(isvfrozen[iat]){
	  *mol_volume += volumep_constant[iat];
	}else{
	  *mol_volume += volumep[iat];
	}
      }
    }

    /* filters surface areas to avoid negative surface areas */
    memset(surf_area,0,natoms*sizeof(float_i));
    if(!agb->do_frozen){
      for(i=0;i<nheavyat;i++){
	iat = iheavyat[i];
	/* use occupancy-unweighted surface area for switching function. This
	   way in the limit of identical alternate conformations surface areas
	   are the same as for one conformation */
	a = cavw->surf_area[iat]/agb->occupancy[iat];
	f = swf_area(a, &fp);
	surf_area[iat] = cavw->surf_area[iat]*f;
	gammap[iat] = gamma[iat]*(f+a*fp);
      }
    }else{
      for(i=0;i<nheavyat;i++){
	iat = iheavyat[i];
	if(isvfrozen[iat]){
	  cavw->surf_area[iat] = surfarea_constant[iat];
	}
	a = cavw->surf_area[iat]/agb->occupancy[iat];
	f = swf_area(a, &fp);
	surf_area[iat] = cavw->surf_area[iat]*f;
	gammap[iat] = gamma[iat]*(f+a*fp);
      }
    }
    /* calculates cavity energy */
    *ecav = 0.0;
    *ecorr = 0.0;
    for(i=0;i<nheavyat;i++){
      iat = iheavyat[i];
      *ecav += igamma[iat]*surf_area[iat];
      *ecorr += sgamma[iat]*surf_area[iat];
    }
  }
#pragma omp barrier

#ifdef MMTIMER
  mmtimer_start( timer3 );
#endif

  /* calculate derivatives */
  if(agb->dopbc || agb->doalt){
    res = agbnp_cavity_ders_pbc(agb, cavw_h, agb->x, agb->y, agb->z);
  }else{
    res = agbnp_cavity_ders(agb, cavw_h, agb->x, agb->y, agb->z);
  }
  if(res != AGBNP_OK){
    fprintf(stderr,"agbnp_cavity_energy(): error in agbnp_cavity_ders()\n");
    error = 1;
    goto ERROR;
  }

  if(agb->do_frozen){
#pragma omp barrier
#pragma omp single
    {
      for(i=0;i<nheavyat;i++){
	iat = iheavyat[i];
	if(isvfrozen[iat]){
	  decav_h[iat][0] = decav_const[iat][0];
	  decav_h[iat][1] = decav_const[iat][1];
	  decav_h[iat][2] = decav_const[iat][2];
	}
      }
    }
#pragma omp barrier
  }

#ifdef MMTIMER
  mmtimer_stop( timer3 );
#endif

 ERROR:
  nop;

} /* #pragma omp parallel */

 if(error){
   return AGBNP_ERR;
 } 

/* return derivatives */
 for(iat=0;iat<natoms;iat++){
    for(ki=0;ki<3;ki++){
      decav[iat][ki] = agb->cavw->decav_h[iat][ki];
    }
  }

#ifdef MMTIMER
 mmtimer_stop( cavity_timer );
 mmtimer_report( timer1, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
	"Cavity NB timer", user_time, system_time, real_time);
 mmtimer_report( timer2, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
        "Cavity Surf timer", user_time, system_time, real_time);
 mmtimer_report( timer3, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
        "Cavity Dervs. timer", user_time, system_time, real_time);
 mmtimer_report( cavity_timer, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
        "Cavity tot timer", user_time, system_time, real_time);
#endif

 return AGBNP_OK;
}

int plop_agb_cavity_energy(AGBNPdata *agb, int natom_asu, float_i *x, float_i *y, float_i *z, 
                           float_i *mol_volume, float_i *surf_area,
                           float_i *ecav, float_i *ecorr){
  
  int i, ki; /* counter */
  int iat; /* atomic counter */
  int res;
  int natoms = agb->natoms;
  float_a *igamma = agb->igamma; /* ideal gamma */
  float_a *gamma = agb->cavw->gamma; /* ideal gamma + corrction gamma */
  float_a *sgamma = agb->sgamma; /* correction gamma */
  float_a *gammap = agb->cavw->gammap; /* gamma corrected by derivative of 
					 surface area switching function 
					 (for derivative calculation) */
  float_a *volumep = agb->cavw->volumep; /* self-volumes */
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  float_a a, f, fp;
  float_a (*decav_h)[3] = agb->cavw->decav_h;
  float_a (*decav_const)[3] = agb->cavw->decav_const;
  CAVworkdata *cavw = agb->cavw;
  CAVworkdata *cavw_h;
#ifdef _OPENMP
  int iproc;
#endif
  int error = 0, nop = 0;
  int *isvfrozen = cavw->isvfrozen;
  float_a *volumep_constant = cavw->volumep_constant;
  float_a *surfarea_constant = cavw->surfarea_constant;

#ifdef MMTIMER
  static MMtimer timer1 = -1, timer2 = -1, timer3 = -1;
  double user_time, system_time, real_time;
#endif

#ifdef MMTIMER
  if(cavity_timer < 0) mmtimer_new( &cavity_timer );
  if(timer1<0) mmtimer_new( &timer1);
  if(timer2<0) mmtimer_new( &timer2);
  if(timer3<0) mmtimer_new( &timer3);
  mmtimer_start( cavity_timer );
#endif

  /* copy coordinates to internal buffers */
  agbnp_cpyd(agb->x,x,agb->natoms);
  agbnp_cpyd(agb->y,y,agb->natoms);
  agbnp_cpyd(agb->z,z,agb->natoms);

#ifdef _OPENMP
#pragma omp parallel private(cavw_h,iproc,error)
#endif
{

#ifdef _OPENMP
  /* get local work structure */
  iproc = omp_get_thread_num();
  cavw_h = agb->cavw_p[iproc];
#else
  cavw_h = cavw;
#endif /* _OPENMP */

#ifdef MMTIMER
  mmtimer_start( timer1 );
#endif

  error = 0;

  /* constructs neighbor list */
  if(agb->dopbc || agb->doalt){
    res = agbnp_neighbor_list_cavity_pbc(agb, cavw_h, agb->x, agb->y, agb->z);
  }else{
    res = agbnp_neighbor_list_cavity(agb, cavw_h, agb->x, agb->y, agb->z);
  }
  if(res != AGBNP_OK){
    fprintf(stderr,"plop_agb_cavity_energy(): error in agbnp_neighbor_list_cavity()\n");
    error = 1;
    goto ERROR;
  }

#ifdef MMTIMER
  mmtimer_stop( timer1 );
#endif

#ifdef MMTIMER
  mmtimer_start( timer2 );
#endif

  /* calculates self volumes and surface areas*/
  if(agb->dopbc || agb->doalt){
    res = agbnp_surface_areas_pbc(agb, cavw_h, agb->x, agb->y, agb->z);
  }else{
    res = agbnp_surface_areas(agb, cavw_h, agb->x, agb->y, agb->z);
  }
  if(res != AGBNP_OK){
    fprintf(stderr,"plop_agb_cavity_energy(): error in agbnp_surface_areas()\n");
    error = 1;
    goto ERROR;
  }

#ifdef MMTIMER
  mmtimer_stop( timer2 );
#endif


#pragma omp barrier
#pragma omp single
  {
    /* volume of molecule */
    *mol_volume = 0.0;
    if(!agb->do_frozen){
      for(i = 0; i < nheavyat; i++){
	iat = iheavyat[i];
        if(iat >= natom_asu) continue;
	*mol_volume += volumep[iat];
      }
    }else{
      for(i = 0; i < nheavyat; i++){
	iat = iheavyat[i];
        if(iat >= natom_asu) continue;

	if(isvfrozen[iat]){
	  *mol_volume += volumep_constant[iat];
	}else{
	  *mol_volume += volumep[iat];
	}
      }
    }

    /* filters surface areas to avoid negative surface areas */
    memset(surf_area,0,natoms*sizeof(float_i));

    if(!agb->do_frozen){
      for(i=0; i<nheavyat; i++){
	iat = iheavyat[i];
	/* use occupancy-unweighted surface area for switching function. This
	   way in the limit of identical alternate conformations surface areas
	   are the same as for one conformation */
	a = cavw->surf_area[iat]/agb->occupancy[iat];
	f = swf_area(a, &fp);
	surf_area[iat] = cavw->surf_area[iat]*f;
	gammap[iat] = gamma[iat]*(f+a*fp);
      }
    }else{
      for(i=0;i<nheavyat;i++){
	iat = iheavyat[i];
	if(isvfrozen[iat]){
	  cavw->surf_area[iat] = surfarea_constant[iat];
	}
	a = cavw->surf_area[iat]/agb->occupancy[iat];
	f = swf_area(a, &fp);
	surf_area[iat] = cavw->surf_area[iat]*f;
	gammap[iat] = gamma[iat]*(f+a*fp);
      }
    }
    /* calculates cavity energy */
    *ecav = 0.0;
    *ecorr = 0.0;
    for(i=0;i<nheavyat;i++){
      iat = iheavyat[i];
      if(iat >= natom_asu) continue;
      *ecav += igamma[iat]*surf_area[iat];
      *ecorr += sgamma[iat]*surf_area[iat];
    }
  }
#pragma omp barrier

 ERROR:
  nop;

} /* #pragma omp parallel */

 if(error){
   return AGBNP_ERR;
 } 

#ifdef MMTIMER
 mmtimer_stop( cavity_timer );
 mmtimer_report( timer1, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
	"Cavity NB timer", user_time, system_time, real_time);
 mmtimer_report( timer2, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
        "Cavity Surf timer", user_time, system_time, real_time);
 mmtimer_report( cavity_timer, &user_time, &system_time,  &real_time );
 printf("%40s: %lf %lf %lf\n",
        "Cavity tot timer", user_time, system_time, real_time);
#endif

 return AGBNP_OK;
}


/* return born radii and scaled radii only (no energies, no derivatives) */
int agbnp_bornr(int tag, float_i *x, float_i *y, float_i *z,
		     float_i *sp, float_i *br){
  AGBNPdata *agb;
  AGBworkdata *agbw;    /* shared work space */
  AGBworkdata *agbw_h;  /* work space for this thread */
  int i, iat;
#ifdef _OPENMP
  int iproc;
#endif
  int error = 0, nop =0;
  int res;

  if(!agbnp_initialized){
    fprintf(stderr,"agbnp_bornr(): agbnp library is not initialized.\n");
    return AGBNP_ERR;
  }
  if(!agbnp_tag_ok(tag)){
    fprintf(stderr,"agbnp_bornr(): invalid tag %d.\n",tag);
    return AGBNP_ERR;
  }

  /* pointer to agb data structure */
  agb = &(agbdata_list[tag]);
  agbw = agb->agbw;

  /* copy coordinates to internal buffers */
  agbnp_cpyd(agb->x,x,agb->natoms);
  agbnp_cpyd(agb->y,y,agb->natoms);
  agbnp_cpyd(agb->z,z,agb->natoms);

#ifdef _OPENMP
#pragma omp parallel private(agbw_h, iproc, error)
#endif
  {
    
#ifdef _OPENMP
    /* get local work structure */
    iproc = omp_get_thread_num();
    agbw_h = agb->agbw_p[iproc];
#else
    agbw_h = agb->agbw;
#endif /* _OPENMP */

    error = 0;

    /* constructs neighbor lists */
    if(agb->dopbc || agb->doalt){
      res = agbnp_neighbor_lists_pbc(agb, agbw_h, agb->x, agb->y, agb->z);
    }else{
      res = agbnp_neighbor_lists(agb, agbw_h, agb->x, agb->y, agb->z);
    }
    if(res != AGBNP_OK){
      fprintf(stderr,"agbnp_bornr(): error in agbnp_neighbor_lists()\n");
      error = 1;
      goto ERROR;
    }

    /* computes self volumes  */
    if(agb->dopbc || agb->doalt){
      res = agbnp_self_volumes_pbc(agb, agbw_h, agb->x, agb->y, agb->z);
    }else{
      res = agbnp_self_volumes(agb, agbw_h, agb->x, agb->y, agb->z);
    }
    if(res != AGBNP_OK){
      fprintf(stderr, "agbnp_bornr(): error in agbnp_self_volumes()\n");
      error = 1;
      goto ERROR;
    }

    /* calculate volume scaling factors */
    if(agbnp_scaling_factors(agb, agbw_h) != AGBNP_OK){
      fprintf(stderr, "agbnp_bornr(): error in agbnp_scaling_factors()\n");
      error = 1;
      goto ERROR;
    }

    /* calculates inverse Born radii */
    if(agb->dopbc || agb->doalt){
      res = agbnp_inverse_born_radii_pbc(agb, agbw_h, agb->x, agb->y, agb->z, 0);
    }else{
      res = agbnp_inverse_born_radii(agb, agbw_h, agb->x, agb->y, agb->z, 0);
    }
    if(res != AGBNP_OK){
      fprintf(stderr, "agbnp_bornr(): error in agbnp_inverse_born_radii()\n");
      error = 1;
      goto ERROR;
    }
    
    /* calculates Born radii and related quantities (brw, etc.) */
    if(agbnp_born_radii(agb, agbw_h) != AGBNP_OK){
      fprintf(stderr, "agb_bornr(): error in agbnp_born_radii()\n");
      error = 1;
      goto ERROR;
    }

  ERROR:
    nop;

  }/* #pragma omp parallel */

  if(error){
    return AGBNP_ERR;
  }

  /* returns scaled volume factors */
  for(iat=0;iat<agb->natoms;iat++){
    sp[iat] = 1.0;
  }
  for(i=0;i<agb->nheavyat;i++){
    iat = agb->iheavyat[i];
    sp[iat] = agbw->sp[iat];
  }

  /* return born radii */
  for(iat=0;iat<agb->natoms;iat++){
    br[iat] = agbw->br[iat];
  }

  return AGBNP_OK;
}

#ifdef AGBNP_HB
int agbnp_create_wsatoms(AGBNPdata *agb){
  int iat, iws, nws;
  int jat, kat, lat;
  int i, incr;
  static const float_a wsrad = 1.4;
  static const int min_incr = 20;

  /* check that connection table exists */
  if(!agb->conntbl){
    fprintf(stderr, "agbnp_create_watoms(): unable to retrieve connection table.\n");
    return AGBNP_ERR;
  }

  /* allocate master list of ws atoms */
  nws = 1*agb->nheavyat; /* initial number of ws atoms */
  agb->wsat = (WSat *)calloc(nws,sizeof(WSat));
  if(!agb->wsat){
    fprintf(stderr, "agbnp_create_watoms(): unable to allocate list of ws atoms (%d WSat)\n",nws);
    return AGBNP_ERR;
  }
  agb->wsat_size = nws;

  /* allocate list of ws atoms for each real atom */
  agb->wsat_list = (NeighList *)calloc(1,sizeof(NeighList));
  if(!agb->wsat_list){
    fprintf(stderr, "agbnp_create_watoms(): unable to allocate msat_list neighbor list structure.\n");
    return AGBNP_ERR;
  }
  nblist_reset_neighbor_list(agb->wsat_list);
  if(nblist_reallocate_neighbor_list(agb->wsat_list,agb->natoms,
				     agb->wsat_size) != NBLIST_OK){
      fprintf(stderr, "agbnp_create_watoms(): error in nblist_reallocate_neighbor_list().\n");
      return AGBNP_ERR;   
  }
  
  /* create ws atoms for each HB active atom */
  agb->nwsat = 0;
  incr = agb->natoms/10;
  if(incr<min_incr) incr = min_incr;
  for(iat = 0; iat < agb->natoms; iat++){
    agb->wsat_list->neighl[iat] = &(agb->wsat_list->neighl1[agb->nwsat]);
    agb->wsat_list->nne[iat] = 0;
    if(agb->hbtype[iat] == AGBNP_HB_INACTIVE) continue;
    if(agb->wsat_size - agb->nwsat < min_incr){
      nws = agb->wsat_size + incr;
      agb->wsat = (WSat *)realloc(agb->wsat,nws*sizeof(WSat));
      if(!agb->wsat){
	fprintf(stderr, "agbnp_create_watoms(): unable to re-allocate list of ws atoms (%d WSat)\n",nws);
	return AGBNP_ERR;
      }
      if(nblist_reallocate_neighbor_list(agb->wsat_list,agb->natoms,
					 nws) != NBLIST_OK){
	fprintf(stderr, "agbnp_create_watoms(): error in nblist_reallocate_neighbor_list().\n");
	return AGBNP_ERR;   
      }
      agb->wsat_size = nws;
    }
    if(agbnp_create_ws_ofatom(agb,iat,agb->nwsat,&nws) != AGBNP_OK){
      fprintf(stderr, "agbnp_create_watoms(): error in agbnp_create_ws_ofatom()\n");
      return AGBNP_ERR;
    }
    agb->nwsat += nws;
  }
  
  /* computes free volume of water sites */
  if(agbnp_ws_free_volumes(agb) != AGBNP_OK){
    fprintf(stderr, "agbnp_create_watoms(): error in agbnp_ws_free_volumes(agb)\n");
    return AGBNP_ERR;
  }


  fprintf(stderr,"%d\n\n",agb->nwsat);
  for(iws = 0; iws < agb->nwsat ; iws++){

    fprintf(stderr,"X %lf %lf %lf   %lf\n",
            agb->wsat[iws].pos[0],agb->wsat[iws].pos[1],agb->wsat[iws].pos[2],
	    agb->wsat[iws].sp);

    /*    fprintf(stderr,"ms %d: %d %d %d %f %f %f\n",ims, agb->msat[ims].parent1,
	    agb->msat[ims].parent2,agb->msat[ims].parent3,
	    agb->msat[ims].pos[0],agb->msat[ims].pos[1],agb->msat[ims].pos[2]);
    */

  }

  {
    float_a ce;
    if(agbnp_ws_correction(agb, &ce) != AGBNP_OK){
      fprintf(stderr, "agbnp_create_watoms(): error in agbnp_ws_correction()\n");
      return AGBNP_ERR;
    }
    fprintf(stderr,"HB correction energy: %lf\n",ce);
  }

  return AGBNP_OK;
}


/* create water sites pseudo atoms for atom iat, stores ws atoms
   starting at location iws and returns the number of added ws atoms
   in nws */
int agbnp_create_ws_ofatom(AGBNPdata *agb, int iat, int iws, int *nws){
  *nws = 0;
  switch (agb->hbtype[iat]){
  case AGBNP_HB_INACTIVE:
    return AGBNP_OK;
  case AGBNP_HB_POLARH:
    if(agbnp_create_ws_atoms_ph(agb, iat, iws, nws) != AGBNP_OK){
      fprintf(stderr, "agbnp_create_ws_atoms(): error in agbnp_create_ws_atoms_ph()\n");
      return AGBNP_ERR;
    }
    return AGBNP_OK;
  case AGBNP_HB_SP2ACCPT:
    if(agbnp_create_ws_atoms_sp2(agb, iat, iws, nws) != AGBNP_OK){
      fprintf(stderr, "agbnp_create_ws_atoms(): error in agbnp_create_ws_atoms_sp2()\n");
      return AGBNP_ERR;
    }
    return AGBNP_OK;
  case AGBNP_HB_SP3ACCPT:
    if(agb->conntbl->nne[iat] == 2){
      if(agbnp_create_ws_atoms_sp3_o(agb, iat, iws, nws) != AGBNP_OK){
	fprintf(stderr, "agbnp_create_ws_atoms(): error in agbnp_create_ws_atoms_sp3_o()\n");
	return AGBNP_ERR;
      }
    }
    return AGBNP_OK;
  default:
    return AGBNP_OK;
  }
}

int agbnp_create_ws_atoms_ph(AGBNPdata *agb, int iat, int iws, int *nws){
  int jat = -1; /* heavy atom the hydrogen is attached to */
  int i;
  float_a xw, yw, zw;
  float_a *x = agb->x;
  float_a *y = agb->y;
  float_a *z = agb->z;
  NeighList *conntbl = agb->conntbl;
  WSat *wsat;
  NeighList *wsat_list = agb->wsat_list;
  

  if(conntbl->nne[iat] <= 0){
    fprintf(stderr, "agbnp_create_ws_atoms_ph(): atom %d is isolated, I can't do isolated atoms yet.\n", iat);
    return AGBNP_ERR;
  }
  /* find parent */
  for(i=0;i<agb->conntbl->nne[iat];i++){
    if(agb->agbw->isheavy[conntbl->neighl[iat][i]]){
      jat = conntbl->neighl[iat][i];
      break;
    }
  }
  if(jat<0){
    fprintf(stderr, "agbnp_create_ws_atoms_ph(): unable to find parent heavy atom of atom %d.\n", iat);
    return AGBNP_ERR;
  }
  agbnp_place_wat_hydrogen(x[jat],y[jat],z[jat],
			   x[iat],y[iat],z[iat],
			   (float_a)AGBNP_HB_LENGTH,
			   &xw, &yw, &zw);
  /* stores water sites */
  wsat = &(agb->wsat[iws]);
  wsat->pos[0] = xw;
  wsat->pos[1] = yw;
  wsat->pos[2] = zw;
  wsat->r = AGBNP_HB_RADIUS;
  wsat->volume = (4./3.)*pi*pow(wsat->r,3);
  wsat->nparents = 2;
  wsat->parent[0] = iat;
  wsat->parent[1] = jat;
  wsat_list->neighl1[iws] = iws; 
  wsat_list->nne[iat] = 1;

  *nws = 1;
  return AGBNP_OK;
}


/* places a water site at distance d from donor along hydrogen-donor bond 
   Input:
    xd: position of heavy atom
    xh: position of hydrogen atom
    d: distance between donor and water site
   Output:
    xw: position of water site
*/
void agbnp_place_wat_hydrogen(float_a xd, float_a yd, float_a zd,
			      float_a xh, float_a yh, float_a zh, 
			      float_a d, 
			      float_a *xw, float_a *yw, float_a *zw){
  float_a dx, dy, dz, w;

  dx = xh - xd;
  dy = yh - yd;
  dz = zh - zd;
  w = d/sqrt(dx*dx+dy*dy+dz*dz);
  *xw = xd + w*dx;
  *yw = yd + w*dy;
  *zw = zd + w*dz;
}

int agbnp_create_ws_atoms_sp2(AGBNPdata *agb, int iat, int iws, int *nws){
  int i;
  int ir, ir1, ir2;
  float_a xw1, yw1, zw1;
  float_a xw2, yw2, zw2;
  float_a *x = agb->x;
  float_a *y = agb->y;
  float_a *z = agb->z;
  NeighList *conntbl = agb->conntbl;
  WSat *wsat;
  NeighList *wsat_list = agb->wsat_list;
  
  /* assumed topology:

       R1
         \
          R-A 
         /
       R2
  */

  if(conntbl->nne[iat] != 1){
    fprintf(stderr, "agbnp_create_ws_atoms_sp2(): acceptor atom (%d) should have one connecting atoms. Found %d.\n", iat,conntbl->nne[iat]);
    return AGBNP_ERR;
  }
  /* R parent */
  ir = conntbl->neighl[iat][0];
  /* R1 and R2 parents */
  if(conntbl->nne[ir] != 3){
    fprintf(stderr, "agbnp_create_ws_atoms_sp2(): atom bound to acceptor atom (%d) should have 3 connecting atoms. Found %d.\n", ir,conntbl->nne[ir]);
    return AGBNP_ERR;
  }
  ir1 = iat;
  i = 0;
  while(ir1 == iat){
    if(i >= conntbl->nne[ir]){
      fprintf(stderr, "agbnp_create_ws_atoms_sp2(): error parsing connection table of R atom (atom %d).\n", ir);
    return AGBNP_ERR;
    }
    ir1 = conntbl->neighl[ir][i++];
  }
  ir2 = iat;
  i = 0;
  while(ir2 == iat || ir2 == ir1){
    if(i >= conntbl->nne[ir]){
      fprintf(stderr, "agbnp_create_ws_atoms_sp2(): error parsing connection table of R atom (atom %d).\n", ir);
    return AGBNP_ERR;
    }
    ir2 = conntbl->neighl[ir][i++];
  } 
  agbnp_place_wat_lp_sp2(x[iat],y[iat],z[iat],
			 x[ir],y[ir],z[ir],
			 x[ir1],y[ir1],z[ir1],
			 x[ir2],y[ir2],z[ir2],
			 (float_a)AGBNP_HB_LENGTH,
			 &xw1, &yw1, &zw1,
			 &xw2, &yw2, &zw2);
  /* stores water sites */
  i = iws;
  wsat = &(agb->wsat[i]);
  wsat->pos[0] = xw1;
  wsat->pos[1] = yw1;
  wsat->pos[2] = zw1;
  wsat->r = AGBNP_HB_RADIUS;
  wsat->volume = (4./3.)*pi*pow(wsat->r,3);
  wsat->nparents = 3;
  wsat->parent[0] = iat;
  wsat->parent[1] = ir;
  wsat->parent[2] = ir1;
  wsat_list->neighl1[i] = i;

  i = iws + 1;
  wsat = &(agb->wsat[i]);
  wsat->pos[0] = xw2;
  wsat->pos[1] = yw2;
  wsat->pos[2] = zw2;
  wsat->r = AGBNP_HB_RADIUS;
  wsat->volume = (4./3.)*pi*pow(wsat->r,3);
  wsat->nparents = 3;
  wsat->parent[0] = iat;
  wsat->parent[1] = ir;
  wsat->parent[2] = ir2;
  wsat_list->neighl1[i] = i;

  wsat_list->nne[iat] = 2;

  *nws = 2;
  return AGBNP_OK;
}


/* places water sites for sp2 acceptors such as the carboxy oxygen:
   Input:
    xa: position of acceptor
    xr: position of atom bound to acceptor
    xr1, xr2: positions of atoms flanking atom r
    d: distance between acceptor and water site
   Output:
    xw1,xw2: positions of water sites 
*/
void agbnp_place_wat_lp_sp2(float_a xa, float_a ya, float_a za, 
			    float_a xr, float_a yr, float_a zr, 
			    float_a xr1, float_a yr1, float_a zr1, 
			    float_a xr2, float_a yr2, float_a zr2,
			    float_a d,
			    float_a *xw1, float_a *yw1, float_a *zw1,
			    float_a *xw2, float_a *yw2, float_a *zw2){
  float_a dx, dy, dz, w;

  dx = xr1 - xr;
  dy = yr1 - yr;
  dz = zr1 - zr;
  w = -d/sqrt(dx*dx+dy*dy+dz*dz);
  *xw1 = xa + w*dx;
  *yw1 = ya + w*dy;
  *zw1 = za + w*dz;

  dx = xr2 - xr;
  dy = yr2 - yr;
  dz = zr2 - zr;
  w = -d/sqrt(dx*dx+dy*dy+dz*dz);
  *xw2 = xa + w*dx;
  *yw2 = ya + w*dy;
  *zw2 = za + w*dz;
}

int agbnp_create_ws_atoms_sp3_o(AGBNPdata *agb, int iat, int iws, int *nws){
  int i;
  int ir1, ir2;
  float_a xw1, yw1, zw1;
  float_a xw2, yw2, zw2;
  float_a *x = agb->x;
  float_a *y = agb->y;
  float_a *z = agb->z;
  NeighList *conntbl = agb->conntbl;
  WSat *wsat;
  NeighList *wsat_list = agb->wsat_list;
  
  /* assumed topology:

       R1
         \
          A 
         /
       R2
  */

  if(conntbl->nne[iat] != 2){
    fprintf(stderr, "agbnp_create_ws_atoms_sp3_o(): acceptor atom (%d) should have two connecting atoms. Found %d.\n", iat,conntbl->nne[iat]);
    return AGBNP_ERR;
  }
  /* R1 and R2 parents */
  ir1 = conntbl->neighl[iat][0];
  ir2 = conntbl->neighl[iat][1];
  agbnp_place_wat_lp_sp3_o(x[iat],y[iat],z[iat],
			   x[ir1],y[ir1],z[ir1],
			   x[ir2],y[ir2],z[ir2],
			   (float_a)AGBNP_HB_LENGTH,
			   &xw1, &yw1, &zw1,
			   &xw2, &yw2, &zw2);
  /* stores water sites */
  i = iws;
  wsat = &(agb->wsat[i]);
  wsat->pos[0] = xw1;
  wsat->pos[1] = yw1;
  wsat->pos[2] = zw1;
  wsat->r = AGBNP_HB_RADIUS;
  wsat->volume = (4./3.)*pi*pow(wsat->r,3);
  wsat->nparents = 3;
  wsat->parent[0] = iat;
  wsat->parent[1] = ir1;
  wsat->parent[2] = ir2;
  wsat_list->neighl1[i] = i;

  i = iws + 1;
  wsat = &(agb->wsat[i]);
  wsat->pos[0] = xw2;
  wsat->pos[1] = yw2;
  wsat->pos[2] = zw2;
  wsat->r = AGBNP_HB_RADIUS;
  wsat->volume = (4./3.)*pi*pow(wsat->r,3);
  wsat->nparents = 3;
  wsat->parent[0] = iat;
  wsat->parent[1] = ir1;
  wsat->parent[2] = ir2;
  wsat_list->neighl1[i] = i;

  wsat_list->nne[iat] = 2;

  *nws = 2;
  return AGBNP_OK;
}

void agbnp_place_wat_lp_sp3_o(float_a xa, float_a ya, float_a za, 
			      float_a xr1, float_a yr1, float_a zr1, 
			      float_a xr2, float_a yr2, float_a zr2,
			      float_a d,
			      float_a *xw1, float_a *yw1, float_a *zw1,
			      float_a *xw2, float_a *yw2, float_a *zw2){
  float_a costh = -0.577350269; /* -1/sqrt(3) */
  float_a sinth =  0.816496581; /* sqrt(2/3) */
  float_a dx, dy, dz, d1[3], d2[3], dpin[3], dout[3], w;
  float_a dist1, u[3];

  /* unit vectors of bonds */
  dx = xr1 - xa;
  dy = yr1 - ya;
  dz = zr1 - za;
  dist1 =  1./sqrt(dx*dx+dy*dy+dz*dz);
  d1[0] = dist1*dx;
  d1[1] = dist1*dy;
  d1[2] = dist1*dz;

  dx = xr2 - xa;
  dy = yr2 - ya;
  dz = zr2 - za;
  dist1 =  1./sqrt(dx*dx+dy*dy+dz*dz);
  d2[0] = dist1*dx;
  d2[1] = dist1*dy;
  d2[2] = dist1*dz;

  /* in plane direction (sum of bond unit vectors) */
  dx = d1[0] + d2[0];
  dy = d1[1] + d2[1];
  dz = d1[2] + d2[2];
  w =  costh/sqrt(dx*dx+dy*dy+dz*dz);
  dpin[0] = w*dx;
  dpin[1] = w*dy;
  dpin[2] = w*dz;

  /* out of plane projection (cross product between bond vectors) */
  dx = d2[1]*d1[2] - d2[2]*d1[1];
  dy = d2[2]*d1[0] - d2[0]*d1[2];
  dz = d2[0]*d1[1] - d2[1]*d1[0];
  w =  sinth/sqrt(dx*dx+dy*dy+dz*dz);
  dout[0] = w*dx;
  dout[1] = w*dy;
  dout[2] = w*dz;

  /* first LP position */
  u[0] = dpin[0] + dout[0];
  u[1] = dpin[1] + dout[1];
  u[2] = dpin[2] + dout[2];
  *xw1 = xa + d*u[0];
  *yw1 = ya + d*u[1];
  *zw1 = za + d*u[2];
  /* second LP position */
  u[0] = dpin[0] - dout[0];
  u[1] = dpin[1] - dout[1];
  u[2] = dpin[2] - dout[2];
  *xw2 = xa + d*u[0];
  *yw2 = ya + d*u[1];
  *zw2 = za + d*u[2];
}


/* computes/updates the free volumes of water site atoms */
int agbnp_ws_free_volumes(AGBNPdata *agb){
  float_a *x = agb->x;
  float_a *y = agb->y;
  float_a *z = agb->z;
  AGBworkdata *agbw = agb->agbw;
  int iws;
  
  int iv2cache = 0;
  /* coordinate buffer for Gaussian overlap calculation */
  float_a gx[AGBNP_MAX_OVERLAP_LEVEL][3];
  /* radius buffer for  Gaussian overlap calculation */
  float_a gr[AGBNP_MAX_OVERLAP_LEVEL];
  /* Gaussian parameters buffers for overlap calculation */
  float_a ga[AGBNP_MAX_OVERLAP_LEVEL];
  float_a gp[AGBNP_MAX_OVERLAP_LEVEL];
  /* holds the atom indexes being overlapped */
  int gnlist[AGBNP_MAX_OVERLAP_LEVEL];
  int gatlist[AGBNP_MAX_OVERLAP_LEVEL];
  float_a gvol; /* gaussian overlap volume */
  int order; /* order of the overlap */
  /* holds overlap Gaussian parameters at each level */
  GParm gparams[AGBNP_MAX_OVERLAP_LEVEL];
  float_a an, pn, cn[3];
  /* coefficients for scaled volume */
  float_a volpcoeff[AGBNP_MAX_OVERLAP_LEVEL];
  float_a sign;

  int i, iat, jat, kat, ii,j;
  float_a u,v,altw=1.0;

  int natoms = agb->natoms;
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  float_a *occupancy = agb->occupancy;
  NeighList *near_nl = agbw->near_nl;
  float_a *r = agb->r;
  float_a *galpha = agbw->galpha;
  float_a *gprefac = agbw->gprefac;
  float_a *volumep = agbw->volumep;
  float *v2cache = agbw->v2cache;
  int do_frozen = agb->do_frozen;
  int *isvfrozen = agb->agbw->isvfrozen;
  int *isvvfrozen = agb->agbw->isvvfrozen;

  const float_a kf_ws = KFC;
  const float_a pf_ws = PFC;

  /*  const float_a kf_ws = 1.0;
  const float_a pf_ws = 0.7522;
  */

  WSat *wsat;

  /* set scaled volume coefficients */
  sign = 1.0;
  for(order=1;order<=AGBNP_MAX_OVERLAP_LEVEL;order++){
    volpcoeff[order-1] = sign;
    sign *= -1.0;
  }

  /* computes self volumes */
  for(iws=0;iws<agb->nwsat;iws++){
    
    wsat = &(agb->wsat[iws]);
    wsat->free_volume = wsat->volume;
    
    gatlist[0] = iws;
    gx[0][0] = wsat->pos[0];
    gx[0][1] = wsat->pos[1];
    gx[0][2] = wsat->pos[2];
    ga[0] = kf_ws/(wsat->r*wsat->r);
    gp[0] = pf_ws;
    gr[0] = wsat->r;
    gparams[0].a = ga[0];
    gparams[0].p = gp[0];
    gparams[0].c[0] = gx[0][0];
    gparams[0].c[1] = gx[0][1];
    gparams[0].c[2] = gx[0][2];

    order = 2;
    gnlist[order-1] = 0;

    while(order>1){

      j = gnlist[order-1];
      jat = agb->iheavyat[j];
      gatlist[order-1] = jat;

      gx[order-1][0] = x[jat];
      gx[order-1][1] = y[jat];
      gx[order-1][2] = z[jat];
      ga[order-1] = galpha[jat];
      gp[order-1] = gprefac[jat];
      gr[order-1] = r[jat];


      gvol = ogauss_incremental(
		  order, gx, ga, gp, gr,
		  gparams[order-2].a, gparams[order-2].p, gparams[order-2].c,
		  &an, &pn, cn,
		  AGBNP_MIN_VOLA,AGBNP_MIN_VOLB,
		  NULL, NULL, NULL);

      gparams[order-1].a = an;
      gparams[order-1].p = pn;
      gparams[order-1].c[0] = cn[0];
      gparams[order-1].c[1] = cn[1];
      gparams[order-1].c[2] = cn[2];

      if(altw*gvol > FLT_MIN && order < AGBNP_MAX_OVERLAP_LEVEL){

	wsat->free_volume += volpcoeff[order-1]*gvol;
	
	/* increment overlap level */
	order += 1;
	/* initialize neighbor index at this level
	   (neighbor index at previous level + 1)
	*/
	gnlist[order-1] = gnlist[order-2]+1;
      }else{
	/* new neighbor at this level */
	gnlist[order-1] += 1;
      }

      /* 
	 if neighbor index is out of bounds there are no more neighbors at
	 this level. Decrement order and go to next neighbor at that order.
	 Do this recursively until an order with available neighbors is found.
      */
      while(order>1 && gnlist[order-1]>=agb->nheavyat){
	order -= 1;
	gnlist[order-1] += 1;
      }

    }
    wsat->sp = wsat->free_volume/wsat->volume;
  }	

  return AGBNP_OK;
}


/* computes the correction energy based on the water sites */
int agbnp_ws_correction(AGBNPdata *agb, float_a *ce){
  float_a khb = 0.0; /* -3.0; */
  float_a xa = 0.75;
  float_a xb = 1.0;

  float_a kws = 1.0;
  float_a wa = 0.5;
  float_a wb = 1.0;
  float_a dp, d, dx, dy, dz;
  WSat *wsat1, *wsat2;

  float_a s;
  int iws, jws;

  *ce = 0.0;
  for(iws=0;iws<agb->nwsat;iws++){
    s = agbnp_pol_switchfunc(agb->wsat[iws].sp, xa, xb, NULL, NULL);
    *ce += khb*s;
  }

  for(iws = 0; iws < agb->nwsat ; iws++){
    wsat1 = &(agb->wsat[iws]);
    for(jws = iws + 1 ; jws <  agb->nwsat ; jws++){
      wsat2 =  &(agb->wsat[jws]);
      dx = wsat2->pos[0] - wsat1->pos[0];
      dy = wsat2->pos[1] - wsat1->pos[1];
      dz = wsat2->pos[2] - wsat1->pos[2];
 
      d = sqrt(dx*dx+dy*dy+dz*dz);
      dp = d/(wsat1->r+wsat2->r);
      s = 1. - agbnp_pol_switchfunc(dp, wa, wb, NULL, NULL);
      *ce += kws*wsat1->sp*wsat2->sp*s;

    }
  }

  return AGBNP_OK;
}

#ifdef NOTNOW

int agbnp_find_wparents(AGBNPdata *agb, int iat, int *jat, int *kat, int *lat){
  int i;
  
  *jat = -1;
  *kat = -1;

  if(agb->conntbl->nne[iat] <= 0){
    fprintf(stderr, "agbnp_find_wparents(): atom %d is isolated, I can't do isolated atoms yet.\n", iat+1);
    return AGBNP_ERR;
  }

  /* first choice for jat: first heavy atom in connection table */
  for(i=0;i<agb->conntbl->nne[iat];i++){
    if(agb->agbw->isheavy[agb->conntbl->neighl[iat][i]]){
      *jat = agb->conntbl->neighl[iat][i];
      break;
    }
  }

  if(*jat >= 0){
    /* first choice for kat: second heavy atom in connection table */
    for(i=0;i<agb->conntbl->nne[iat];i++){
      if(agb->agbw->isheavy[agb->conntbl->neighl[iat][i]]){
	if(agb->conntbl->neighl[iat][i] != *jat){
	  *kat = agb->conntbl->neighl[iat][i];
	  break;
	}
      }
    }
    if(*kat < 0){
      /* second choice for kat: first heavy atom in connection table of jat */
      for(i=0;i<agb->conntbl->nne[*jat];i++){
	if(agb->agbw->isheavy[agb->conntbl->neighl[*jat][i]] &&
	   agb->conntbl->neighl[*jat][i] != iat){
	  *kat = agb->conntbl->neighl[*jat][i];
	  break;
	}
      }
    }
    if(*kat < 0){
      /* third choice for kat: any other atom in connection table of iat */
      for(i=0;i<agb->conntbl->nne[iat];i++){
	if(agb->conntbl->neighl[iat][i] != *jat){
	  *kat = agb->conntbl->neighl[iat][i];
	  break;
	}
      }
    }
    if(*kat < 0){
      /* give up */
      fprintf(stderr, "agbnp_find_wparents(): 1 - can't find any kat for atom %d (jat = %d).\n", iat+1, *jat+1);
      return AGBNP_ERR;
    }
    return AGBNP_OK;
  }

  if(*jat < 0){
    /* second choice for jat: first atom in connection table
       (must be a hydrogen)
    */
    *jat = agb->conntbl->neighl[iat][0];
  }
  /* if jat is the only atom give up */
  if(agb->conntbl->nne[iat] < 2){
    fprintf(stderr, "agbnp_find_wparents(): 2 - can't find any kat for atom %d (jat = %d).\n", iat+1, *jat+1);
    return AGBNP_ERR;
  }
  /* otherwise kat is set as the second hydrogen */
  *kat = agb->conntbl->neighl[iat][1];

  return AGBNP_OK;
}

 int agbnp_calc_waxes(AGBNPdata *agb, int iat, int jat, int kat,
			   float_a origin[3], 
			   float_a x[3], float_a y[3], float_a z[3]){
  float_a d,dotp;
  float_a dr[3];
  int i;

  /* the origin is the location of iat */
  origin[0] = agb->x[iat];
  origin[1] = agb->y[iat];
  origin[2] = agb->z[iat];
  
  /* z axis: versor from iat to jat */
  z[0] = agb->x[jat] - origin[0];
  z[1] = agb->y[jat] - origin[1];
  z[2] = agb->z[jat] - origin[2];
  d = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
  if(d < FLT_MIN){
    fprintf(stderr, "agbnp_calc_waxes(): atoms %d and %d have identical positions.\n", iat, jat);
    return AGBNP_ERR;
  }
  z[0] /= d; z[1] /= d; z[2] /= d;
  
  /* x axis: vector from kat to iat minus projection along z axis */
  dr[0] = agb->x[kat] - origin[0];
  dr[1] = agb->y[kat] - origin[1];
  dr[2] = agb->z[kat] - origin[2];
  dotp = dr[0]*z[0] + dr[1]*z[1] + dr[2]*z[2];
  dr[0] -= dotp*z[0];
  dr[1] -= dotp*z[1];
  dr[2] -= dotp*z[2];
  d = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
  if(d < FLT_MIN){
    fprintf(stderr, "agbnp_calc_waxes(): atoms %d, %d, and %d are colinear.\n", iat, jat, kat);
    return AGBNP_ERR;
  }
  x[0] = dr[0]/d; x[1] = dr[1]/d; x[2] = dr[2]/d;

  /* y axis: cross product between z and x axes */
  y[0] = z[1]*x[2] - z[2]*x[1];
  y[1] = z[2]*x[0] - z[0]*x[2];
  y[2] = z[0]*x[1] - z[1]*x[0];

  return AGBNP_OK;
}

 int agbnp_init_watoms(AGBNPdata *agb){
  int iat;


  return AGBNP_OK;
}

 int agbnp_position_watoms(AGBNPdata *agb){
  int iat,jat,kat,ims,i;
  float_a origin[3];
  float_a x[3], y[3], z[3];
  static float_a c[6][3] = {
    { 0.0,  0.0,  1.8},
    { 0.0,  0.0, -1.8},
    { 1.8,  0.0,  0.0},
    {-1.8,  0.0,  0.0},
    { 0.0,  1.8,  0.0},
    { 0.0, -1.8,  0.0},
  };

  for(iat = 0; iat < agb->natoms; iat++){
    if(!agb->agbw->isheavy[iat]) continue;
    ims = agb->msat_list->neighl[iat][0]; /* the first ms sphere for iat */
    jat = agb->msat[ims].parent2;
    kat = agb->msat[ims].parent3;
    if(agbnp_calc_waxes(agb,iat,jat,kat,origin,x,y,z) != AGBNP_OK){
      fprintf(stderr, "agbnp_position_watoms(): error in agbnp_calc_waxes()\n");
      return AGBNP_ERR;
    }
    for(i=0;i<agb->msat_list->nne[iat];i++){
      ims = agb->msat_list->neighl[iat][i];
      agb->msat[ims].pos[0] = c[i][0]*x[0]+c[i][1]*y[0]+c[i][2]*z[0] + origin[0];
      agb->msat[ims].pos[1] = c[i][0]*x[1]+c[i][1]*y[1]+c[i][2]*z[1] + origin[1];
      agb->msat[ims].pos[2] = c[i][0]*x[2]+c[i][1]*y[2]+c[i][2]*z[2] + origin[2];
    }
  }

  return AGBNP_OK;
}
#endif /* NOTNOW */
#endif /* AGBNP_HB */

 int agbnp_create_msatoms(AGBNPdata *agb){
  AGBworkdata *agbw = agb->agbw;
  AGBworkdata *agbw_h;
  float_a *x = agb->x, *y = agb->y, *z = agb->z;
  int error = 0;
  int jat, iat, i, j, k, kat;
  float_a dx, dy, dz, d2, d, r, rsum, u, v;
  float_a msx, msy, msz;
  static const float_a msrad = 1.4;
  int nop = 0;
  int overlap;

  agbw_h = agbw;
  
  {

    /* constructs neighbor lists */
    if(agbnp_neighbor_lists(agb, agbw_h, x, y, z)!= AGBNP_OK){
      fprintf(stderr,"agbnp_create_msatoms(): error in agbnp_neighbor_lists()\n");
      error = 1;
      goto ERROR;
    }

    /* create ms atoms for each suitable atom pair */
    agb->nmsat;
    agb->msat_size = 10*agb->natoms; /* first guess for list of msatoms */
    /* allocate master list of ms atoms */
    agb->msat = (MSat *)realloc(agb->msat,agb->msat_size*sizeof(MSat));
    if(!agb->msat){
      fprintf(stderr, "agbnp_create_msatoms(): unable to allocate list of ms atoms (%d MSat)\n",agb->msat_size);
      error = 1;
      goto ERROR;
    }
    /* go through far neighbor list */
    for(i=0;i<agb->nheavyat;i++){
      iat = agb->iheavyat[i];
      for(j=0;j<agbw_h->far_nl->nne[iat];j++){
	jat = agbw_h->far_nl->neighl[iat][j];
	if(!agbw_h->isheavy[jat]) continue;
	dx = x[jat] - x[iat];
	dy = y[jat] - y[iat];
	dz = z[jat] - z[iat];
	d2 = dx*dx + dy*dy + dz*dz;
	d = sqrt(d2);
	rsum = agb->r[iat] + agb->r[jat];
	u = (d - rsum)/2.*msrad;
	if(u < 0.0 || u > 0.75) continue;
	/* pair supports a ms atom */
	v = (agb->r[iat]+0.5*(d-rsum))/d; 
	msx = x[iat] + v*(x[jat] - x[iat]);
	msy = y[iat] + v*(y[jat] - y[iat]);
	msz = z[iat] + v*(z[jat] - z[iat]);
	/* reject if it overlaps with any heavy atom */
	overlap = 0;
	for(k=0;k<agb->nheavyat;k++){
	  kat = agb->iheavyat[k];
	  dx = x[kat] - msx;
	  dy = y[kat] - msy;
	  dz = z[kat] - msz;
	  d2 = dx*dx + dy*dy + dz*dz;
	  r = agb->r[kat];
	  if(d2 < r*r) {
	    overlap = 1;
	    break;
	  }
	}
	if(overlap) continue;
	/* add ms atom */
	if(agb->nmsat >= agb->msat_size){
	  agb->msat_size = 2.*agb->msat_size;
	  agb->msat = (MSat *)realloc(agb->msat,agb->msat_size*sizeof(MSat));
	  if(!agb->msat){
	    fprintf(stderr, "agbnp_create_msatoms(): unable to allocate list of ms atoms (%d MSat)\n",agb->msat_size);
	    error = 1;
	    goto ERROR;
	  }
	}
	agb->msat[agb->nmsat].pos[0] = 0.5*(x[jat] + x[iat]);
	agb->msat[agb->nmsat].pos[1] = 0.5*(y[jat] + y[iat]);
	agb->msat[agb->nmsat].pos[2] = 0.5*(z[jat] + z[iat]);
	agb->msat[agb->nmsat].r = msrad;
	agb->msat[agb->nmsat].parent1 = iat;
	agb->msat[agb->nmsat].parent1 = jat;
	agb->msat[agb->nmsat].s = 1.;
	agb->nmsat += 1;
      }
    }
  
  ERROR:
    nop;
  }

  if(error){
    return AGBNP_ERR;
  }

#ifdef AGBNP_MS_DEBUG
  {
    int ims;
    fprintf(stderr,"%d\n",agb->nmsat);
    for(ims = 0; ims < agb->nmsat ; ims++){

      fprintf(stderr,"X %f %f %f\n",
	      agb->msat[ims].pos[0],agb->msat[ims].pos[1],agb->msat[ims].pos[2]);
    }
  }
#endif

  return AGBNP_OK;
}

/* Solves a linear system using Gaussian elimination
 return 1 if system not solving
 nDim - system dimension
 pfMatr - matrix with coefficients
 pfVect - vector with free members
 pfSolution - vector with system solution
 pfMatr becames trianglular after function call
 pfVect changes after function call

 Credits: Henry Guennadi Levkin www.hlevkin.com/NumAlg/LinearEquations.c
*/
static int LinearEquationsSolving(int nDim, float_a *pfMatr, float_a *pfVect,
				  float_a *pfSolution)
{
  float_a fMaxElem;
  float_a fAcc;

  int i , j, k, m;


  for(k=0; k<(nDim-1); k++) // base row of matrix
  {
    // search of line with max element
    fMaxElem = fabs( pfMatr[k*nDim + k] );
    m = k;
    for(i=k+1; i<nDim; i++)
    {
      if(fMaxElem < fabs(pfMatr[i*nDim + k]) )
      {
        fMaxElem = pfMatr[i*nDim + k];
        m = i;
      }
    }
    
    // permutation of base line (index k) and max element line(index m)
    if(m != k)
    {
      for(i=k; i<nDim; i++)
      {
        fAcc               = pfMatr[k*nDim + i];
        pfMatr[k*nDim + i] = pfMatr[m*nDim + i];
        pfMatr[m*nDim + i] = fAcc;
      }
      fAcc = pfVect[k];
      pfVect[k] = pfVect[m];
      pfVect[m] = fAcc;
    }

    if( pfMatr[k*nDim + k] == 0.) return 1; // needs improvement !!!

    // triangulation of matrix with coefficients
    for(j=(k+1); j<nDim; j++) // current row of matrix
    {
      fAcc = - pfMatr[j*nDim + k] / pfMatr[k*nDim + k];
      for(i=k; i<nDim; i++)
      {
        pfMatr[j*nDim + i] = pfMatr[j*nDim + i] + fAcc*pfMatr[k*nDim + i];
      }
      pfVect[j] = pfVect[j] + fAcc*pfVect[k]; // free member recalculation
    }
  }

  for(k=(nDim-1); k>=0; k--)
  {
    pfSolution[k] = pfVect[k];
    for(i=(k+1); i<nDim; i++)
    {
      pfSolution[k] -= (pfMatr[k*nDim + i]*pfSolution[i]);
    }
    pfSolution[k] = pfSolution[k] / pfMatr[k*nDim + k];
  }

  return 0;
}

int agbnp_fill_ctable(int n, float_a *x, float_a *y, float_a *yp,
		      C1Table *c1table){
  float_a *data;
  int i;
  float_a am[4*4],bv[4], f, xv[4];
  int pv[4], info, idx;
  int nam=4,nrhs=1,lda=4,ldb=4;
  int res;

  for(i=0;i<(n-1);i++){
    bv[0] = y[i];
    bv[1] = y[i+1];
    bv[2] = yp[i];
    bv[3] = yp[i+1];

    idx = 0;
    am[idx++] = pow(x[i],3);
    am[idx++] = pow(x[i],2);
    am[idx++] = x[i];
    am[idx++] = 1.;

    am[idx++] = pow(x[i+1],3);
    am[idx++] = pow(x[i+1],2);    
    am[idx++] = x[i+1];
    am[idx++] = 1.;

    am[idx++] = 3.*x[i]*x[i];
    am[idx++] = 2.*x[i];
    am[idx++] = 1.;
    am[idx++] = 0.;

    am[idx++] = 3.*x[i+1]*x[i+1];
    am[idx++] = 2.*x[i+1];
    am[idx++] = 1.;
    am[idx++] = 0.;


#ifdef NOTNOW
    {
      int j,k,idx;
      idx = 0;
      for(j=0;j<4;j++){
	for(k=0;k<4;k++){
	  printf("%f ",am[idx++]);
	}
	printf("   %f\n",bv[j]);
      }
    }
#endif

    /* printf("agbnp_fill_ctable(2): n=%d c1table=%p\n",n,c1table); */

    res = LinearEquationsSolving(nam,am,bv,xv);
    if(res){
      fprintf(stderr, "agbnp_fill_ctable(): error in LinearEquationsSolving()\n");
      return AGBNP_ERR;
    }
#ifdef NOTNOW
    {
      int j,k,idx;
      idx = 0;
      for(j=0;j<4;j++){
	for(k=0;k<4;k++){
	  printf("%f ",am[idx++]);
	}
	printf("   %f %f\n",bv[j],xv[j]);
      }
    }
#endif
    bv[0] = xv[0];
    bv[1] = xv[1];
    bv[2] = xv[2];
    bv[3] = xv[3];
    data = c1table->data[i];
    data[0] = bv[0]; /* a3 */
    data[1] = bv[1]; /* a2 */
    data[2] = bv[2]; /* a1 */
    data[3] = bv[3]; /* a0 */
  }

  return AGBNP_OK;
}


int agbnp_create_ctablef4(int n, float_a amax, float_a b, 
			  C1Table **c1table){
  C1Table *tbl;
  int i;
  float_a da = amax/(n-1), dx, dxinv;
  float_a u, a, u1, q1, u2, q2, dr, Rj=1.0, x1, x2, f1, fp, fpp, s;
  float_a ca3, ca2, ca1, ca0;
  float_a y2, y1, y1p, y2p;
  float_a *x, *y, *yp;
  float_a yinf = 1.0;
  float_a *data;

  x =  (float_a *)malloc(n*sizeof(float_a));
  y =  (float_a *)malloc(n*sizeof(float_a));
  yp =  (float_a *)malloc(n*sizeof(float_a));
  if(!(x && y && yp)){
    fprintf(stderr, "agbnp_create_ctablef4(): unable to allocate work buffers (%d floats)\n", 3*n);
    return AGBNP_ERR;
  }

  tbl = malloc(sizeof(C1Table));
  if(!tbl){
    fprintf(stderr, "agbnp_create_ctablef4(): unable to allocate table structure.\n");
    return AGBNP_ERR;
  }
  tbl->n = n;
  tbl->data = malloc(4*n*sizeof(float_a));
  if(!tbl->data){
    fprintf(stderr, "agbnp_create_ctablef4(): unable to allocate storage for data array (%d floats)\n",4*n);
    return AGBNP_ERR;
  }
  tbl->dx = da;
  tbl->dxinv = 1./da;
  tbl->yinf = yinf;

  a = 0.0;
  i = 0;
  x[i] = a;
  y[i] = 0.0;
  yp[i] = 0.0;
  a += da;
  for(i=1;i<n-1;i++){
    x[i] = a;
    u1 = 3.*Rj*pow(a,4)/(4.*pi);
    q1 = i4(a,b,Rj,&dr);
    s = 1. - agbnp_pol_switchfunc(a, F4LOOKUP_SMOOTHA, amax, &fp, &fpp);
    y[i] = 1.+ s*(u1*q1-1.);
    yp[i] = -fp*(u1*q1-1.) + s*(3./(4.*pi))*pow(a,3)*(4.*q1+a*dr);
    a += da;
  }
  x[n-1] = a;
  y[n-1] = yinf;
  yp[n-1] = 0.0;
  if(agbnp_fill_ctable(n, x, y, yp, tbl) != AGBNP_OK){
    fprintf(stderr, "agbnp_create_ctablef4(): error in agbnp_fill_ctable()\n");
    return AGBNP_ERR;
  }

  *c1table = tbl;
  free(x); free(y); free(yp);
  return AGBNP_OK;
}

int agbnp_create_ctablef42d(int na, float_a amax, 
			    int nb, float_a bmax, 
			    C1Table2D **table2d){
  C1Table2D *tbl2d;
  float_a b, bmin = 1./bmax, db;
  float_a a, u, q, dr, Rj=1.0;
  int i;
  float_a bmax1;
  
  db = bmax/(nb-1);
  bmax1 = bmax + 2.*db;    /* add buffer to bmax */
  db = bmax/(nb-1);
  

  /*  printf("agbnp_create_ctablef42d: na=%d, amax=%f, nb=%d, bmax=%f\n",na, amax, nb, bmax); */

  tbl2d = malloc(sizeof(C1Table2D));
  tbl2d->ny = nb;
  tbl2d->dy = db;
  tbl2d->dyinv = 1./db;
  tbl2d->table = (C1Table **)malloc(nb*sizeof(C1Table *));
  i = 0;
  b = db;
  if(agbnp_create_ctablef4(na,amax,b,&(tbl2d->table[i]))!=AGBNP_OK){
    fprintf(stderr, "agbnp_create_ctablef42d(): error in agbnp_create_ctablef4()\n");
    return AGBNP_ERR;
  }
  for(i=1;i<nb;i++){
    if(agbnp_create_ctablef4(na,amax,b,&(tbl2d->table[i]))!=AGBNP_OK){
      fprintf(stderr, "agbnp_create_ctablef42d(): error in agbnp_create_ctablef4()\n");
      return AGBNP_ERR;
    }
    b += db;
  }
  *table2d = tbl2d;
  return AGBNP_OK;
}

int agbnp_interpolate_ctable(C1Table *c1table, float_a x, 
			     float_a *f, float_a *fp){
  int i;
  float_a *data;
  float_a u, up, ca0, ca1, ca2, ca3;

  i = c1table->dxinv*x;
  if(i>(c1table->n-2)) {
    *f = c1table->yinf;
    *fp = 0.0;
    return AGBNP_OK;
  }
  data = c1table->data[i];
  ca3 = data[0];
  ca2 = data[1];
  ca1 = data[2];
  ca0 = data[3];
  u = ca3;       up = 3.*ca3;
  u = u*x + ca2; up = up*x + 2.*ca2;
  u = u*x + ca1; up = up*x + ca1;
  u = u*x + ca0;
  *f = u; *fp = up;
  return AGBNP_OK;
}

int agbnp_interpolate_ctablef42d(C1Table2D *table2d, float_a x, float_a y,
				 float_a *f, float_a *fp){
  int i, n, iy, j1, j2;
  float_a u1, u2, v2, v1, x2, y2, *values1, *values2; 
  float_a dx, dxinv, dy, dyinv;
  float_a a, b, f1, f2, fp1, fp2;
  C1Table *table1, *table2;

  dyinv = table2d->dyinv;
  dy = table2d->dy;
  iy = dyinv*y;
  if(iy>table2d->ny-2) { return AGBNP_ERR;}

  table1 = table2d->table[iy];
  table2 = table2d->table[iy+1];
  agbnp_interpolate_ctable(table1, x, &f1, &fp1);
  agbnp_interpolate_ctable(table2, x, &f2, &fp2);

  y2 = (iy+1)*dy;

  a = dyinv*(f2-f1);
  b = f2 - a*y2;
  *f = a*y + b;

  a = dyinv*(fp2-fp1);
  b = fp2 - a*y2;
  *fp = a*y + b;

  return AGBNP_OK;
}

/* initializes i4p(), the lookup table version of i4 */
int agbnp_init_i4p(void){
  if(agbnp_f4c1table2d != NULL){
    return AGBNP_OK;
  }
  if(agbnp_create_ctablef42d(F4LOOKUP_NA, F4LOOKUP_MAXA,
			     F4LOOKUP_NB, F4LOOKUP_MAXB, 
			     &agbnp_f4c1table2d) != AGBNP_OK){
    fprintf(stderr,"agbnp_init_i4p(): error in agbnp_create_ctablef42d()\n");
    return AGBNP_ERR;
  }
  return AGBNP_OK;
}

float_a i4p(float_a rij, float_a Ri, float_a Rj, 
	    float_a *dr){
  float_a a, b, rj1, f, fp, u, ainv, ainv2;
  static const float_a pf = (4.*pi)/3.;

  rj1 = 1./Rj;
  a = rj1*rij;
  ainv = 1./a;
  b = rj1*Ri;
  ainv2 = ainv*ainv;
  u = pf*rj1*ainv2*ainv2;
  agbnp_interpolate_ctablef42d(agbnp_f4c1table2d, a, b, &f, &fp);
  *dr = rj1*u*(fp-4.*ainv*f);

  return u*f;
}

static float_a *fextval;
#pragma omp threadprivate(fextval)
/* float comparison function for qsort */
int agbnp_fcompare( const void *val1, const void *val2 ) {
  int indx1 = *(int*)val1;
  int indx2 = *(int*)val2;

  /* primary key, sort by fextval */
  if( fextval[indx1] < fextval[indx2] ) {
    return -1;
  } else if( fextval[indx1] > fextval[indx2] ) {
    return 1;
  }

  /* secondary key, if fextvals identical, sort by index number */
  if( indx1 < indx2 )
    return -1;
  else if( indx1 > indx2 )
    return 1;

  return 0;
}

/* index sort of a float array */
void agbnp_fsortindx( int pnval, float_a *val, int *indx ) {
  int i;
  for (i=0; i< pnval; i++) indx[i] = i;
  fextval = val;
  qsort( indx, pnval, sizeof(int), agbnp_fcompare );
}

/* reorder neighbor list of atom iat based on mapping index indx */
int agbnp_nblist_reorder(AGBworkdata *agbw, NeighList *nl, int iat, int *indx){
  int j, nn = nl->nne[iat];
  int *js = agbw->js;
  NeighVector  *pbcs = agbw->pbcs;
  void **datas = agbw->datas;
  int error = 0;

  for(j=0;j<nn;j++){
    js[j] = nl->neighl[iat][indx[j]];
  }
  for(j=0;j<nn;j++){
    nl->neighl[iat][j] = js[j];
  }

  /* reorder pbc_trans */
  if(nl->pbc_trans){

    for(j=0;j<nn;j++){
      pbcs[j] = nl->pbc_trans[iat][indx[j]];
    }
    for(j=0;j<nn;j++){
      nl->pbc_trans[iat][j] = pbcs[j];
    }
  }

  /* reorder data */
  if(nl->data){

    for(j=0;j<nn;j++){
      datas[j] = nl->data_index[iat][indx[j]];
    }
    for(j=0;j<nn;j++){
      nl->data_index[iat][j] = datas[j];
    }
  }

  return AGBNP_OK;
}

/* reorder neighbor list of atom iat based on mapping index indx */
int agbnp_nblist_reorder_cavity( CAVworkdata *cavw, NeighList *nl, int iat, int *indx){
  int j, nn = nl->nne[iat];
  int *js = cavw->js;
  NeighVector  *pbcs = cavw->pbcs;
  void **datas = cavw->datas;
  int error = 0;

  for(j=0;j<nn;j++){
    js[j] = nl->neighl[iat][indx[j]];
  }
  for(j=0;j<nn;j++){
    nl->neighl[iat][j] = js[j];
  }

  /* reorder pbc_trans */
  if(nl->pbc_trans){

    for(j=0;j<nn;j++){
      pbcs[j] = nl->pbc_trans[iat][indx[j]];
    }
    for(j=0;j<nn;j++){
      nl->pbc_trans[iat][j] = pbcs[j];
    }
  }

  /* reorder data */
  if(nl->data){

    for(j=0;j<nn;j++){
      datas[j] = nl->data_index[iat][indx[j]];
    }
    for(j=0;j<nn;j++){
      nl->data_index[iat][j] = datas[j];
    }
  }

  return AGBNP_OK;
}


#endif
/* binds #ifndef USE_PBC */


/*
   routines below have USE_PBC and USE_ALT conditional compilations to
   to build slow and fast routine versions. 
*/

#ifdef USE_PBC
/* the PBC version is identical to the standard version. It's repeated
   so that inlining can work */
float_a i4_pbc(
#else
float_a i4(
#endif
	  float_a rij, float_a Ri, float_a Rj, 
	  float_a *dr){
  float_a u1, u2, u3,u4,u5,u6;
  float_a u4sq,u5sq;
  float_a rij2 = rij*rij;
  float_a q;
  static const float_a twopi = 2.0*pi;

  if(rij>(Ri+Rj)){
    u1 = rij+Rj;
    u2 = rij-Rj;
    u3 = u1*u2;
    u4 = 0.5*log(u1/u2);
    q = twopi*(Rj/u3 - u4/rij);
    *dr =  twopi*( (Rj/(rij*u3))*(1. - 2.*rij2/u3 ) + u4/rij2 );
  }else{
    u1 = Rj-Ri;
    if (rij2 > u1*u1){
      /* overlap */
      u1 = rij+Rj;
      u2 = rij-Rj;
      u3 = u1*u2;
      u4 = 1./u1;
      u4sq = u4*u4;
      u5 = 1./Ri;
      u5sq = u5*u5;
      u6 = 0.5*log(u1/Ri);
      q = twopi*(-(u4-u5) + (0.25*u3*(u4sq-u5sq) - u6)/rij);
      *dr = twopi*(0.5*(1.-0.5*u3/rij2)*(u4sq-u5sq)+u6/rij2);
    }else{
      /* inclusion */
      if(Ri>Rj){
	q = 0.0;
	*dr = 0.0;
      }else{
	u1 = rij+Rj;
	u2 = Rj - rij;
	u3 = -u1*u2; /* rij^2 - Rj^2 */
	u6 =  0.5*log(u1/u2);
	q = twopi*(2./Ri + Rj/u3 - u6/rij); 
	*dr = twopi*(-(Rj/u3)*(2.*rij/u3 - 1./rij) + u6/rij2);
      }
    }

  }

  return q;
}

#ifdef USE_PBC
/* the PBC/ALT version is identical to the standard version. It's repeated
   so that inlining can work */
 float_a ogauss_incremental_pbc(
#else
 float_a ogauss_incremental(
#endif
			   int n, float_a (*x)[3], float_a *a, float_a *p,
			   float_a *r,
			   float_a anm1, float_a pnm1, float_a cnm1[3],
			   float_a *an, float_a *pn, float_a cn[3],
			   float_a volmina, float_a volminb,
			   float_a (*dr)[3], float_a *dR, 
			   float_a (*d2rR)[AGBNP_MAX_OVERLAP_LEVEL][3]){
  /* calculated n-order gaussian overlap given 
     a(n-1): gaussian exponent of overlap if first n-1 gaussians
     p(n-1): gaussian prefactor of overlap if first n-1 gaussians
     c(n-1): center of overlap gaussian of first n-1 gaussians

     The new gaussian to overlap to the n-1 gaussian overlap is defined by the
     n-th member of x, a, p, that is x[n-1][3], a[n-1] and p[n-1]

     Returns the a, p, and c (gaussian center) of the nth-order 
     overlap gaussian
  */
  int i,j;
  float_a delta, deltai, kappa, vol, volp, fp, fpp;
  register float_a v0,v1,v2,d2,u;
  register float_a v, u1, u2;

  if(n<2){
    return 0.0;
  }

  /* new gaussian exponent */
  *an = delta = anm1 + a[n-1];
  deltai = 1./delta;

  /* new center */
  cn[0] = deltai*(anm1*cnm1[0] + a[n-1]*x[n-1][0]);
  cn[1] = deltai*(anm1*cnm1[1] + a[n-1]*x[n-1][1]);
  cn[2] = deltai*(anm1*cnm1[2] + a[n-1]*x[n-1][2]);

  /* square distance from old center */
  v0 = x[n-1][0] - cnm1[0];
  v1 = x[n-1][1] - cnm1[1];
  v2 = x[n-1][2] - cnm1[2];
  d2 = v0*v0 + v1*v1 + v2*v2;
  
  /* kappa */
#ifdef USE_SSE
  kappa = expbf(-anm1*a[n-1]*d2*deltai);
#else
  kappa = exp(-anm1*a[n-1]*d2*deltai);
#endif
  /* new prefactor */
  *pn = pnm1*p[n-1]*kappa;

  /* gaussian overlap volume */
  u = pi*deltai;
  vol = (*pn)*u*sqrt(u);
  volp = swf_vol3(vol, &fp, &fpp, volmina, volminb);
  /* if volume is less than minimum volume can quit now */
  if(volp < FLT_MIN){
    /* note that if ogauss_incremental returns 0 the derivatives
       are undefined */
    return 0.0;
  }


  /*                      *
   *        derivatives   *
   *                      */

  if(dr && dR){

    /* positional and radius derivatives */
    for(i=0;i<n;i++){
      u = -2.*a[i]*vol;
      dr[i][0] = v0 = u*(x[i][0]-cn[0]);
      dr[i][1] = v1 = u*(x[i][1]-cn[1]);
      dr[i][2] = v2 = u*(x[i][2]-cn[2]);
      d2 = v0*v0 + v1*v1 + v2*v2;
      u = a[i]*vol;
      dR[i] = (3.*u*deltai + 0.5*d2/u)/r[i];
    }

    if(d2rR){

      /* pos-radius derivatives */
      /* diagonal terms */
      for(i=0;i<n;i++){
	u = -2.*fp*(1.-a[i]*deltai)/r[i] + (fpp+fp/vol)*dR[i];
	d2rR[i][i][0] = u*dr[i][0];
	d2rR[i][i][1] = u*dr[i][1];
	d2rR[i][i][2] = u*dr[i][2];
      }
      /* cross terms */
      for(i=0;i<n;i++){
	for(j=i+1;j<n;j++){
	  u = fpp+fp/vol;
	  v = 2.*fp*deltai;
	  u1 = v*a[i]/r[j];
	  u2 = u*dR[j];
	  d2rR[i][j][0] = u1*dr[j][0]+u2*dr[i][0];
	  d2rR[i][j][1] = u1*dr[j][1]+u2*dr[i][1];
	  d2rR[i][j][2] = u1*dr[j][2]+u2*dr[i][2];
	  u1 = v*a[j]/r[i];
	  u2 = u*dR[i];
	  d2rR[j][i][0] = u1*dr[i][0]+u2*dr[j][0];
	  d2rR[j][i][1] = u1*dr[i][1]+u2*dr[j][1];
	  d2rR[j][i][2] = u1*dr[i][2]+u2*dr[j][2];
	}
      }
      
    }
  
    /* switch first derivatives */
    for(i=0;i<n;i++){
      dr[i][0] *= fp;
      dr[i][1] *= fp;
      dr[i][2] *= fp;            
      dR[i] *= fp;
    }
    
  }
  return volp;
}



/* calculates inverse Born radii */
#ifdef USE_PBC
int agbnp_inverse_born_radii_pbc(
#else
int agbnp_inverse_born_radii(
#endif
				 AGBNPdata *agb, AGBworkdata *agbw_h,
				    float_a *x, float_a *y, float_a *z,
				    int init_frozen){
  float_a fourpi1 = 1./(4.*pi);
  int iq4cache=0, iv2cache = 0; /* cache counters */
  int i, iat, j, jat;
  float_a volume2, spiat, spjat, dx, dy, dz, d, q, dr4;
  float_a riat, rjat;
  float_a altw = 1.0, awi = 1.0, awj = 1.0;
  /* needed for I4_FAR() macro */
  float_a _i4u1, _i4u2, _i4u3, _i4u4;
  float_a _i4rij2;
  float_a _i4dr4;
  float_a _i4twopi = 2.0*pi;

  int natoms = agb->natoms;
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  int *isheavy = agbw_h->isheavy;
  int nhydrogen = agb->nhydrogen;
  int *ihydrogen = agb->ihydrogen;
  NeighList *near_nl = agbw_h->near_nl;
  NeighList *far_nl = agbw_h->far_nl;
  float_a *r = agb->r;
  float_a *sp = agbw_h->sp;
  float_a *vols = agbw_h->vols;
  float_a *br1 = agbw_h->br1;
  float_a *br1_const = agb->agbw->br1_const;
  int *isvfrozen = agb->agbw->isvfrozen;
  float *v2cache = agbw_h->v2cache;
  float *q4cache = agbw_h->q4cache;
  int do_frozen = agb->do_frozen;
#ifdef USE_PBC
  float_a *xj = agbw_h->xj;
  float_a *yj = agbw_h->yj;
  float_a *zj = agbw_h->zj;
#endif
#ifdef USE_ALT
  int doalt = agb->doalt;
  float_a *occupancy = agb->occupancy;
#endif

  /* reset inverse born radii */
  memset(br1, 0, natoms*sizeof(float_a));
  if(init_frozen){
#pragma omp single
    memset(br1_const,0,natoms*sizeof(float_a));
#pragma omp barrier
  }

  /* Loop over near neighbors first, these need scaled volume correction */
  for(i=0;i<nheavyat;i++){
    iat = iheavyat[i];
#ifdef USE_PBC
    if(agbnp_get_neigh_coords(agb, iat, near_nl, xj, yj, zj) != AGBNP_OK){
      fprintf(stderr, "agbnp_inverse_born_radii(): error in agbnp_get_neigh_coords()\n");
      return AGBNP_ERR;
    } 
#endif
    for(j= 0; j< near_nl->nne[iat]; j++){
      jat = near_nl->neighl[iat][j] % natoms;
      if(do_frozen){
	/* skip vfrozen pairs */
	if(isvfrozen[iat] && isvfrozen[jat]) continue;
      }
#ifdef USE_ALT
      /* alt. conf weight */
      if(doalt){
	altw = altweight_pair(agb, iat, jat);
	awi = altw/occupancy[iat];
	awj = altw/occupancy[jat];
      }
#endif
      /* correct scaled volumes */
      volume2 = v2cache[iv2cache++];
      spiat = sp[iat] + 0.5*volume2/vols[iat];
      spjat = sp[jat] + 0.5*volume2/vols[jat];
      /* computes q4 integrals and subtract them from inverse Born radius */
#ifdef USE_PBC
      dx = xj[j] - x[iat];
      dy = yj[j] - y[iat];
      dz = zj[j] - z[iat];
#else
      dx = x[jat] - x[iat];
      dy = y[jat] - y[iat];
      dz = z[jat] - z[iat];
#endif
      d = mysqrt(dx*dx + dy*dy + dz*dz);
#ifdef USE_PBC
      q = awj*i4_pbc(d,r[jat],r[iat], &dr4);
#else
      q = i4p(d,r[jat],r[iat], &dr4);
#endif
      br1[jat] -= fourpi1*q*spiat;
      q4cache[iq4cache++] = q;
      q4cache[iq4cache++] = awj*dr4;
      if(init_frozen){
	if(isvfrozen[iat] && isvfrozen[jat]){
#ifdef _OPENMP
#pragma omp atomic
#endif
	  br1_const[jat] -= fourpi1*q*spiat;
	}
      }
#ifdef USE_PBC
      q = awi*i4_pbc(d,r[iat],r[jat], &dr4);
#else
      q = i4p(d,r[iat],r[jat], &dr4);
#endif
      br1[iat] -= fourpi1*q*spjat;
      q4cache[iq4cache++] = q;
      q4cache[iq4cache++] = awi*dr4;
      if(init_frozen){
	if(isvfrozen[iat] && isvfrozen[jat]){
#ifdef _OPENMP
#pragma omp atomic
#endif
	  br1_const[iat] -= fourpi1*q*spjat;
	}
      }
    }
  }

  /* loop over far heavy atom neighbors */
  for(i = 0; i < nheavyat; i++){
    iat = iheavyat[i];
    riat = r[iat];
#ifdef USE_PBC
    if(agbnp_get_neigh_coords(agb, iat, far_nl, xj, yj, zj) != AGBNP_OK){
      fprintf(stderr, "agbnp_inverse_born_radii(): error in agbnp_get_neigh_coords()\n");
      return AGBNP_ERR;
    } 
#endif
    for(j= 0; j< far_nl->nne[iat]; j++){
      jat = far_nl->neighl[iat][j] % natoms;
      if(do_frozen){
	if(isvfrozen[iat] && isvfrozen[jat]) continue;
      }
#ifdef USE_ALT
      if(doalt){
	altw = altweight_pair(agb, iat, jat);
	awi = altw/occupancy[iat];
	awj = altw/occupancy[jat];
      }
#endif
#ifdef USE_PBC
      dx = xj[j] - x[iat];
      dy = yj[j] - y[iat];
      dz = zj[j] - z[iat];
#else
      dx = x[jat] - x[iat];
      dy = y[jat] - y[iat];
      dz = z[jat] - z[iat];
#endif
      d = mysqrt(dx*dx + dy*dy + dz*dz);
      rjat = r[jat];
      if(isheavy[jat]){              /* Born radii change for iat and jat*/
	q = awj*I4_FAR(d,rjat,riat);
	br1[jat] -= fourpi1*q*sp[iat];
	q4cache[iq4cache++] = q;
	q4cache[iq4cache++] = awj*_i4dr4;
	if(init_frozen){
	  if(isvfrozen[iat] && isvfrozen[jat]){
#ifdef _OPENMP
#pragma omp atomic
#endif
	    br1_const[jat] -= fourpi1*q*sp[iat];
	  }
	}
	q = awi*I4_FAR(d,riat,rjat);
	br1[iat] -= fourpi1*q*sp[jat];
	q4cache[iq4cache++] = q;
	q4cache[iq4cache++] = awi*_i4dr4;
	if(init_frozen){
	  if(isvfrozen[iat] && isvfrozen[jat]){
#ifdef _OPENMP
#pragma omp atomic
#endif
	    br1_const[iat] -= fourpi1*q*sp[jat];
	  }
	}
      }else{           /* jat is hydrogen, Born radii change only for jat */
#ifdef USE_PBC
	q = awj*i4_pbc(d,rjat,riat, &dr4);
#else
	q = i4p(d,rjat,riat, &dr4);
#endif
	br1[jat] -= fourpi1*q*sp[iat];
	q4cache[iq4cache++] = q;
	q4cache[iq4cache++] = awj*dr4;
	if(init_frozen){
	  if(isvfrozen[iat] && isvfrozen[jat]){
#ifdef _OPENMP
#pragma omp atomic
#endif
	    br1_const[jat] -= fourpi1*q*sp[iat];
	  }
	}
      }
    }
  }

  for(i = 0; i < nhydrogen; i++){ /* loop for hydrogen-heavy interactions */
    iat = ihydrogen[i];
#ifdef USE_PBC
    if(agbnp_get_neigh_coords(agb, iat, far_nl, xj, yj, zj) != AGBNP_OK){
      fprintf(stderr, "agbnp_inverse_born_radii(): error in agbnp_get_neigh_coords()\n");
      return AGBNP_ERR;
    } 
#endif
    for(j = 0; j < far_nl->nne[iat]; j++){
      jat = far_nl->neighl[iat][j] % natoms;
      if(do_frozen){
	if(isvfrozen[jat] && isvfrozen[iat]) continue;
      }
      if(isheavy[jat]){
#ifdef USE_ALT 
	if(doalt){
	  altw = altweight_pair(agb, iat, jat);
	  awi = altw/occupancy[iat];
	}
#endif
#ifdef USE_PBC
	dx = xj[j] - x[iat];
	dy = yj[j] - y[iat];
	dz = zj[j] - z[iat];
#else
	dx = x[jat] - x[iat];
	dy = y[jat] - y[iat];
	dz = z[jat] - z[iat];
#endif
	d = mysqrt(dx*dx + dy*dy + dz*dz);
#ifdef USE_PBC
	q = awi*i4_pbc(d,r[iat],r[jat], &dr4);
#else
	q = i4p(d,r[iat],r[jat], &dr4);
#endif
	br1[iat] -= fourpi1*q*sp[jat];
	q4cache[iq4cache++] = q;
	q4cache[iq4cache++] = altw*dr4;
	if(init_frozen){
	  if(isvfrozen[jat] && isvfrozen[iat]){
#ifdef _OPENMP
#pragma omp atomic
#endif
	    br1_const[iat] -= fourpi1*q*sp[jat];
	  }
	}
      }
    }
  }
#pragma omp barrier

  return AGBNP_OK;
}


#if defined(USE_SSE) && !defined(USE_PBC)
float pow2lt1(float x){
  float emz;
  int j;
  static const int n=8;
  emz = twoexpa[n-1];
  for(j=n-2;j>=0;j--){
    emz = emz*x + twoexpa[j];
  }
  return emz;
}

float expbf(float y){
  unsigned int n, tn;
  float dx;
  int inv = 0;
  static const float xmax = 31.9999;
  float x = ln2inv*y;
  float u;

  if(x < 0){
    x = -x;
    inv = 1;
  }
  if(x > xmax) {x = xmax;}
  n = x;
  tn = (1<<n);
  dx = x - n;
  u = tn*pow2lt1(dx);
  if(inv){
    return 1./u;
  }
  return u;
}


__m128 vexpf(__m128 vy){
  unsigned int n[4], tn[4];
  __m128 vx, vtn, vn,vdx, vemz;
  float *sx = (float *)&vx;
  float *sn = (float *)&vn;
  float *stn = (float *)&vtn;
  int j, nt=8;
  static const float xmax = 31.9999;

  vx =_mm_mul_ps(vln2inv,vy);
  if(sx[0]>xmax) {sx[0] = xmax;}
  if(sx[1]>xmax) {sx[1] = xmax;}
  if(sx[2]>xmax) {sx[2] = xmax;}
  if(sx[3]>xmax) {sx[3] = xmax;}
  n[0] = sx[0];
  n[1] = sx[1];
  n[2] = sx[2];
  n[3] = sx[3];
  sn[0] = n[0];
  sn[1] = n[1];
  sn[2] = n[2];
  sn[3] = n[3];
  stn[0] = (1<<n[0]);
  stn[1] = (1<<n[1]);
  stn[2] = (1<<n[2]);
  stn[3] = (1<<n[3]);
  vdx = _mm_sub_ps(vx,vn);
  vemz = vtwoexpa[nt-1];
  for(j=nt-2;j>=0;j--){
    vemz = _mm_mul_ps(vemz,vdx);
    vemz = _mm_add_ps(vemz,vtwoexpa[j]);
  }
  vemz = _mm_mul_ps(vtn,vemz);
  return vemz;
}
#endif

#ifdef USE_PBC
 int agbnp_gb_energy_inner_pbc(
#else
 int agbnp_gb_energy_inner(
#endif
		    AGBNPdata *agb, int iat, int natoms, 
		    int nb, int *nbiat, float_a *x, float_a *y, float_a *z,
		    float_a qiat, float_a biat,
		    float_a *charge, float_a *br, float_a *ab, float_a *dera,
		    float_a (*dgbdr)[3], float_a *egb_pair, 
		    float_a dielectric_factor,
		    int init_frozen, float_a (*dgbdr_c)[3],
		    float_a *gbpair_const, float_a *dera_const,
		    int *isbfrozen
#ifdef USE_PBC
		    , float_a *xj, float_a *yj, float_a *zj,
		    float_a (*rot)[3][3]
#endif
		    ){
   
   int j, jat, ib;
   float_a dx, dy, dz, d2, qjat, qq, ab2ij, bb, etij, fgb, gbij, atij, htij;
   float_a c1gb = 1./4.; 
   float_a u[3], ur[3];
#ifdef USE_PBC
   int isym;
#endif
#ifdef USE_ALT
   float_a altw, fsame;
#endif
#ifdef USE_SSE
   int j0,j1,j2,j3;
   __m128 vxi, vyi, vzi, vxj, vyj, vzj;
   __m128 vqiat, vqjat,vbiat,vbjat;
   float two = 2.0; __m128 vtwo;
   float pt25 = 0.25; __m128 vpt25;
   float pt5 = 0.5; __m128 vpt5;
   float one = 1.0; __m128 vone;
   float three = 3.0;__m128 vthree;
   float del = 2.*dielectric_factor;  __m128 vdiel;
   __m128 vr2, vetij, vu, vf, vfgb, vgbij, vatij, vhtij;
   __m128 vdx, vdy, vdz, vux, vuy, vuz;
   float *sgbij = (float *)&vgbij, *satij = (float *)&vatij, 
     *setij = (float *)&vetij, *sr2 = (float *)&vr2;
   float *sux = (float *)&vux, *suy = (float *)&vuy, *suz = (float *)&vuz;
   float se[4];
#endif

   ib = 0;
#if defined(USE_SSE) && !defined(USE_PBC)
   /* vectorized loop (no pbc for the moment)*/
   vtwo = _mm_set_ps(two,two,two,two);
   vpt25 = _mm_set_ps(pt25,pt25,pt25,pt25);
   vpt5 = _mm_set_ps(pt5,pt5,pt5,pt5); 
   vone = _mm_set_ps(one,one,one,one);
   vthree = _mm_set_ps(three,three,three,three);
   vdiel = _mm_set_ps(del,del,del,del);

   vxi = _mm_set_ps(x[iat],x[iat],x[iat],x[iat]);
   vyi = _mm_set_ps(y[iat],y[iat],y[iat],y[iat]);
   vzi = _mm_set_ps(z[iat],z[iat],z[iat],z[iat]);
   vqiat = _mm_set_ps(qiat,qiat,qiat,qiat);
   vbiat = _mm_set_ps(biat,biat,biat,biat);
   for(j=0;j<nb/4;j++){
     j0 = nbiat[ib++];
     j1 = nbiat[ib++];
     j2 = nbiat[ib++];
     j3 = nbiat[ib++];
     vxj = _mm_set_ps(x[j0],x[j1],x[j2],x[j3]);
     vyj = _mm_set_ps(y[j0],y[j1],y[j2],y[j3]);
     vzj = _mm_set_ps(z[j0],z[j1],z[j2],z[j3]);
     vqjat = _mm_set_ps(charge[j0],charge[j1],charge[j2],charge[j3]);
     vbjat = _mm_set_ps(br[j0],br[j1],br[j2],br[j3]);
     vdx = _mm_sub_ps(vxj,vxi);
     vdy = _mm_sub_ps(vyj,vyi);
     vdz = _mm_sub_ps(vzj,vzi);
     vxj = _mm_mul_ps(vdx,vdx);
     vyj = _mm_mul_ps(vdy,vdy);
     vzj = _mm_mul_ps(vdz,vdz);
     vxj = _mm_add_ps(vxj,vyj);
     vr2 = _mm_add_ps(vxj,vzj);/* r^2 */
     vqjat = _mm_mul_ps(vqiat,vqjat);
     vbjat = _mm_mul_ps(vbiat,vbjat);

     vetij = _mm_mul_ps(vpt25,vr2);
     vetij = _mm_div_ps(vetij,vbjat);

     vetij = vexpf(vetij);
     vetij = _mm_rcp_ps(vetij);

     vu = _mm_mul_ps(vbjat,vetij);
     vu = _mm_add_ps(vr2,vu);

     vf = _mm_rsqrt_ps(vu); /* 1/sqrt(...) first guess */
     /* 1 step Newton-Raphson vfgb = pt5*vf*(three-vf*vf*vu) */
     vfgb = _mm_mul_ps(vf,vf);
     vfgb = _mm_mul_ps(vfgb,vu);
     vfgb = _mm_sub_ps(vthree,vfgb);
     vfgb = _mm_mul_ps(vf,vfgb);
     vfgb = _mm_mul_ps(vpt5,vfgb);

     vu = _mm_mul_ps(vtwo,vfgb);
     vgbij = _mm_mul_ps(vqjat,vu);

     *egb_pair += sgbij[0] +sgbij[1] +sgbij[2] +sgbij[3];
     /* Ai's */
     vu = _mm_mul_ps(vfgb,vfgb);
     vfgb = _mm_mul_ps(vu,vfgb);

     vatij = _mm_mul_ps(vpt25,vr2);
     vatij = _mm_add_ps(vbjat, vatij);
     vatij = _mm_mul_ps(vqjat,vatij);
     vatij = _mm_mul_ps(vetij, vatij);
     vatij = _mm_mul_ps(vfgb, vatij);

     dera[iat] += satij[0] + satij[1] + satij[2] + satij[3];
     dera[j0] += satij[3];
     dera[j1] += satij[2];
     dera[j2] += satij[1];
     dera[j3] += satij[0];
     
     /* derivative contribution if Born radii are constant */
     vhtij = _mm_mul_ps(vpt25,vetij);
     vhtij = _mm_sub_ps(vone,vhtij);
     vhtij = _mm_mul_ps(vqjat,vhtij);
     vhtij = _mm_mul_ps(vfgb,vhtij);
     vhtij = _mm_mul_ps(vdiel,vhtij); /* vdiel includes 2 scale factor */

     vux = _mm_mul_ps(vhtij,vdx);
     vuy = _mm_mul_ps(vhtij,vdy);
     vuz = _mm_mul_ps(vhtij,vdz);

     dgbdr[iat][0] +=   sux[3]+sux[2]+sux[1]+sux[0];
     dgbdr[iat][1] +=   suy[3]+suy[2]+suy[1]+suy[0];
     dgbdr[iat][2] +=   suz[3]+suz[2]+suz[1]+suz[0];
     
     dgbdr[j0][0] -=   sux[3];
     dgbdr[j0][1] -=   suy[3];
     dgbdr[j0][2] -=   suz[3];

     dgbdr[j1][0] -=   sux[2];
     dgbdr[j1][1] -=   suy[2];
     dgbdr[j1][2] -=   suz[2];

     dgbdr[j2][0] -=   sux[1];
     dgbdr[j2][1] -=   suy[1];
     dgbdr[j2][2] -=   suz[1];

     dgbdr[j3][0] -=   sux[0];
     dgbdr[j3][1] -=   suy[0];
     dgbdr[j3][2] -=   suz[0];

   }
#endif
   /* regular loop */
   for(j=ib;j<nb;j++){
     jat = nbiat[j];
#ifdef USE_PBC
     isym = jat/natoms;
     jat = jat % natoms;
#endif
#ifdef USE_ALT
     altw = altweight_pair(agb, iat, jat);
#endif
#ifdef USE_PBC
     dx = xj[j] - x[iat];
     dy = yj[j] - y[iat];
     dz = zj[j] - z[iat];
#else
     dx = x[jat] - x[iat];
     dy = y[jat] - y[iat];
     dz = z[jat] - z[iat];
#endif
     d2 = dx*dx + dy*dy + dz*dz;
     qjat = charge[jat];
#ifdef USE_ALT
     fsame = (jat==iat) ? 0.5 : 1.0;
     qq = fsame*altw*qiat*qjat;
#else
     qq = qiat*qjat;
#endif
#ifndef LORENTIAN_FIJ
     ab2ij = ab[iat]*ab[jat];
     bb = biat*br[jat];
     /*
     if(qq<0.0) {
       c1gb = 1./2.;
     }else{
       c1gb = 1./4.;
     }
     */
     etij = exp(-c1gb*d2/bb);
     fgb = 1./mysqrt(d2 + ab2ij*bb*etij);
     gbij = 2.*qq*fgb;
     *egb_pair += gbij;
     /* Ai's */
     fgb = fgb*fgb*fgb;
     atij = qq*ab2ij*(bb+0.25*d2)*etij*fgb;
     dera[iat] += atij;
     dera[jat] += atij;
     /* derivative contribution if Born radii are constant */
     htij = 2.*qq*(1.-ab2ij*0.25*etij)*fgb;
#else
     /* Modified effective GB pair distance
	fgbij = sqrt[rij^2 + Bij^2*g(xij)]
	where
	xij = rij/Bij
	and
	g(xij) = a^4/[ (a^4 - Sij^4)+(xij^2+Sij^2)^2]
	
	Sij = "screening factor"
	Sij = Si + Sj
     */
     bb = biat*br[jat];
     abij = fabs(ab[iat]*qiat)+fabs(ab[jat]*qjat);
     ab2ij = abij*abij;
     xij2 = d2/bb;
     xaij2 = xij2 + ab2ij;
     gij = a4/(a4-ab2ij*ab2ij+xaij2*xaij2);
     fgb = 1./mysqrt(d2+bb*gij);
     gbij = 2.*qq*fgb;
     *egb_pair += gbij;
     
     /* Ai's (Yi's in the paper)*/
     fgb = fgb*fgb*fgb;
     hij = 2.*gij*gij*xaij2/a4;
     atij = qq*bb*(gij+xij2*hij)*fgb;
     dera[iat] += atij;
     dera[jat] += atij;
     
     /* derivative contribution if Born radii are constant */
     htij = 2.*qq*(1.-hij)*fgb;
     
#endif
     htij = dielectric_factor*htij;
     u[0] = ur[0] = htij*dx;
     u[1] = ur[1] = htij*dy;
     u[2] = ur[2] = htij*dz;
#ifdef USE_PBC
     if(agb->docryst){
       /* rotate gradient */
       rtvec(ur,rot[isym],u);
     }
#endif
     dgbdr[iat][0] +=   u[0];
     dgbdr[iat][1] +=   u[1];
     dgbdr[iat][2] +=   u[2];
     dgbdr[jat][0] -=   ur[0];
     dgbdr[jat][1] -=   ur[1];
     dgbdr[jat][2] -=   ur[2];
     
     if(init_frozen){
       if(isbfrozen[iat] && isbfrozen[jat]){
#pragma omp critical
	 {
	   gbpair_const[iat] += gbij;
	   dera_const[iat] += atij;
	   dgbdr_c[iat][0] +=  u[0];
	   dgbdr_c[iat][1] +=  u[1];
	   dgbdr_c[iat][2] +=  u[2];
	   dgbdr_c[jat][0] -=  ur[0];
	   dgbdr_c[jat][1] -=  ur[1];
	   dgbdr_c[jat][2] -=  ur[2];
	 }
       }
     }
     
   }
   
   return AGBNP_OK;
}

/* Evaluates solvation energy, Ai's and derivatives of GB energy 
   at constant Born radii */
#ifdef USE_PBC
 int agbnp_gb_energy_pbc(
#else
 int agbnp_gb_energy(
#endif
		     AGBNPdata *agb, AGBworkdata *agbw_h,
		     float_a *x, float_a *y, float_a *z,
		     float_a *egb_selfs, float_a *egb_pairs,
		     int init_frozen){
  float_a egb_self = 0.0;
  float_a egb_pair = 0.0;
  int iat,j,jat;
#ifdef USE_ALT
  float_a fsame = 1.0, altw = 1.0;
#endif
  float_a qiat,qjat,biat,dx,dy,dz,d2,qq,ab2ij,bb,etij,fgb,atij,htij;
  float_a gbij;
  float_a dielectric_factor = 
    -0.5*(1./agb->dielectric_in - 1./agb->dielectric_out);
#ifdef LORENTIAN_FIJ
  static const float_a a4 = 9.0;
  float_a gij,hij,xij2,xaij2,abij;
#endif
  int natoms = agb->natoms;
  NeighList *near_nl = agbw_h->near_nl;
  NeighList *far_nl = agbw_h->far_nl;
  float_a *charge = agb->charge;
  float_a *ab = agb->ab;
  float_a *br = agbw_h->br;
  float_a *dera = agbw_h->dera;
  float_a (*dgbdr)[3] = agbw_h->dgbdr_h;
  int *isbfrozen = agb->agbw->isbfrozen;
  float_a *gbpair_const = agb->agbw->gbpair_const;
  float_a *dera_const = agb->agbw->dera_const;
  int do_frozen = agb->do_frozen;
  float_a *dera_m = agb->agbw->dera;
  float_a (*dgbdr_c)[3] = agb->agbw->dgbdr_c;
#ifdef USE_PBC
  float_a *xj = agbw_h->xj;
  float_a *yj = agbw_h->yj;
  float_a *zj = agbw_h->zj;
  float_a (*rot)[3][3] = agb->rot;
  int isym;
#endif
  float_a u[3], ur[3];


  if(init_frozen){
#pragma omp single
    memset(gbpair_const,0,natoms*sizeof(float_a));
#pragma omp barrier
  }

  for(iat=0;iat<natoms;iat++){
    qiat = charge[iat];
    biat = br[iat];
#pragma omp master
#ifdef USE_ALT
    egb_self += agb->occupancy[iat]*qiat*qiat/biat;
#else
    egb_self += qiat*qiat/biat;
#endif
#ifdef USE_PBC
    if(agbnp_get_neigh_coords(agb, iat, near_nl, xj, yj, zj) != AGBNP_OK){
      fprintf(stderr, "agbnp_gb_energy(): error in agbnp_get_neigh_coords()\n");
      return AGBNP_ERR;
    } 
#endif
#ifdef USE_PBC
    agbnp_gb_energy_inner_pbc(agb, iat, natoms, 
			      near_nl->nne[iat],  near_nl->neighl[iat],
			      x, y, z, qiat, biat,
			      charge, br, ab, dera, dgbdr, &egb_pair,
			      dielectric_factor,
			      init_frozen, dgbdr_c, gbpair_const, dera_const,
			      isbfrozen,
			      xj, yj, zj, rot);
			      
#else
    
    agbnp_gb_energy_inner(agb, iat, natoms, 
			  near_nl->nne[iat],  near_nl->neighl[iat],
			  x, y, z, qiat, biat,
			  charge, br, ab, dera, dgbdr, &egb_pair,
			  dielectric_factor,
			  init_frozen, dgbdr_c, gbpair_const, dera_const,
			  isbfrozen);

#endif

    /* loop over far neighbors */
#ifdef USE_PBC
    if(agbnp_get_neigh_coords(agb, iat, far_nl, xj, yj, zj) != AGBNP_OK){
      fprintf(stderr, "agbnp_gb_energy(): error in agbnp_get_neigh_coords()\n");
      return AGBNP_ERR;
    } 
#endif
#ifdef USE_PBC
    agbnp_gb_energy_inner_pbc(agb, iat, natoms, 
			      far_nl->nne[iat], far_nl->neighl[iat],
			      x, y, z, qiat, biat,
			      charge, br, ab, dera, dgbdr, &egb_pair,
			      dielectric_factor,
			      init_frozen, dgbdr_c, gbpair_const, dera_const,
			      isbfrozen,
			      xj, yj, zj, rot);
			      
#else
    
    agbnp_gb_energy_inner(agb, iat, natoms, 
			  far_nl->nne[iat],  far_nl->neighl[iat],
			  x, y, z, qiat, biat,
			  charge, br, ab, dera, dgbdr, &egb_pair,
			  dielectric_factor,
			  init_frozen, dgbdr_c, gbpair_const, dera_const,
			  isbfrozen);

#endif

  }

#ifdef _OPENMP

  /* reduce energies and dera */
#pragma omp single
  {
    *egb_selfs = 0.0;
    *egb_pairs = 0.0;
  }
#pragma omp barrier
#pragma omp critical
  {
    for(iat=0;iat<natoms;iat++){
      dera_m[iat] += dera[iat];
    }
    *egb_selfs += egb_self;
    *egb_pairs += egb_pair;
  }
#pragma omp barrier


#else

 *egb_selfs = egb_self;
 *egb_pairs = egb_pair;

#endif

 if(do_frozen){
   /* add constant terms to GB pair energy */
#pragma omp single
   {
     for(iat=0;iat<natoms;iat++){
       *egb_pairs += gbpair_const[iat];
       dera_m[iat] += dera_const[iat];
     }
   }
 }

#ifdef _OPENMP
 /* copy total dera back to the threads */
#pragma omp barrier
  memcpy(dera,dera_m,natoms*sizeof(float_a));
#endif


 /* compute some auxiliary arrays that depend on dera and 
    born radii */
#pragma omp barrier
#pragma omp single
  {
    float_a *q2ab = agb->agbw->q2ab;
    float_a *br1_swf_der = agb->agbw->br1_swf_der;
    float_a *alpha = agb->agbw->alpha;
    float_a *brw = agb->agbw->brw;
    float_a *abrw = agb->agbw->abrw;
    float_a *br_m = agb->agbw->br;
    /* q2ab[], brw[], abrw[], br1_swf_der[], and dera_m[] 
       are master copies */
    for(iat=0;iat<natoms;iat++){
      q2ab[iat] = agb->occupancy[iat]*charge[iat]*charge[iat]+dera_m[iat]*br_m[iat];
      q2ab[iat] *= br1_swf_der[iat];
      abrw[iat] = agb->occupancy[iat]*alpha[iat]*brw[iat];
      abrw[iat] *= br1_swf_der[iat];
    }
  }
#pragma omp barrier

#ifdef _OPENMP
  /* copy auxiliary arrays to threads */
  memcpy(agbw_h->q2ab,agb->agbw->q2ab,natoms*sizeof(float_a));
  memcpy(agbw_h->abrw,agb->agbw->abrw,natoms*sizeof(float_a));
#endif


  return AGBNP_OK;
}

#ifdef USE_PBC
 int agbnp_self_volumes_pbc(
#else
 int agbnp_self_volumes(
#endif
			AGBNPdata *agb, AGBworkdata *agbw,
				float_a *x, float_a *y, float_a *z){
  int iv2cache = 0;
  /* coordinate buffer for Gaussian overlap calculation */
  float_a gx[AGBNP_MAX_OVERLAP_LEVEL][3];
  /* radius buffer for  Gaussian overlap calculation */
  float_a gr[AGBNP_MAX_OVERLAP_LEVEL];
  /* Gaussian parameters buffers for overlap calculation */
  float_a ga[AGBNP_MAX_OVERLAP_LEVEL];
  float_a gp[AGBNP_MAX_OVERLAP_LEVEL];
  /* holds the atom indexes being overlapped */
  int gnlist[AGBNP_MAX_OVERLAP_LEVEL];
  int gatlist[AGBNP_MAX_OVERLAP_LEVEL];
  float_a gvol; /* gaussian overlap volume */
  int order; /* order of the overlap */
  /* holds overlap Gaussian parameters at each level */
  GParm gparams[AGBNP_MAX_OVERLAP_LEVEL];
  float_a an, pn, cn[3];
  /* coefficients for scaled volume */
  float_a volpcoeff[AGBNP_MAX_OVERLAP_LEVEL];
  float_a sign;

  int i, iat, jat, kat, ii,j;
  float_a u,v,altw=1.0;
#ifdef USE_ALT
  float_a gaw[AGBNP_MAX_OVERLAP_LEVEL];
#endif

  int natoms = agb->natoms;
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  float_a *occupancy = agb->occupancy;
  NeighList *near_nl = agbw->near_nl;
  float_a *r = agb->r;
  float_a *galpha = agbw->galpha;
  float_a *gprefac = agbw->gprefac;
  float_a *volumep = agbw->volumep;
  float *v2cache = agbw->v2cache;
  int do_frozen = agb->do_frozen;
  int *isvfrozen = agb->agbw->isvfrozen;
  int *isvvfrozen = agb->agbw->isvvfrozen;
#ifdef USE_PBC
  float_a *xj = agbw->xj;
  float_a *yj = agbw->yj;
  float_a *zj = agbw->zj;
#endif

  /* reset volumep array */
  memset(volumep,0,natoms*sizeof(float_a));

  /* set scaled volume coefficients */
  sign = 1.0;
  for(order=1;order<=AGBNP_MAX_OVERLAP_LEVEL;order++){
    volpcoeff[order-1] = sign/((float_a)order);
    sign *= -1.0;
  }

  /* computes self volumes */
  for(i=0;i<nheavyat;i++){
    iat = iheavyat[i];

    if(near_nl->nne[iat]<=0) continue;

    if(do_frozen){
      if(isvvfrozen[iat]) continue;
    }

    gatlist[0] = iat;
    gx[0][0] = x[iat];
    gx[0][1] = y[iat];
    gx[0][2] = z[iat];
    ga[0] = galpha[iat];
    gp[0] = gprefac[iat];
    gr[0] = r[iat];
    gparams[0].a = galpha[iat];
    gparams[0].p = gprefac[iat];
    gparams[0].c[0] = x[iat];
    gparams[0].c[1] = y[iat];
    gparams[0].c[2] = z[iat];
#ifdef USE_ALT
    gaw[0] = agb->occupancy[iat];
#endif

#ifdef USE_PBC
    if(agbnp_get_neigh_coords(agb, iat, near_nl, xj, yj, zj) != AGBNP_OK){
      fprintf(stderr, "agbnp_self_volumes(): error in agbnp_get_neigh_coords()\n");
      return AGBNP_ERR;
    } 
#endif

    order = 2;
    gnlist[order-1] = 0;

    while(order>1){

      j = gnlist[order-1];
#ifdef USE_PBC
      jat = near_nl->neighl[iat][j] % natoms;
#else
      jat = near_nl->neighl[iat][j];
#endif
      gatlist[order-1] = jat;

#ifdef USE_PBC      
      gx[order-1][0] = xj[j];
      gx[order-1][1] = yj[j];
      gx[order-1][2] = zj[j];
#else
      gx[order-1][0] = x[jat];
      gx[order-1][1] = y[jat];
      gx[order-1][2] = z[jat];
#endif
      ga[order-1] = galpha[jat];
      gp[order-1] = gprefac[jat];
      gr[order-1] = r[jat];

#ifdef USE_PBC
      gvol = ogauss_incremental_pbc(
#else
      gvol = ogauss_incremental(
#endif
		  order, gx, ga, gp, gr,
		  gparams[order-2].a, gparams[order-2].p, gparams[order-2].c,
		  &an, &pn, cn,
		  AGBNP_MIN_VOLA,AGBNP_MIN_VOLB,
		  NULL, NULL, NULL);

#ifdef USE_ALT
      if(agb->doalt){
	altw = altweight_incremental(agb, order-1, gatlist,
				     gaw[order-2],jat);
      }
#endif

      gparams[order-1].a = an;
      gparams[order-1].p = pn;
      gparams[order-1].c[0] = cn[0];
      gparams[order-1].c[1] = cn[1];
      gparams[order-1].c[2] = cn[2];
#ifdef USE_ALT
      gaw[order-1] = altw;
#endif

      /* store in cache overlap volume and positional derivatives */
      if(order == 2){
	if(!(do_frozen && isvfrozen[iat] && isvfrozen[jat])){ 
	  v2cache[iv2cache++] = altw*gvol;
	}
      }

      if(altw*gvol > FLT_MIN && order < AGBNP_MAX_OVERLAP_LEVEL){


	/* collect scaled volumes */
#ifdef USE_ALT
	u = altw*volpcoeff[order-1];
	v = altw*order*volpcoeff[order-1];
#else
	u = volpcoeff[order-1];
	v = order*volpcoeff[order-1];
#endif
	for(ii=0;ii<order;ii++){
	  kat = gatlist[ii];
	  volumep[kat] += u*gvol;
	}       

	
	/* increment overlap level */
	order += 1;
	/* initialize neighbor index at this level
	   (neighbor index at previous level + 1)
	*/
	gnlist[order-1] = gnlist[order-2]+1;
      }else{
	/* new neighbor at this level */
	gnlist[order-1] += 1;
      }

      /* 
	 if neighbor index is out of bounds there are no more neighbors at
	 this level. Decrement order and go to next neighbor at that order.
	 Do this recursively until an order with available neighbors is found.
      */
      while(order>1 && gnlist[order-1]>=near_nl->nne[iat]){
	order -= 1;
	gnlist[order-1] += 1;
      }

    }

  }	

  return AGBNP_OK;
}


#ifdef USE_PBC
 int agbnp_surface_areas_pbc(
#else
 int agbnp_surface_areas(
#endif
			 AGBNPdata *agb, CAVworkdata *cavw,
			 float_a *x, float_a *y, float_a *z){
  /* coordinate buffer for Gaussian overlap calculation */
  float_a gx[AGBNP_MAX_OVERLAP_LEVEL][3];
  /* radius buffer for  Gaussian overlap calculation */
  float_a gr[AGBNP_MAX_OVERLAP_LEVEL];
  /* Gaussian parameters buffers for overlap calculation */
  float_a ga[AGBNP_MAX_OVERLAP_LEVEL];
  float_a gp[AGBNP_MAX_OVERLAP_LEVEL];
  /* derivatives buffers for Gaussian overlaps */
  float_a gdr[AGBNP_MAX_OVERLAP_LEVEL][3];
  float_a gdR[AGBNP_MAX_OVERLAP_LEVEL];
  /* holds the atom indexes being overlapped */
  int gnlist[AGBNP_MAX_OVERLAP_LEVEL];
  int gatlist[AGBNP_MAX_OVERLAP_LEVEL];
  float_a gvol; /* gaussian overlap volume */
  int order; /* order of the overlap */
  /* holds overlap Gaussian parameters at each level */
  GParm gparams[AGBNP_MAX_OVERLAP_LEVEL];
  float_a an = 0.0, pn = 0.0, cn[3] = {0.0,0.0,0.0};
  /* coefficients for scaled volume */
  static float_a volpcoeff[AGBNP_MAX_OVERLAP_LEVEL];
  float_a sign;
  int i,j,iat,jat,kat,ii;
  float_a u,v, altw=1.0;
  float_a sumv1 = 0.0;

  int natoms = agb->natoms;
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  NeighList *near_nl = cavw->near_nl;
  float_a *r = agb->rcav;
  float_a *galpha = cavw->galpha;
  float_a *gprefac = cavw->gprefac;
  float_a *volumep = cavw->volumep;
  float_a *surf_area = cavw->surf_area;
  float_a *vols = cavw->vols;
#ifdef _OPENMP
  float_a *volumep_m = agb->cavw->volumep;
  float_a *surf_area_m = agb->cavw->surf_area;
#endif
#ifdef USE_PBC
  float_a gaw[AGBNP_MAX_OVERLAP_LEVEL];
  float_a *xj = cavw->xj;
  float_a *yj = cavw->yj;
  float_a *zj = cavw->zj;
#endif

  /* reset self volume and surface area accumulators */
  memset(volumep,0,natoms*sizeof(float_a));
  memset(surf_area,0,natoms*sizeof(float_a));

  /* set scaled volume coefficients */
  sign = 1.0;
  for(order=1;order<=AGBNP_MAX_OVERLAP_LEVEL;order++){
    volpcoeff[order-1] = sign/((float_a)order);
    sign *= -1.0;
  }

  /* computes n-body volumes */
  for(i=0;i<nheavyat;i++){
    iat = iheavyat[i];

    if(near_nl->nne[iat]<=0) continue;

    gatlist[0] = iat;
    gx[0][0] = x[iat];
    gx[0][1] = y[iat];
    gx[0][2] = z[iat];
    ga[0] = galpha[iat];
    gp[0] = gprefac[iat];
    gr[0] = r[iat];
    gparams[0].a = galpha[iat];
    gparams[0].p = gprefac[iat];
    gparams[0].c[0] = x[iat];
    gparams[0].c[1] = y[iat];
    gparams[0].c[2] = z[iat];
#ifdef USE_ALT
    gaw[0] = agb->occupancy[iat];
#endif

#ifdef USE_PBC
    if(agbnp_get_neigh_coords(agb, iat, near_nl, xj, yj, zj) != AGBNP_OK){
      fprintf(stderr, "agbnp_surface_areas(): error in agbnp_get_neigh_coords()\n");
      return AGBNP_ERR;
    } 
#endif

    order = 2;
    gnlist[order-1] = 0;

    while(order>1){

      j = gnlist[order-1];
#ifdef USE_PBC
      jat = near_nl->neighl[iat][j] % natoms;
#else
      jat = near_nl->neighl[iat][j];
#endif

      gatlist[order-1] = jat;
#ifdef USE_PBC
      gx[order-1][0] = xj[j];
      gx[order-1][1] = yj[j];
      gx[order-1][2] = zj[j];
#else
      gx[order-1][0] = x[jat];
      gx[order-1][1] = y[jat];
      gx[order-1][2] = z[jat];
#endif
      ga[order-1] = galpha[jat];
      gp[order-1] = gprefac[jat];
      gr[order-1] = r[jat];
     
#ifdef USE_PBC
      gvol = ogauss_incremental_pbc(
#else
      gvol = ogauss_incremental(
#endif
		order, gx, ga, gp, gr,
		gparams[order-2].a, gparams[order-2].p, gparams[order-2].c,
		&an, &pn, cn,
		AGBNP_MINC_VOLA, AGBNP_MINC_VOLB,
		gdr, gdR, NULL);
#ifdef USE_ALT
      if(agb->doalt){
	altw = altweight_incremental(agb, order-1, gatlist,
					   gaw[order-2],jat);
      }
#endif

      gparams[order-1].a = an;
      gparams[order-1].p = pn;
      gparams[order-1].c[0] = cn[0];
      gparams[order-1].c[1] = cn[1];
      gparams[order-1].c[2] = cn[2];
#ifdef USE_ALT
      gaw[order-1] = altw;
#endif

      if(altw*gvol > FLT_MIN && order < AGBNP_MAX_OVERLAP_LEVEL){

	/* collect scaled volumes, surface areas, 
	   and derivatives of cavity energy*/
#ifdef USE_ALT
	u = altw*volpcoeff[order-1];
	v = altw*order*volpcoeff[order-1];
#else
	u = volpcoeff[order-1];
	v = order*volpcoeff[order-1];
#endif
	for(ii=0;ii<order;ii++){
	  kat = gatlist[ii];
	  volumep[kat] += u*gvol;
	  surf_area[kat] += v*gdR[ii];
	}


	/* increment overlap level */
	order += 1;
	/* initialize neighbor index at this level
	   (neighbor index at previous level + 1)
	*/
	gnlist[order-1] = gnlist[order-2]+1;
      }else{
	/* new neighbor at this level */
	gnlist[order-1] += 1;
      }

      /* 
	 if neighbor index is out of bounds there are no more neighbors at
	 this level. Decrement order and go to next neighbor at that order.
	 Do this recursively until an order with available neighbors is found.
      */
      while(order>1 && gnlist[order-1]>=near_nl->nne[iat]){
	order -= 1;
	gnlist[order-1] += 1;
      }

    }

  }

#ifndef _OPENMP

  for(iat=0;iat<natoms;iat++){
    volumep[iat] += vols[iat];
  }
  for(i=0;i<nheavyat;i++){
    iat = iheavyat[i];
    surf_area[iat] +=  agb->occupancy[iat]*4.*pi*r[iat]*r[iat];
  }

#else

#pragma omp barrier
#pragma omp single
  {
    /* initialize master copy of self-volumes and surface areas */
    for(iat=0;iat<natoms;iat++){
      volumep_m[iat] = vols[iat];
    }
    memset(surf_area_m,0,natoms*sizeof(float_a));
    for(i=0;i<nheavyat;i++){
      iat = iheavyat[i];
      surf_area_m[iat] +=  agb->occupancy[iat]*4.*pi*r[iat]*r[iat];
    }
  }
#pragma omp barrier
#pragma omp critical
  {
    /* reduce surface areas and self-volumes */
    for(i=0;i<nheavyat;i++){
      iat = iheavyat[i];
      volumep_m[iat] += volumep[iat];
      surf_area_m[iat] += surf_area[iat];
    }
  }
#endif /* _OPENMP */

  return AGBNP_OK;
}


#ifdef USE_PBC
 int agbnp_cavity_ders_pbc(
#else
 int agbnp_cavity_ders(
#endif
		       AGBNPdata *agb, CAVworkdata *cavw,
		       float_a *x, float_a *y, float_a *z){
  /* coordinate buffer for Gaussian overlap calculation */
  float_a gx[AGBNP_MAX_OVERLAP_LEVEL][3];
  /* radius buffer for  Gaussian overlap calculation */
  float_a gr[AGBNP_MAX_OVERLAP_LEVEL];
  /* Gaussian parameters buffers for overlap calculation */
  float_a ga[AGBNP_MAX_OVERLAP_LEVEL];
  float_a gp[AGBNP_MAX_OVERLAP_LEVEL];
  /* derivatives buffers for Gaussian overlaps */
  float_a gdr[AGBNP_MAX_OVERLAP_LEVEL][3];
  float_a gdR[AGBNP_MAX_OVERLAP_LEVEL];
  float_a gd2rR[AGBNP_MAX_OVERLAP_LEVEL][AGBNP_MAX_OVERLAP_LEVEL][3];
  /* holds the atom indexes being overlapped */
  int gnlist[AGBNP_MAX_OVERLAP_LEVEL];
  int gatlist[AGBNP_MAX_OVERLAP_LEVEL];
  float_a gvol; /* gaussian overlap volume */
  int order; /* order of the overlap */
  /* holds overlap Gaussian parameters at each level */
  GParm gparams[AGBNP_MAX_OVERLAP_LEVEL];
  float_a an, pn, cn[3];
  /* coefficients for scaled volume */
  static float_a volpcoeff[AGBNP_MAX_OVERLAP_LEVEL];
  float_a sign;
  int i,j,iat,jat,kat,lat,ii,jj;
  float_a u,v, altw=1.0;
#ifdef USE_ALT
  float_a gaw[AGBNP_MAX_OVERLAP_LEVEL];
#endif
  int natoms = agb->natoms;
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  NeighList *near_nl = cavw->near_nl;
  float_a *r = agb->rcav;
  float_a *galpha = cavw->galpha;
  float_a *gprefac = cavw->gprefac;
  float_a *gammap = agb->cavw->gammap;
  float_a (*decav)[3] = cavw->decav_h;
#ifdef _OPENMP
  float_a (*decav_m)[3] = agb->cavw->decav_h;
#endif
#ifdef USE_PBC
  float_a *xj = cavw->xj;
  float_a *yj = cavw->yj;
  float_a *zj = cavw->zj;
  /* symmetry indexes */
  int ksym[AGBNP_MAX_OVERLAP_LEVEL];
  float_a (*rot)[3][3] = agb->rot;
#endif
  float_a ur[3];

  /* reset derivatives */
  memset(decav,0,3*natoms*sizeof(float_a));

  /* set scaled volume coefficients */
  sign = 1.0;
  for(order=1;order<=AGBNP_MAX_OVERLAP_LEVEL;order++){
    volpcoeff[order-1] = sign/((float_a)order);
    sign *= -1.0;
  }

  for(i=0;i<nheavyat;i++){
    iat = iheavyat[i];
    if(near_nl->nne[iat]<=0) continue;

    gatlist[0] = iat;
    gx[0][0] = x[iat];
    gx[0][1] = y[iat];
    gx[0][2] = z[iat];
    ga[0] = galpha[iat];
    gp[0] = gprefac[iat];
    gr[0] = r[iat];
    gparams[0].a = galpha[iat];
    gparams[0].p = gprefac[iat];
    gparams[0].c[0] = x[iat];
    gparams[0].c[1] = y[iat];
    gparams[0].c[2] = z[iat];
#ifdef USE_PBC
    ksym[0] = 0;
#endif
#ifdef USE_ALT
    gaw[0] = agb->occupancy[iat];
#endif

#ifdef USE_PBC
    if(agbnp_get_neigh_coords(agb, iat, near_nl, xj, yj, zj) != AGBNP_OK){
      fprintf(stderr, "agbnp_cavity_ders(): error in agbnp_get_neigh_coords()\n");
      return AGBNP_ERR;
    } 
#endif

    order = 2;
    gnlist[order-1] = 0;

    while(order>1){

      j = gnlist[order-1];
#ifdef USE_PBC
      jat = near_nl->neighl[iat][j] % natoms;
#else
      jat = near_nl->neighl[iat][j];
#endif

      gatlist[order-1] = jat;
#ifdef USE_PBC
      gx[order-1][0] = xj[j];
      gx[order-1][1] = yj[j];
      gx[order-1][2] = zj[j];
#else
      gx[order-1][0] = x[jat];
      gx[order-1][1] = y[jat];
      gx[order-1][2] = z[jat];
#endif
      ga[order-1] = galpha[jat];
      gp[order-1] = gprefac[jat];
      gr[order-1] = r[jat];
#ifdef USE_PBC
      ksym[order-1] = near_nl->neighl[iat][j]/natoms;
#endif
      
#ifdef USE_PBC
      gvol = ogauss_incremental_pbc(
#else
      gvol = ogauss_incremental(
#endif
		 order, gx, ga, gp, gr,
	         gparams[order-2].a, gparams[order-2].p, gparams[order-2].c,
		 &an, &pn, cn,
		 AGBNP_MINC_VOLA, AGBNP_MINC_VOLB,
		 gdr, gdR, gd2rR);

#ifdef USE_ALT
      if(agb->doalt){
	altw = altweight_incremental(agb, order-1, gatlist,
				     gaw[order-2],jat);
      }
#endif

      gparams[order-1].a = an;
      gparams[order-1].p = pn;
      gparams[order-1].c[0] = cn[0];
      gparams[order-1].c[1] = cn[1];
      gparams[order-1].c[2] = cn[2];
#ifdef USE_ALT
      gaw[order-1] = altw;
#endif

      if(altw*gvol > FLT_MIN && order < AGBNP_MAX_OVERLAP_LEVEL){

	/* derivatives of cavity energy*/
#ifdef USE_ALT
	u = altw*volpcoeff[order-1];
	v = altw*order*volpcoeff[order-1];
#else
	u = volpcoeff[order-1];
	v = order*volpcoeff[order-1];
#endif
	for(ii=0;ii<order;ii++){
	  kat = gatlist[ii];
	  for(jj=0;jj<order;jj++){
	    lat = gatlist[jj];
#ifdef USE_PBC
	    if(agb->docryst){
	      rtvec(ur,rot[ksym[ii]],gd2rR[ii][jj]);
	    }else{
	      ur[0] = gd2rR[ii][jj][0];
	      ur[1] = gd2rR[ii][jj][1];
	      ur[2] = gd2rR[ii][jj][2];
	    }
#else
	    ur[0] = gd2rR[ii][jj][0];
	    ur[1] = gd2rR[ii][jj][1];
	    ur[2] = gd2rR[ii][jj][2];
#endif
	    decav[kat][0] += v*gammap[lat]*ur[0];
	    decav[kat][1] += v*gammap[lat]*ur[1];
	    decav[kat][2] += v*gammap[lat]*ur[2];
	  }
	}      

	/* increment overlap level */
	order += 1;
	/* initialize neighbor index at this level
	   (neighbor index at previous level + 1)
	*/
	gnlist[order-1] = gnlist[order-2]+1;
      }else{
	/* new neighbor at this level */
	gnlist[order-1] += 1;
      }

      /* 
	 if neighbor index is out of bounds there are no more neighbors at
	 this level. Decrement order and go to next neighbor at that order.
	 Do this recursively until an order with available neighbors is found.
      */
      while(order>1 && gnlist[order-1]>=near_nl->nne[iat]){
	order -= 1;
	gnlist[order-1] += 1;
      }

    }
  }

#ifdef _OPENMP
#pragma omp single
  memset(decav_m,0,3*natoms*sizeof(float_a));
#pragma omp barrier
#pragma omp critical
  /* reduction of derivatives */
  {
    for(iat=0;iat<natoms;iat++){
      for(i=0;i<3;i++){
	decav_m[iat][i] += decav[iat][i];
      }
    }
  }
#pragma omp barrier
#endif

  return AGBNP_OK;
}


/* GB and vdw derivatives contribution at constant self volumes */
#ifdef USE_PBC
 int agbnp_gb_ders_constvp_pbc(
#else
 int agbnp_gb_ders_constvp(
#endif

			       AGBNPdata *agb, AGBworkdata *agbw_h,
			       float_a *x, float_a *y, float_a *z,
			       int init_frozen){

  float_a dielectric_factor = 
    -0.5*(1./agb->dielectric_in - 1./agb->dielectric_out);
  int iq4cache=0, iv2cache = 0; /* cache counters */
  int i, iat, j, jat;
  float_a dx, dy, dz, d2, d, volume2, q, dr4, spiat, htij, utij, spjat;
  float_a fourpi1 = 1./(4.*pi);
  
  int natoms = agb->natoms;
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  int *isheavy = agbw_h->isheavy;
  int nhydrogen = agb->nhydrogen;
  int *ihydrogen = agb->ihydrogen;
  NeighList *near_nl = agbw_h->near_nl;
  NeighList *far_nl = agbw_h->far_nl;
  float_a *sp = agbw_h->sp;
  float_a *vols = agbw_h->vols;
  float_a *q2ab = agbw_h->q2ab;
  float_a *abrw = agbw_h->abrw;
  float_a (*dgbdr)[3] = agbw_h->dgbdr_h;
  float_a (*dvwdr)[3] = agbw_h->dvwdr_h;
  float *v2cache = agbw_h->v2cache;
  float *q4cache = agbw_h->q4cache;

  int do_frozen = agb->do_frozen;
  int *isvfrozen = agb->agbw->isvfrozen;
  float_a (*dgbdr_c)[3] = agb->agbw->dgbdr_c;
  float_a (*dvwdr_c)[3] = agb->agbw->dvwdr_c;
#ifdef USE_PBC
  float_a *xj = agbw_h->xj;
  float_a *yj = agbw_h->yj;
  float_a *zj = agbw_h->zj;
  float_a (*rot)[3][3] = agb->rot;
  int jsym;
#endif
  float_a u[3], ur[3];
  float_a w[3], wr[3];

  /* loop over near heavy-heavy interactions */
  for(i=0;i<nheavyat;i++){
    iat = iheavyat[i];
#ifdef USE_PBC
    if(agbnp_get_neigh_coords(agb, iat, near_nl, xj, yj, zj) != AGBNP_OK){
      fprintf(stderr, "agbnp_gb_ders_constvp(): error in agbnp_get_neigh_coords()\n");
      return AGBNP_ERR;
    } 
#endif
    for(j=0;j<near_nl->nne[iat];j++){
      jat = near_nl->neighl[iat][j];
#ifdef USE_PBC
      jsym = jat/natoms;
      jat = jat % natoms;
#endif
      if(do_frozen){
	/* skip vfrozen pairs */
	if(isvfrozen[iat] && isvfrozen[jat]) continue;
      }
#ifdef USE_PBC
      dx = xj[j] - x[iat];
      dy = yj[j] - y[iat];
      dz = zj[j] - z[iat];
#else
      dx = x[jat] - x[iat];
      dy = y[jat] - y[iat];
      dz = z[jat] - z[iat];
#endif
      d2 = dx*dx + dy*dy + dz*dz;
      d = mysqrt(d2);
      volume2 = v2cache[iv2cache++];
      q = q4cache[iq4cache++];
      dr4 = q4cache[iq4cache++];
      spiat = sp[iat]+0.5*volume2/vols[iat];
      htij = q2ab[jat]*dr4*spiat;
      utij = abrw[jat]*dr4*spiat;
      q = q4cache[iq4cache++];
      dr4 = q4cache[iq4cache++];
      spjat = sp[jat]+0.5*volume2/vols[jat];
      htij += q2ab[iat]*dr4*spjat;
      utij += abrw[iat]*dr4*spjat;

      htij = fourpi1*dielectric_factor*htij/d;	
      u[0] = ur[0] = htij*dx;
      u[1] = ur[1] = htij*dy;
      u[2] = ur[2] = htij*dz;
#ifdef USE_PBC
      if(agb->docryst){
        /* rotate gradient */
        rtvec(ur,rot[jsym],u);
      }
#endif
      dgbdr[iat][0] +=  u[0]; 
      dgbdr[iat][1] +=  u[1]; 
      dgbdr[iat][2] +=  u[2]; 
      dgbdr[jat][0] -=  ur[0];
      dgbdr[jat][1] -=  ur[1];
      dgbdr[jat][2] -=  ur[2];
      
      utij = fourpi1*utij/d;
      w[0] = wr[0] = utij*dx;
      w[1] = wr[1] = utij*dy;
      w[2] = wr[2] = utij*dz;
#ifdef USE_PBC
      if(agb->docryst){
        /* rotate gradient */
        rtvec(wr,rot[jsym],w);
      }
#endif
      dvwdr[iat][0] +=  w[0]; 
      dvwdr[iat][1] +=  w[1]; 
      dvwdr[iat][2] +=  w[2]; 
      dvwdr[jat][0] -=  wr[0];
      dvwdr[jat][1] -=  wr[1];
      dvwdr[jat][2] -=  wr[2];

      if(init_frozen){
	if(isvfrozen[iat] && isvfrozen[jat]) {
#pragma omp critical
	  {
	    dgbdr_c[iat][0] +=  u[0];   
	    dgbdr_c[iat][1] +=  u[1];  
	    dgbdr_c[iat][2] +=  u[2];  
	    dgbdr_c[jat][0] -=  ur[0]; 
	    dgbdr_c[jat][1] -=  ur[1]; 
	    dgbdr_c[jat][2] -=  ur[2]; 
	    
	    dvwdr_c[iat][0] +=  w[0];   
	    dvwdr_c[iat][1] +=  w[1];  
	    dvwdr_c[iat][2] +=  w[2];  
	    dvwdr_c[jat][0] -=  wr[0]; 
	    dvwdr_c[jat][1] -=  wr[1]; 
	    dvwdr_c[jat][2] -=  wr[2]; 
	  }
	}
      }


    }
  }
  /* loop over far heavy-heavy interactions */
  for(i=0;i<nheavyat;i++){
    iat = iheavyat[i];
#ifdef USE_PBC
    if(agbnp_get_neigh_coords(agb, iat, far_nl, xj, yj, zj) != AGBNP_OK){
      fprintf(stderr, "agbnp_gb_ders_constvp(): error in agbnp_get_neigh_coords()\n");
      return AGBNP_ERR;
    } 
#endif
    for(j=0;j<far_nl->nne[iat];j++){
      jat = far_nl->neighl[iat][j];
#ifdef USE_PBC
      jsym = jat/natoms;
      jat = jat % natoms;
#endif
      if(do_frozen){
	if(isvfrozen[iat] && isvfrozen[jat]) continue;
      }
#ifdef USE_PBC
      dx = xj[j] - x[iat];
      dy = yj[j] - y[iat];
      dz = zj[j] - z[iat];
#else
      dx = x[jat] - x[iat];
      dy = y[jat] - y[iat];
      dz = z[jat] - z[iat];
#endif
      d2 = dx*dx + dy*dy + dz*dz;
      d = mysqrt(d2);
      /* get from cache */
      if(isheavy[jat]){
	iq4cache += 1;
	dr4 = q4cache[iq4cache++];
	htij = q2ab[jat]*dr4*sp[iat];
	utij = abrw[jat]*dr4*sp[iat];
	iq4cache += 1;
	dr4 = q4cache[iq4cache++];
	htij += q2ab[iat]*dr4*sp[jat];
	utij += abrw[iat]*dr4*sp[jat];
      }else{
	iq4cache += 1;
	dr4 = q4cache[iq4cache++];
	htij = q2ab[jat]*dr4*sp[iat];
	utij = abrw[jat]*dr4*sp[iat];
      }


      htij = fourpi1*dielectric_factor*htij/d;	
      u[0] = ur[0] = htij*dx;
      u[1] = ur[1] = htij*dy;
      u[2] = ur[2] = htij*dz;
#ifdef USE_PBC
      if(agb->docryst){
        /* rotate gradient */
        rtvec(ur,rot[jsym],u);
      }
#endif
      dgbdr[iat][0] +=  u[0]; 
      dgbdr[iat][1] +=  u[1]; 
      dgbdr[iat][2] +=  u[2]; 
      dgbdr[jat][0] -=  ur[0];
      dgbdr[jat][1] -=  ur[1];
      dgbdr[jat][2] -=  ur[2];
      
      utij = fourpi1*utij/d;
      w[0] = wr[0] = utij*dx;
      w[1] = wr[1] = utij*dy;
      w[2] = wr[2] = utij*dz;
#ifdef USE_PBC
      if(agb->docryst){
        /* rotate gradient */
        rtvec(wr,rot[jsym],w);
      }
#endif
      dvwdr[iat][0] +=  w[0]; 
      dvwdr[iat][1] +=  w[1]; 
      dvwdr[iat][2] +=  w[2]; 
      dvwdr[jat][0] -=  wr[0];
      dvwdr[jat][1] -=  wr[1];
      dvwdr[jat][2] -=  wr[2];

      if(init_frozen){
	if(isvfrozen[iat] && isvfrozen[jat]) {
#pragma omp critical
	  {
	    dgbdr_c[iat][0] +=  u[0];   
	    dgbdr_c[iat][1] +=  u[1];  
	    dgbdr_c[iat][2] +=  u[2];  
	    dgbdr_c[jat][0] -=  ur[0]; 
	    dgbdr_c[jat][1] -=  ur[1]; 
	    dgbdr_c[jat][2] -=  ur[2]; 
	    
	    dvwdr_c[iat][0] +=  w[0];   
	    dvwdr_c[iat][1] +=  w[1];  
	    dvwdr_c[iat][2] +=  w[2];  
	    dvwdr_c[jat][0] -=  wr[0]; 
	    dvwdr_c[jat][1] -=  wr[1]; 
	    dvwdr_c[jat][2] -=  wr[2]; 
	  }
	}
      }

    }
  }
  /* loop for hydrogen-heavy interactions */
  for(i=0;i<nhydrogen;i++){
    iat = ihydrogen[i];
#ifdef USE_PBC
    if(agbnp_get_neigh_coords(agb, iat, far_nl, xj, yj, zj) != AGBNP_OK){
      fprintf(stderr, "agbnp_gb_ders_constvp(): error in agbnp_get_neigh_coords()\n");
      return AGBNP_ERR;
    } 
#endif
    for(j=0;j< far_nl->nne[iat];j++){
      jat = far_nl->neighl[iat][j];
#ifdef USE_PBC
      jsym = jat/natoms;
      jat = jat % natoms;
#endif
      if(do_frozen){
	if(isvfrozen[jat] && isvfrozen[iat]) continue;
      }
      if(isheavy[jat]){
#ifdef USE_PBC
	dx = xj[j] - x[iat];
	dy = yj[j] - y[iat];
	dz = zj[j] - z[iat];
#else
	dx = x[jat] - x[iat];
	dy = y[jat] - y[iat];
	dz = z[jat] - z[iat];
#endif
	d2 = dx*dx + dy*dy + dz*dz;
	d = mysqrt(d2);
	iq4cache += 1;
	dr4 = q4cache[iq4cache++];
	htij = q2ab[iat]*dr4*sp[jat];
	utij = abrw[iat]*dr4*sp[jat];

	htij = fourpi1*dielectric_factor*htij/d;	
	u[0] = ur[0] = htij*dx;
	u[1] = ur[1] = htij*dy;
	u[2] = ur[2] = htij*dz;
#ifdef USE_PBC
	if(agb->docryst){
	  /* rotate gradient */
	  rtvec(ur,rot[jsym],u);
	}
#endif
	dgbdr[iat][0] +=  u[0]; 
	dgbdr[iat][1] +=  u[1]; 
	dgbdr[iat][2] +=  u[2]; 
	dgbdr[jat][0] -=  ur[0];
	dgbdr[jat][1] -=  ur[1];
	dgbdr[jat][2] -=  ur[2];
	
	utij = fourpi1*utij/d;
	w[0] = wr[0] = utij*dx;
	w[1] = wr[1] = utij*dy;
	w[2] = wr[2] = utij*dz;
#ifdef USE_PBC
	if(agb->docryst){
	  /* rotate gradient */
	  rtvec(wr,rot[jsym],w);
	}
#endif
	dvwdr[iat][0] +=  w[0]; 
	dvwdr[iat][1] +=  w[1]; 
	dvwdr[iat][2] +=  w[2]; 
	dvwdr[jat][0] -=  wr[0];
	dvwdr[jat][1] -=  wr[1];
	dvwdr[jat][2] -=  wr[2];
	
	if(init_frozen){
	  if(isvfrozen[iat] && isvfrozen[jat]) {
#pragma omp critical
	    {
	      dgbdr_c[iat][0] +=  u[0];   
	      dgbdr_c[iat][1] +=  u[1];  
	      dgbdr_c[iat][2] +=  u[2];  
	      dgbdr_c[jat][0] -=  ur[0]; 
	      dgbdr_c[jat][1] -=  ur[1]; 
	      dgbdr_c[jat][2] -=  ur[2]; 
	      
	      dvwdr_c[iat][0] +=  w[0];   
	      dvwdr_c[iat][1] +=  w[1];  
	      dvwdr_c[iat][2] +=  w[2];  
	      dvwdr_c[jat][0] -=  wr[0]; 
	      dvwdr_c[jat][1] -=  wr[1]; 
	      dvwdr_c[jat][2] -=  wr[2]; 
	    }
	  }
	}


      }
    }
  }

  return AGBNP_OK;
}

/* derivatives due to changes in self volumes */
#ifdef USE_PBC
 int agbnp_der_vp_pbc(
#else
 int agbnp_der_vp(
#endif
		  AGBNPdata *agb, AGBworkdata *agbw_h,
		  float_a *x, float_a *y, float_a *z,
		  int init_frozen){

  int natoms = agb->natoms;
  int nheavyat = agb->nheavyat;
  int *iheavyat = agb->iheavyat;
  NeighList *near_nl = agbw_h->near_nl;
  float_a *r = agb->r;
  float_a *galpha = agbw_h->galpha;
  float_a *gprefac = agbw_h->gprefac;
  float_a *vols = agbw_h->vols;
  float_a *q2ab = agbw_h->q2ab;
  float_a *abrw = agbw_h->abrw;
  float_a *deru = agb->agbw->deru;
  float_a *derv = agb->agbw->derv;
  float_a (*dgbdr)[3] = agbw_h->dgbdr_h;
  float_a (*dvwdr)[3] = agbw_h->dvwdr_h;
  float *q4cache = agbw_h->q4cache;

  int iq4cache = 0;
  float_a dielectric_factor = 
    -0.5*(1./agb->dielectric_in - 1./agb->dielectric_out);
  float_a fourpi1 = 1./(4.*pi);
  /* coordinate buffer for Gaussian overlap calculation */
  float_a gx[AGBNP_MAX_OVERLAP_LEVEL][3];
  /* radius buffer for  Gaussian overlap calculation */
  float_a gr[AGBNP_MAX_OVERLAP_LEVEL];
  /* Gaussian parameters buffers for overlap calculation */
  float_a ga[AGBNP_MAX_OVERLAP_LEVEL];
  float_a gp[AGBNP_MAX_OVERLAP_LEVEL];
  /* holds the atom indexes being overlapped */
  int gnlist[AGBNP_MAX_OVERLAP_LEVEL];
  int gatlist[AGBNP_MAX_OVERLAP_LEVEL];
  float_a gvol; /* gaussian overlap volume */
  int order; /* order of the overlap */
  /* holds overlap Gaussian parameters at each level */
  GParm gparams[AGBNP_MAX_OVERLAP_LEVEL];
  float_a an, pn, cn[3];
  /* coefficients for scaled volume */
  float_a volpcoeff[AGBNP_MAX_OVERLAP_LEVEL];
  float_a sign;
  /* derivatives buffers for Gaussian overlaps */
  float_a gdr[AGBNP_MAX_OVERLAP_LEVEL][3];
  float_a gdR[AGBNP_MAX_OVERLAP_LEVEL];
#ifdef USE_PBC
  /* symmetry indexes */
  int ksym[AGBNP_MAX_OVERLAP_LEVEL];
  float_a (*rot)[3][3] = agb->rot;
#endif
  float_a ur[3], altw=1.0;
#ifdef USE_ALT
  float_a gaw[AGBNP_MAX_OVERLAP_LEVEL]; 
#endif

  int i, j, iat, jat, kat, ii;
  float_a u,v,deruij,deruji,dervij,dervji,q;
  int do_frozen = agb->do_frozen;
  int *isvvfrozen = agb->agbw->isvvfrozen;
  float_a (*dgbdr_c)[3] = agb->agbw->dgbdr_c;
  float_a (*dvwdr_c)[3] = agb->agbw->dvwdr_c;
#ifdef USE_PBC
  float_a *xj = agbw_h->xj;
  float_a *yj = agbw_h->yj;
  float_a *zj = agbw_h->zj;
#endif

  /* set scaled volume coefficients */
  sign = 1.0;
  for(order=1;order<=AGBNP_MAX_OVERLAP_LEVEL;order++){
    volpcoeff[order-1] = sign/((float_a)order);
    sign *= -1.0;
  }

  for(i=0;i<nheavyat;i++){
    iat = iheavyat[i];

    if(near_nl->nne[iat]<=0) continue;

    if(do_frozen){
      if(isvvfrozen[iat]) continue;
    }

    gatlist[0] = iat;
    gx[0][0] = x[iat];
    gx[0][1] = y[iat];
    gx[0][2] = z[iat];
    ga[0] = galpha[iat];
    gp[0] = gprefac[iat];
    gr[0] = r[iat];
    gparams[0].a = galpha[iat];
    gparams[0].p = gprefac[iat];
    gparams[0].c[0] = x[iat];
    gparams[0].c[1] = y[iat];
    gparams[0].c[2] = z[iat];
#ifdef USE_PBC
    ksym[0] = 0;
#endif
#ifdef USE_ALT
    gaw[0] = agb->occupancy[iat];
#endif

#ifdef USE_PBC
    if(agbnp_get_neigh_coords(agb, iat, near_nl, xj, yj, zj) != AGBNP_OK){
      fprintf(stderr, "agbnp_der_vp(): error in agbnp_get_neigh_coords()\n");
      return AGBNP_ERR;
    }
#endif

    order = 2;
    gnlist[order-1] = 0;

    while(order>1){

      j = gnlist[order-1];
      jat = near_nl->neighl[iat][j];
#ifdef USE_PBC
      ksym[order-1] = jat/natoms;
      jat = jat % natoms;
#endif

      gatlist[order-1] = jat;
#ifdef USE_PBC
      gx[order-1][0] = xj[j];
      gx[order-1][1] = yj[j];
      gx[order-1][2] = zj[j];
#else
      gx[order-1][0] = x[jat];
      gx[order-1][1] = y[jat];
      gx[order-1][2] = z[jat];
#endif
      ga[order-1] = galpha[jat];
      gp[order-1] = gprefac[jat];
      gr[order-1] = r[jat];
      
#ifdef USE_PBC
      gvol = ogauss_incremental_pbc(
#else
      gvol = ogauss_incremental(
#endif
		 order, gx, ga, gp, gr,
	         gparams[order-2].a, gparams[order-2].p, gparams[order-2].c,
		 &an, &pn, cn, 
		 AGBNP_MIN_VOLA, AGBNP_MIN_VOLB,
		 gdr, gdR, NULL);

#ifdef USE_ALT
      if(agb->doalt){
	altw = altweight_incremental(agb, order-1, gatlist,
				     gaw[order-2],jat);
      }
#endif

      gparams[order-1].a = an;
      gparams[order-1].p = pn;
      gparams[order-1].c[0] = cn[0];
      gparams[order-1].c[1] = cn[1];
      gparams[order-1].c[2] = cn[2];
#ifdef USE_ALT
      gaw[order-1] = altw;
#endif

      if(order == 2){
	q = q4cache[iq4cache++];
	iq4cache += 1;
	deruij = q2ab[jat]*q/vols[iat];
	dervij = abrw[jat]*q/vols[iat];
	q = q4cache[iq4cache++];
	iq4cache += 1;
	deruji = q2ab[iat]*q/vols[jat];
	dervji = abrw[iat]*q/vols[jat];
	u = -(deruij+deruji);
	v = -(dervij+dervji);
      }else{    
	u = 0.0; 
	v = 0.0;
      }

      if(altw*gvol > FLT_MIN && order < AGBNP_MAX_OVERLAP_LEVEL){

      /* contribution to GB energy derivatives due to the change in scaled
	 volumes */
      for(ii=0;ii<order;ii++){
	kat = gatlist[ii];
	u += deru[kat];
	v += derv[kat];
      }

#ifdef USE_ALT
      u *= -altw*dielectric_factor*fourpi1*volpcoeff[order-1];
      v *= -altw*fourpi1*volpcoeff[order-1];
#else
      u *= -dielectric_factor*fourpi1*volpcoeff[order-1];
      v *= -fourpi1*volpcoeff[order-1];
#endif
      for(ii=0;ii<order;ii++){
	kat = gatlist[ii];
#ifdef USE_PBC
	if(agb->docryst){
	  rtvec(ur,rot[ksym[ii]],gdr[ii]);
	}else{
	  ur[0] = gdr[ii][0];
	  ur[1] = gdr[ii][1];
	  ur[2] = gdr[ii][2];
	}
#else
	ur[0] = gdr[ii][0];
	ur[1] = gdr[ii][1];
	ur[2] = gdr[ii][2];
#endif

	dgbdr[kat][0] += u*ur[0];
	dgbdr[kat][1] += u*ur[1];
	dgbdr[kat][2] += u*ur[2];
	
	dvwdr[kat][0] += v*ur[0];
	dvwdr[kat][1] += v*ur[1];
	dvwdr[kat][2] += v*ur[2];

	if(init_frozen){
	  if(isvvfrozen[iat]){
#pragma omp critical
	    {
	      dgbdr_c[kat][0] += u*ur[0];
	      dgbdr_c[kat][1] += u*ur[1];
	      dgbdr_c[kat][2] += u*ur[2];
	      
	      dvwdr_c[kat][0] += v*ur[0];
	      dvwdr_c[kat][1] += v*ur[1];
	      dvwdr_c[kat][2] += v*ur[2];
	    }
	  }

	}

      }



	/* increment overlap level */
	order += 1;
	/* initialize neighbor index at this level
	   (neighbor index at previous level + 1)
	*/
	gnlist[order-1] = gnlist[order-2]+1;
      }else{
	/* new neighbor at this level */
	gnlist[order-1] += 1;
      }

      /* 
	 if neighbor index is out of bounds there are no more neighbors at
	 this level. Decrement order and go to next neighbor at that order.
	 Do this recursively until an order with available neighbors is found.
      */
      while(order>1 && gnlist[order-1]>=near_nl->nne[iat]){
	order -= 1;
	gnlist[order-1] += 1;
      }

    }

  }

  return AGBNP_OK;
}

#ifdef USE_PBC
 int agbnp_getatomneighbors_pbc(
#else
 int agbnp_getatomneighbors(
#endif

			    AGBNPdata *agb, int iat, int *nbnum, 
			    int bsize, int *nbiat
#ifdef USE_PBC
			    ,float_a *xj, float_a *yj, float_a *zj,
			    void **nbdata){
#else
                            ){
#endif

  int i,j, jat,jatom;
  int excl_nbnum = 0, excl_bsize;
  NeighList *neigh_list = agb->neigh_list;
  NeighList *excl_neigh_list = agb->excl_neigh_list;
  int *excl_nbiat;
  int *ext2int;
  int *int2ext;
#ifdef USE_PBC
  float_a *xa = agb->x;
  float_a *ya = agb->y;
  float_a *za = agb->z;
  float_i *xs = agb->xs;
  float_i *ys = agb->ys;
  float_i *zs = agb->zs;
  float_a *xje, *yje, *zje;
  NeighVector **data_index;
  NeighVector *pbc_trans, *pbc_trans_source;
#endif

  *nbnum = 0;

  if(neigh_list){
    /* retrieve neighbors from Verlet neighbor lists */

    /* Verlet neighbor list uses internal atom indeces starting from 0 as 
       AGBNP. The Verlet neighbor list however stores the external atom 
       indeces of the neighboring atoms. The ext2int[] array maps 
       from external to internal atom indeces. */
    if(iat >= neigh_list->natoms){
      fprintf(stderr, "agbnp_getatomneighbors(): atom %d is not found in Verlet neighbor list\n",iat);
      return AGBNP_ERR;
    }
    if(neigh_list->nne[iat] > bsize){
      fprintf(stderr, "agbnp_getatomneighbors(): number of neighbors (%d) exceeds provided buffer size (%d).\n",neigh_list->nne[iat], bsize);
      return AGBNP_ERR;
    }
    *nbnum = neigh_list->nne[iat];
    memcpy(nbiat,neigh_list->neighl[iat],*nbnum*sizeof(int));
    /* neigh_list holds the external atom indeces of neighbors as
       defined by the application that created it.  Therefore need to
       convert external atom indeces to internal atom indeces. */
    if(neigh_list->idx_remap){
      ext2int = neigh_list->ext2int;
      for(i=0;i<*nbnum;i++){
        nbiat[i] = ext2int[nbiat[i]];
      }
    }

#ifdef USE_PBC
    if(agb->dopbc){
      
      if(!(nbdata && xj && yj && zj)){
	fprintf(stderr, "agbnp_getatomneighbors(): PBC is on but NULL coordinate and PBC translation vectors buffers provided.\n");
	return AGBNP_ERR;
      }
      if(!(neigh_list->data || neigh_list->pbc)){
	fprintf(stderr, "agbnp_getatomneighbors(): PBC is on but cannot find PBC translation vectors in neighbor list.\n");
	return AGBNP_ERR;
      }

      /* return pointers to PBC translation vectors */
      if(neigh_list->data){
	memcpy(nbdata,neigh_list->data_index[iat],*nbnum*sizeof(void *));
      }else{
	for(j=0;j<*nbnum;j++){
	  nbdata[j] = &(neigh_list->pbc_trans[iat][j]);
	}
      }
      
      if(neigh_list->data){
	data_index = (NeighVector **)neigh_list->data_index[iat];
      }else{
	pbc_trans_source = (NeighVector *)neigh_list->pbc_trans[iat];
      }

      /* return PBC-adjusted coordinates of neighbors */
      if(agb->nsym == 1){

	/* apply PBC's but use standard coordinates */
	for(j=0;j<*nbnum;j++){
	  if(neigh_list->data){
	    pbc_trans = data_index[j];
	  }else{
	    pbc_trans = &(pbc_trans_source[j]);
	  }
	  jat = nbiat[j];
	  xj[j] = xa[jat] - pbc_trans->x;
	  yj[j] = ya[jat] - pbc_trans->y;
	  zj[j] = za[jat] - pbc_trans->z;
	}
	
      }else{
	
	/* crystal PBC's, use external coordinate buffers with external 
	   atom indexes*/
	int2ext = neigh_list->int2ext;
	for(j=0;j<*nbnum;j++){
	  if(neigh_list->data){
            pbc_trans = data_index[j];
          }else{
	    /* fprintf(stderr,"&pbc_trans = %d %d %x\n",iat,j,&(pbc_trans_source[j]));
	       fflush(stderr); */
            pbc_trans = &(pbc_trans_source[j]);
          }
	  jat = nbiat[j];
	  jatom = neigh_list->idx_remap ? int2ext[jat] : jat;
	  /*	  fprintf(stderr,"agbnp-pbc: %d %d %lf %lf %lf\n",iat,jat,pbc_trans->x,pbc_trans->y,pbc_trans->z);
		  fflush(stderr); */
	  xj[j] = xs[jatom] - pbc_trans->x;
	  yj[j] = ys[jatom] - pbc_trans->y;
	  zj[j] = zs[jatom] - pbc_trans->z;
	}
      }
      
    }else{
      
      /* not doing PBC's, just copy coordinates */
      for(j=0;j<*nbnum;j++){
	jat = nbiat[j];
	xj[j] = xa[jat];
	yj[j] = ya[jat];
	zj[j] = za[jat];
      }	

    }
#endif

    if(neigh_list && excl_neigh_list){
      /* retrieve additional neighbors from Verlet list of excluded 
	 neighbors */
      
      excl_nbiat = &(nbiat[*nbnum]);
      excl_bsize = bsize - *nbnum;
      
      if(iat >= excl_neigh_list->natoms){
	fprintf(stderr, "agbnp_getatomneighbors(): atom %d is not found in exclided atoms Verlet neighbor list\n",iat);
	return AGBNP_ERR;
      }
      if(excl_neigh_list->nne[iat] > excl_bsize){
	fprintf(stderr, "agbnp_getatomneighbors(): number of excluded neighbors (%d) exceeds provided buffer size (%d).\n",excl_neigh_list->nne[iat], excl_bsize);
	return AGBNP_ERR;
      }
      excl_nbnum = excl_neigh_list->nne[iat];
      memcpy(excl_nbiat,excl_neigh_list->neighl[iat],
	     excl_nbnum*sizeof(int));
      if(neigh_list->idx_remap){
	ext2int = neigh_list->ext2int;
	for(i=0;i<excl_nbnum;i++){
	  excl_nbiat[i] = ext2int[excl_nbiat[i]];
	}
      }

#ifdef USE_PBC
      xje = &(xj[*nbnum]);
      yje = &(yj[*nbnum]);
      zje = &(zj[*nbnum]);
      if(agb->dopbc){

	if(!(nbdata && xj && yj && zj)){
	  fprintf(stderr, "agbnp_getatomneighbors(): PBC is on but NULL coordinate and PBC translation vectors buffers provided.\n");
	  return AGBNP_ERR;
	}
	if(!(excl_neigh_list->data || excl_neigh_list->pbc)){
	  fprintf(stderr, "agbnp_getatomneighbors(): PBC is on but cannot find PBC translation vectors in excluded neighbor list.\n");
	  return AGBNP_ERR;
	}

	/* return pointers to PBC translation vectors */
	if(excl_neigh_list->data){
	  memcpy(&(nbdata[*nbnum]),excl_neigh_list->data_index[iat],excl_nbnum*sizeof(void *));
	}else{
	  for(j=0;j<excl_nbnum;j++){
	    nbdata[j+(*nbnum)] = &(excl_neigh_list->pbc_trans[iat][j]);
	  }
	}
	
	if(excl_neigh_list->data){
	  data_index = (NeighVector **)excl_neigh_list->data_index[iat];
	}else{
	  pbc_trans_source = (NeighVector *)excl_neigh_list->pbc_trans[iat];
	}

	/* return PBC-adjusted coordinates of neighbors */
	if(agb->nsym == 1){
	  
	  /* apply PBC's but use standard coordinates */
	  for(j=0;j<excl_nbnum;j++){
	    if(excl_neigh_list->data){
	      pbc_trans = data_index[j];
	    }else{
	      pbc_trans = &(pbc_trans_source[j]);
	    }
	    jat = excl_nbiat[j];
	    xje[j] = xa[jat] - pbc_trans->x;
	    yje[j] = ya[jat] - pbc_trans->y;
	    zje[j] = za[jat] - pbc_trans->z;
	  }
	  
	}else{
	  
	  /* crystal PBC's, use external coordinate buffers with external 
	     atom indexes*/
	  int2ext = neigh_list->int2ext;
	  for(j=0;j<excl_nbnum;j++){
	    if(excl_neigh_list->data){
	      pbc_trans = data_index[j];
	    }else{
	      pbc_trans = &(pbc_trans_source[j]);
	    }
	    jat = excl_nbiat[j];
	    jatom = neigh_list->idx_remap ? int2ext[jat] : jat;
	    xje[j] = xs[jatom] - pbc_trans->x;
	    yje[j] = ys[jatom] - pbc_trans->y;
	    zje[j] = zs[jatom] - pbc_trans->z;
	  }
	}
	  
      }else{
	  
	/* not doing PBC's, just copy coordinates */
	for(j=0;j<excl_nbnum;j++){
	  jat = excl_nbiat[j];
	  xje[j] = xa[jat];
	  yje[j] = ya[jat];
	  zje[j] = za[jat];
	}
	
      }
#endif
      
    }
    
    *nbnum = *nbnum + excl_nbnum;

    /*
    printf("iat=%d nbnum=%d: ",iat, *nbnum);
    for(j=0;j<*nbnum;j++){
      jat = nbiat[j];
      if(!nbdata[j]){
	printf("NULL nbdata for %d %d\n",iat,jat);
      }
      pbc_trans = (NeighVector *)nbdata[j];
      printf("%p ",pbc_trans);
    }
    printf("\n");
    */
    
    return AGBNP_OK;
  }
    
  /* if no neighbor list do not apply cutoffs */
  *nbnum = 0;
  for(j=iat+1;j<agb->natoms;j++){
    nbiat[*nbnum] = j;
    *nbnum += 1;
  }
  
  return AGBNP_OK;
} 


#ifdef USE_PBC
 int agbnp_neighbor_lists_pbc(
#else
 int agbnp_neighbor_lists(
#endif
			  AGBNPdata *agb, AGBworkdata *agbw,
			  float_a *x, float_a *y, float_a *z){

  int nnl; /* neighbor list counter for near heavy-heavy (d<r1+r2)*/
  int nnlrc; /* neighbor list counter for far heavy-heavy (r1+r2<d<rc) */
  int iat,jat,j,nlsize=0,nbnum;
  float_a dx, dy, dz, d2, u;
  float_a nboffset = 1.6; /* offset for neighbor list distance test */
  float_a nlsize_increment = 1.2;
  int error = 0, res;

  int natoms = agb->natoms;
  float_a *r = agb->r;
  NeighList *near_nl = agbw->near_nl;
  NeighList *far_nl = agbw->far_nl;
  int *nbiat = agbw->nbiat;
#ifdef USE_PBC
  void **nbdata = agbw->nbdata;
  float_a *xj = agbw->xj;
  float_a *yj = agbw->yj;
  float_a *zj = agbw->zj;
#endif
  int *isheavy = agbw->isheavy;
  int do_frozen = agb->do_frozen;
  int *isbfrozen = agb->agbw->isbfrozen;
  float_a *galpha = agb->agbw->galpha;
  float_a *gprefac = agb->agbw->gprefac;
#ifdef USE_ALT
  int doalt = agb->doalt;
  int *alt_site = agb->alt_site;
  int *alt_id = agb->alt_id;
#endif
  int nsym = agb->nsym;
  int nearflag;
  float_a gvol;
  float_a *nl_r2v = agbw->nl_r2v;
  int *nl_indx = agbw->nl_indx;

  /* reset neighbor lists */
  memset(near_nl->nne, 0, natoms*sizeof(int));
  memset(far_nl->nne, 0, natoms*sizeof(int));

#ifdef USE_PBC
  if(agb->dopbc){
    if(!(near_nl->data && far_nl->data)){
      fprintf(stderr, "agbnp_neighbor_lists(): internal error: PBC on but neighbor lists lack PBC data handles.\n");
      return AGBNP_ERR;
    }
  }
#endif

  /* constructs  near and far neighbor lists */
  nnl = 0;
  nnlrc = 0;
#pragma omp for schedule(static,1)
  for(iat=0;iat<natoms;iat++){
    if(error) continue;
#ifdef USE_PBC
    if(agbnp_getatomneighbors_pbc(agb, iat, &nbnum, nsym*natoms,
			      nbiat, xj, yj, zj, nbdata)< 0){
      error = 1;
      continue;
    }
#else
    if(agbnp_getatomneighbors(agb, iat, &nbnum, natoms, nbiat)< 0){
      error = 1;
      continue;
    }
#endif

    if(isheavy[iat]){
      while(nnl + nsym*natoms >= near_nl->neighl_size){
	nlsize = mymax(nlsize_increment*near_nl->neighl_size, nnl + nsym*natoms);
	if(nblist_reallocate_neighbor_list(near_nl,natoms,nlsize) != NBLIST_OK){
	  error = 2;
	  continue;
	}
      }
      while(nnlrc + nsym*natoms >= far_nl->neighl_size){
	nlsize = mymax(nlsize_increment*far_nl->neighl_size,nnlrc + nsym*natoms);
	if(nblist_reallocate_neighbor_list(far_nl,natoms,nlsize) != NBLIST_OK){
	  error = 3;
	  continue;
	}
      }
      /* constructs  neighbor lists for atom iat */
      near_nl->nne[iat] = 0;  /* reset number of near neighbors for atom iat */
      far_nl->nne[iat] = 0; /* reset number of far neighbors for atom iat */
      /*set pntr to beg.of near neigh.list of atm iat*/
      near_nl->neighl[iat] = &(near_nl->neighl1[nnl]); 
      /*set pntr to beg.of far neigh.list of atm iat*/
      far_nl->neighl[iat] = &(far_nl->neighl1[nnlrc]); 
#ifdef USE_PBC
      if(agb->dopbc){
	near_nl->data_index[iat] = &(near_nl->data_index1[nnl]);
	far_nl->data_index[iat] = &(far_nl->data_index1[nnlrc]);
      }
#endif
      for(j=0;j<nbnum;j++){
#ifdef USE_PBC
	jat = nbiat[j] % natoms;
#else
	jat = nbiat[j];
#endif
	if(isheavy[jat]){
	  if(do_frozen && isbfrozen[iat] && isbfrozen[jat]) continue;
#ifdef USE_ALT
	  /*remove pairs in the same site in different alternate confrmatins */
	  if(doalt && (alt_site[iat]==alt_site[jat]) && 
	     (alt_id[iat]!=alt_id[jat])) continue;
#endif
	  /* insert in either far or near neighbor list */
	  nearflag = 1;
#ifdef USE_PBC
	  dx = xj[j] - x[iat];
	  dy = yj[j] - y[iat];
	  dz = zj[j] - z[iat];
#else
	  dx = x[jat] - x[iat];
	  dy = y[jat] - y[iat];
	  dz = z[jat] - z[iat];
#endif
	  d2 = dx*dx + dy*dy + dz*dz;
	  u = (r[iat]+r[jat])*nboffset;
	  /* include only near neighbors with non-zero overlap. First do
	     a simple distance test to save time, then calculate 2-body 
	     Gaussian overlap. Note that ogauss_2body_smpl() does not apply 
	     volume switching function, so that gvol < AGBNP_MIN_VOL means
	     actually a zero overlap volume */
	  if(d2>u*u){
	    nearflag = 0;
	  }else{
	    gvol = ogauss_2body_smpl(d2, gprefac[iat], gprefac[jat], 
				galpha[iat], galpha[jat]);
	    if(gvol<AGBNP_MIN_VOLA){
	      nearflag = 0;
	    }
	  }
	  if(nearflag){
	    /* jat is a near neighbor */
	    near_nl->neighl1[nnl] = nbiat[j];  /* place jat in neighbor list */
#ifdef USE_PBC
	    if(agb->dopbc){
	      near_nl->data_index1[nnl] = nbdata[j];
	    }
#endif
	    nnl += 1;                       /* updates neighbor list counter */
	    near_nl->nne[iat] += 1;     /* updates counter for iat neighbors */
	  }else{
	    /* jat is a far neighbor */
	    far_nl->neighl1[nnlrc] = nbiat[j]; /* place jat in neighbor list */
#ifdef USE_PBC
	    if(agb->dopbc){
              far_nl->data_index1[nnlrc] = nbdata[j];
            }
#endif
	    nnlrc += 1;                     /* updates neighbor list counter */
	    far_nl->nne[iat] += 1;      /* updates counter for iat neighbors */
	  }
	}else{
	  if(do_frozen && isbfrozen[jat] && isbfrozen[iat]) continue;
	  /* insert hydrogen in far neighbor list */
	  far_nl->neighl1[nnlrc] = nbiat[j];  /* place jat in neighbor list */
#ifdef USE_PBC
	  if(agb->dopbc){
	    far_nl->data_index1[nnlrc] = nbdata[j];
	  }
#endif
	  nnlrc += 1;                       /* updates neighbor list counter */
	  far_nl->nne[iat] += 1;        /* updates counter for iat neighbors */
	}
      }

    }else{
      /* neighbor list for hydrogens (no near list for hydrogens) */
      if(nnlrc + nsym*natoms > far_nl->neighl_size){
	nlsize = mymax(nlsize_increment*far_nl->neighl_size,nnlrc + nsym*natoms);
	if(nblist_reallocate_neighbor_list(far_nl,natoms,nlsize) != NBLIST_OK){
	  error = 3;
	  continue;
	}
      }
      /* constructs  neighbor lists for hydrogen atom iat */
      near_nl->nne[iat] = 0; 
      far_nl->nne[iat] = 0;  /* reset number of neighbors for atom iat */
      /*set pntr to beg.of neigh.list of atm iat*/
      near_nl->neighl[iat] = &(near_nl->neighl1[nnl]); 
      far_nl->neighl[iat] = &(far_nl->neighl1[nnlrc]); 
#ifdef USE_PBC
      if(agb->dopbc){
	near_nl->data_index[iat] = &(near_nl->data_index1[nnl]);
	far_nl->data_index[iat] = &(far_nl->data_index1[nnlrc]);
      }
#endif
      for(j=0;j<nbnum;j++){
#ifdef USE_PBC
	jat = nbiat[j] % natoms;
#else
	jat = nbiat[j];
#endif
	if(do_frozen && isbfrozen[iat] && isbfrozen[jat]) continue;
#ifdef USE_ALT
	if(doalt && (alt_site[iat]==alt_site[jat]) && 
	   (alt_id[iat]!=alt_id[jat])) continue;
#endif
	/* jat is a neighbor */
	far_nl->neighl1[nnlrc] = nbiat[j];   /* place jat in neighbor list */
#ifdef USE_PBC
	if(agb->dopbc){
	  far_nl->data_index1[nnlrc] = nbdata[j];
	}
#endif
	nnlrc += 1;                       /* updates neighbor list counter */
	far_nl->nne[iat] += 1;        /* updates counter for iat neighbors */
      }

#define AGBNP_NBL_SORT
#ifdef AGBNP_NBL_SORT
      if(near_nl->nne[iat] > 0){
	/* order near_nl in ascending order of distance */
	agbnp_fsortindx(near_nl->nne[iat], nl_r2v, nl_indx );
	agbnp_nblist_reorder(agbw, near_nl, iat, nl_indx);
      }
#endif
      


    }
  }

  if(error==1){
    fprintf(stderr, "agbnp_neighbor_lists(): error in agbnp_getatomneighbors().\n");
    return AGBNP_ERR;
  }else  if(error==2){
    fprintf(stderr,"agbnp_neighbor_lists(): unable to (re)allocate near_nl neighbor list (natoms=%d, size=%d)\n",natoms, nlsize);
    return AGBNP_ERR;
  }else if(error==3){
    fprintf(stderr,"agbnp_neighbor_lists(): unable to (re)allocate far_nl neighbor list (natoms=%d, size=%d)\n",natoms, nlsize);
    return AGBNP_ERR;
  }

  /* (re)allocation of 2-body volume cache */
  if(!agbw->v2cache || nnl > agbw->nv2cache){
    agbw->v2cache = (float  *)realloc(agbw->v2cache, nnl*sizeof(float ));
    if(!agbw->v2cache){
      fprintf(stderr, "agbnp_neighbor_lists(): fatal error: can't allocate memory for v2cache (%d float_as)!\n", nnl);
      return AGBNP_ERR;
    }
    agbw->nv2cache = nnl;
  }

  /* (re)allocation of i4() memory cache */
  if(!agbw->q4cache || 4*(nnl+nnlrc) > agbw->nq4cache){
    agbw->nq4cache = 4*(nnl+nnlrc);
    agbw->q4cache = (float  *)realloc(agbw->q4cache, agbw->nq4cache*sizeof(float ));
    if(!agbw->q4cache){
      fprintf(stderr, "agbnp_neighbor_lists(): fatal error: can't allocate memory for q4cache (%d float_as)!\n", agbw->nq4cache);
      return AGBNP_ERR;
    }
  }

  /*
  for(iat=0;iat<agb->natoms;iat++){
    NeighVector *pbc_trans;
    printf("iat=%d n=%d:",iat,far_nl->nne[iat]);
    for(j=0;j<far_nl->nne[iat];j++){
      jat = far_nl->neighl[iat][j];
      pbc_trans = (NeighVector *)far_nl->data_index[iat][j];
      printf("%d ",jat);
      printf("(%lf %lf %lf) ",pbc_trans->x,pbc_trans->y,pbc_trans->z);
    }
    printf("\n");
  }
  */
  
  return AGBNP_OK;
}

/* constructs  neighbor list  for surface area calculation */
#ifdef USE_PBC
 int agbnp_neighbor_list_cavity_pbc(
#else
 int agbnp_neighbor_list_cavity(
#endif
				AGBNPdata *agb, CAVworkdata *cavw,
				float_a *x, float_a *y, float_a *z){
  int nnl = 0, iat, nbnum, nlsize=0, j, jat;
  float_a nlsize_increment = 1.2;
  float_a dx, dy, dz, d2, u;
  float_a nboffset = 2.2; /* offset for neighbor list distance test */
  int error = 0, res;

  int natoms = agb->natoms;
  float_a *r = agb->rcav;
  NeighList *near_nl = cavw->near_nl;
  int *nbiat = cavw->nbiat;
#ifdef USE_PBC
  void **nbdata = cavw->nbdata;
  float_a *xj = cavw->xj;
  float_a *yj = cavw->yj;
  float_a *zj = cavw->zj;
#endif
  int *isheavy = cavw->isheavy;
  int do_frozen = agb->do_frozen;
  int *isvfrozen = agb->cavw->isvfrozen;
  int *isvvfrozen = agb->cavw->isvvfrozen;
  float_a *galpha = agb->cavw->galpha;
  float_a *gprefac = agb->cavw->gprefac;
#ifdef USE_ALT
  int doalt = agb->doalt;
  int *alt_site = agb->alt_site;
  int *alt_id = agb->alt_id;
#endif
  int nearflag;
  float_a gvol;
  float_a *nl_r2v = cavw->nl_r2v;
  int *nl_indx = cavw->nl_indx;

  /* reset neighbor list */
  memset(near_nl->nne, 0, natoms*sizeof(int));

#pragma omp for schedule(static,1)
  for(iat=0;iat<natoms;iat++){
    if(isheavy[iat]){

      if(error) continue;

      /* constructs  neighbor lists for atom iat */
      near_nl->nne[iat] = 0;  /* reset number of near neighbors for atom iat */
      /*set pntr to beg.of near neigh.list of atm iat*/
      near_nl->neighl[iat] = &(near_nl->neighl1[nnl]);
#ifdef USE_PBC 
      if(agb->dopbc){
	near_nl->data_index[iat] = &(near_nl->data_index1[nnl]);
      }
#endif
      /* can skip this atom if vvfrozen */
      if(do_frozen && isvvfrozen[iat]) continue;

#ifdef USE_PBC
      res = agbnp_getatomneighbors_pbc(agb, iat, &nbnum, agb->nsym*natoms,
				       nbiat, xj, yj, zj, nbdata);
#else
      res = agbnp_getatomneighbors(agb, iat, &nbnum, agb->natoms, nbiat);
#endif
      if(res < 0){
	error = 1;
	continue;
      }

      while(nnl + agb->nsym*natoms >= near_nl->neighl_size){
	nlsize = mymax(nlsize_increment*near_nl->neighl_size,nnl + agb->nsym*natoms);
	if(nblist_reallocate_neighbor_list(near_nl,natoms,nlsize) !=NBLIST_OK){
	  error = 2;
	  continue;
	}
      }
      for(j=0;j<nbnum;j++){
#ifdef USE_PBC
	jat = nbiat[j] % natoms;
#else
	jat = nbiat[j];
#endif
#ifdef USE_ALT
	/*remove pairs in the same site in different alternate conformations */
	if(doalt && (alt_site[iat]==alt_site[jat]) && 
	   (alt_id[iat]!=alt_id[jat])) continue;
#endif
	if(isheavy[jat]){
	  nearflag = 1;
#ifdef USE_PBC
	  dx = xj[j] - x[iat];
	  dy = yj[j] - y[iat];
	  dz = zj[j] - z[iat];
#else
	  dx = x[jat] - x[iat];
	  dy = y[jat] - y[iat];
	  dz = z[jat] - z[iat];
#endif
	  d2 = dx*dx + dy*dy + dz*dz;
	  u = (r[iat]+r[jat])*nboffset;
	  /* include only near neighbors with non-zero overlap. First do
	     a simple distance test to save time, then calculate 2-body 
	     Gaussian overlap. Note that ogauss_2body_smpl() does not apply 
	     volume switching function, so that gvol < AGBNP_MIN_VOL means
	     actually a zero overlap volume */
	  if(d2>u*u){
	    nearflag = 0;
	  }else{
	    gvol = ogauss_2body_smpl(d2, gprefac[iat], gprefac[jat], 
				galpha[iat], galpha[jat]);
	    if(gvol<AGBNP_MINC_VOLA){
	      nearflag = 0;
	    }
	  }
	  if(nearflag){
	    /* insert in neighbor list */
	    near_nl->neighl1[nnl] = nbiat[j]; /* place jat in neighbor list */
#ifdef USE_PBC
	    if(agb->dopbc){
	      near_nl->data_index1[nnl] = nbdata[j];
	    }
#endif
	    nnl += 1;                       /* updates neighbor list counter */
	    nl_r2v[near_nl->nne[iat]] = d2; /* store distance for reordering */
	    near_nl->nne[iat] += 1;     /* updates counter for iat neighbors */
	  }
	}

#ifdef AGBNP_NBL_SORT
	if(near_nl->nne[iat] > 0){
	  /* order near_nl in ascending order of distance */
	  agbnp_fsortindx(near_nl->nne[iat], nl_r2v, nl_indx );
	  agbnp_nblist_reorder_cavity(cavw, near_nl, iat, nl_indx);
	}
#endif

      }
    }
  }

  if(error==1){
    fprintf(stderr, "agbnp_neighbor_list_cavity(): error in neighbors agbnp_getatomneighbors().\n");
    return AGBNP_ERR;
  }else if(error==2){
    fprintf(stderr,"agbnp_neighbor_list_cavity(): unable to (re)allocate near_nl neighbor list (natoms=%d, size=%d)\n",natoms, nlsize);
    return AGBNP_ERR;
  }

  return AGBNP_OK;
}
