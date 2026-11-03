#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

#include "libnblist.h"
#include "agbnp.h"

/* append a trailing underscore to symbols only if Fortran compiler does */
#ifdef FC_UNDERSCORES
  #define FORTRAN(s) s##_
#else
  #define FORTRAN(s) s
#endif


#if defined(BINDCRAY)
#define agbnpfdelete AGBNPFDELETE
#define agbnpf_init AGBNPF_INIT
#define agbnpf_reinit_on AGBNPF_REINIT_ON
#define agbnpfinitfrozen AGBNPFINITFROZEN
#define agbnpf AGBNPF
#define agbnpf_de AGBNPF_DE
#endif

#ifndef AGBNP_MAXSTRING
#define AGBNP_MAXSTRING 1024
#endif

#ifndef AGBNP_MAXTYPE
#define AGBNP_MAXTYPE 300
#endif

typedef struct {
  int ntypes;            /* number of atom types */
  int *type_index;       /* agbnp fitting type index */
  char **symbol;         /* atom type symbols */
  int *typeid;           /* atom type ids */
  double *radius;        /* vdw radius (defines solute cavity) */
  double *igamma;        /* gamma for ideal np energy */
  double *ialpha;        /* alpha for ideal np energy */
  double *idelta;        /* delta for ideal np energy */
  double *sgamma;        /* gamma for correction energy */
  double *salpha;        /* alpha for correction energy */
  double *sdelta;        /* delta for correction energy */
  double *ab;            /* GB distance correction factor */
  int *hbtype;           /* type of HB donor/acceptor (=0 is HB inactive) */
  double *hbcorr;        /* strength of HB to solvent correction */ 
  double dielectric_in;  /* solute dielectric constant */
  double dielectric_out; /* solvent dielectric constant */
} np_param_t;

typedef struct agbnpfw {
  /* work data for agbnpf */
  int natoms;
  int version;   /* agbnp version (1 or 2) */
  int agbnp_tag; /* libagbnp tag */
  int mmnblist_tag, mmnblistg_tag; /* mmnblist tag */
  NeighList *neigh_list, *excl_neigh_list; /* ptrs to Verlet neighbor lists */
  int nheavyat, *iheavyat; /* number of and list of heavy atoms */
  int nhydrogen, *ihydrogen; /* number of and list of hydrogen atoms */
  int ndummy, *idummy; /* number of and list of dummy atoms */
  int *isfrozen; /* 1 if atom is frozen, 0 otherwise */
  int dopbc;     /* > 0 if doing PBC's */
  double *x, *y, *z; /* atom coordinates */
  double *r;         /* atom radii */
  double *charge;    /* atom charges */
  double *sp;        /* scaled volume factors returned by agbnp() */
  double *br;        /* born radii returned by agbnp() */
  double mol_volume;
  double *surf_area; /* atoms vdw surface areas returned by agbnp */
  double (*dgbdr)[3]; /* positional derivatives of GB energy */
  double (*dvwdr)[3]; /* positional derivatives of NP energy */
  double (*decav)[3]; /* positional derivatives of cavity energy */
  double (*dehb)[3];  /* positional derivatives of HB correction energy */
  double *gamma;      /* gamma np parameters */
  double *alpha;      /* alpha np parameters */
  double *delta;      /* delta np parameters */
  double *sgamma;      /* gamma np correction parameters */
  double *salpha;      /* alpha np correction parameters */
  double *sdelta;      /* delta np correction parameters */
  double *ab;          /* GB pair distance correction factor */
  int *hbtype;         /* type of HB donor/acceptor (=0 is HB inactive) */
  double *hbcorr;      /* strength of HB to solvent correction */ 
  np_param_t *np;      /* non-polar parameters by type */
  NeighList *conntbl;  /* connection table */
  double *occupancy;   /* occupancy of alternate conformation */
  int *alt_site;       /* id of alternate conf. site (starting from 1)*/
  int *alt_id;         /* id of alternate conf. in site (starting from 1)*/
  int ener_init;       /* flag set to > 0 to reinitialize AGBNP energy evaluation */
  int e_count;         /* number of calls to agbnpf() */
  int e_count_reinit;  /* number of calls in between re-initializations (see e_count) */
  int *typeid;           /* atom type ids */
} AGBNPfw;


void agbnpf_delete_np(np_param_t *np){
  int i;
  if(!np){
    return;
  }
  if(np->type_index) { free(np->type_index); np->type_index = NULL; }
  if(np->symbol){
    for(i=0;i<AGBNP_MAXTYPE;i++){
      if(np->symbol[i]) free(np->symbol[i]);
    }
    free(np->symbol);
    np->symbol = NULL;
  }
  if(np->radius) { free(np->radius); np->radius = NULL; }
  if(np->typeid) { free(np->typeid); np->typeid = NULL; }
  if(np->igamma) { free(np->igamma); np->igamma = NULL; }
  if(np->ialpha) { free(np->ialpha); np->ialpha = NULL; }
  if(np->idelta) { free(np->idelta); np->idelta = NULL; }
  if(np->sgamma) { free(np->sgamma); np->sgamma = NULL; }
  if(np->salpha) { free(np->salpha); np->salpha = NULL; }
  if(np->sdelta) { free(np->sdelta); np->sdelta = NULL; }
  if(np->ab) { free(np->ab); np->ab = NULL; }
  if(np->hbtype) { free(np->hbtype); np->hbtype = NULL; }
  if(np->hbcorr) { free(np->hbcorr); np->hbcorr = NULL; }
  return;
}

static np_param_t *agbnp_param( AGBNPfw *w, char *filename ){

  np_param_t *np;
  FILE *infile;
  int n = 0;
  char *ch;
  char line[AGBNP_MAXSTRING];
  double igamma, ialpha, idelta;
  double sgamma, salpha, sdelta;
  double ab = 0.0;
  double hbcorr = 0.0;
  int hbtype = 0;
  double radius;
  char atom_symbl[AGBNP_MAXSTRING];
  int type_index;
  double diel;
  int found_diel_in = 0;
  int found_diel_out = 0;
  int nitems, nitems_expected;
  static const char dielectric_out_tag[] = "dielectric_out";
  static const char dielectric_in_tag[] = "dielectric_in";
  char tag[AGBNP_MAXSTRING];
  extern FILE *fopen_npsolv(const char* filename, const char* mode);

  /* allocate memory for np structure */
  np = (np_param_t *)malloc(sizeof(np_param_t));
  np->type_index = (int *)calloc(AGBNP_MAXTYPE, sizeof(int));
  np->symbol = (char **)calloc(AGBNP_MAXTYPE, sizeof(char *));
  np->typeid = (int *)calloc(AGBNP_MAXTYPE, sizeof(int));
  np->radius = (double *)calloc(AGBNP_MAXTYPE, sizeof(double));
  np->igamma = (double *)calloc(AGBNP_MAXTYPE, sizeof(double));
  np->ialpha = (double *)calloc(AGBNP_MAXTYPE, sizeof(double));
  np->idelta = (double *)calloc(AGBNP_MAXTYPE, sizeof(double));
  np->sgamma = (double *)calloc(AGBNP_MAXTYPE, sizeof(double));
  np->salpha = (double *)calloc(AGBNP_MAXTYPE, sizeof(double));
  np->sdelta = (double *)calloc(AGBNP_MAXTYPE, sizeof(double));
  np->ab = (double *)calloc(AGBNP_MAXTYPE, sizeof(double));
  np->hbtype = (int *)calloc(AGBNP_MAXTYPE, sizeof(int));
  np->hbcorr = (double *)calloc(AGBNP_MAXTYPE, sizeof(double));

  //  for(int i = 0 ; i < AGBNP_MAXTYPE ; i++){
  //  np->symbol[i] = (char *)calloc(5, sizeof(char));
  // }
  
  /* **** location? */
  if((infile = fopen(filename, "r")) == NULL){
    fprintf(stderr,"agbnp_param: cannot open file %s for reading. \n",
	    filename);
    return((np_param_t *)NULL);
  }

  /* read dielectric constants */
  while(!found_diel_in || !found_diel_out){
    if( (ch = (char *)fgets(line, AGBNP_MAXSTRING, infile)) == (char *)NULL){
      fprintf(stderr,"agbnp_param: unexpected end of file %s.\n",
	      filename);
      return((np_param_t *)NULL);
    }
    if(line[0] == '#') continue; /* skip comments */
    if(sscanf(line, "%s %lf", tag, &diel) == 2){
      if(strcmp(tag, dielectric_in_tag) == 0){
	found_diel_in = 1;
	np->dielectric_in = diel;
      }else if(strcmp(tag, dielectric_out_tag) == 0){
        found_diel_out = 1;
        np->dielectric_out = diel;
      }
    }
  }

  /*  printf("Reading file %s\n",filename); */
  while((ch = (char *)fgets(line, AGBNP_MAXSTRING, infile))!= (char *) NULL){
    if(line[0] == '#') continue; /* skip comments */
    if(w->version == 2){
      nitems_expected = 11;
      nitems = 
	sscanf(line, "%d %s %lf %lf %lf %lf %lf %lf %lf %d %lf", 
	       &type_index, 
	       atom_symbl, &radius, 
	       &igamma, &ialpha, &idelta, &sgamma, &salpha, &sdelta, 
	       &hbtype, &hbcorr);
    }else{
      nitems_expected = 10;
      nitems =
	sscanf(line, "%d %s %lf %lf %lf %lf %lf %lf %lf %lf", 
	       &type_index, 
	       atom_symbl, &radius, 
	       &igamma, &ialpha, &idelta, &sgamma, &salpha, &sdelta, 
	       &ab);
    }
    if(nitems == nitems_expected){
      np->type_index[n] = type_index;
      np->symbol[n] = strdup(atom_symbl);
      np->typeid[n] = type_index;
      np->radius[n] = radius;      /* vdw radius */
      np->igamma[n] = igamma;
      np->ialpha[n] = ialpha;
      np->idelta[n] = idelta;
      np->sgamma[n] = sgamma;
      np->salpha[n] = salpha;
      np->sdelta[n] = sdelta;
      np->ab[n] = ab;              /* GB distance correction factor */
      np->hbtype[n] = hbtype;      /* HB corrections */
      np->hbcorr[n] = hbcorr;
      n += 1;
    }else{
      fprintf(stderr,"agbnp_param: error reading data from file %s.\n",
	      filename);
      return((np_param_t *)NULL);
    }
  }

  np->ntypes = n;

#ifdef  AGBNP_WRITE_APARAMS
  for(n=0;n<np->ntypes;n++){
    if(w->version == 2){
      printf("PAR: %d %25s %lf %lf %lf %lf %lf %lf %lf %d %lf\n",
	     np->type_index[n], np->symbol[n],
	     np->radius[n], 
	     np->igamma[n], np->ialpha[n], np->idelta[n],
	     np->sgamma[n], np->salpha[n], np->sdelta[n],
	     np->hbtype[n], np->hbcorr[n]);
    }else{
      printf("PAR: %d %25s %lf %lf %lf %lf %lf %lf %lf %lf\n",
	     np->type_index[n], np->symbol[n],
	     np->radius[n], 
	     np->igamma[n], np->ialpha[n], np->idelta[n],
	     np->sgamma[n], np->salpha[n], np->sdelta[n],
	     np->ab[n]);
    }
  }
#endif

  
  fclose(infile);
  return(np);
}


void agbnpf_delete_w(AGBNPfw *w){
  if(!w){
    return;
  }
  if(w->agbnp_tag >= 0) {
    agbnp_delete(w->agbnp_tag);
    w->agbnp_tag = -1;
  }
  if(w->x) { free(w->x); w->x = NULL; }
  if(w->y) { free(w->y); w->y = NULL; }
  if(w->z) { free(w->z); w->z = NULL; }
  if(w->r) { free(w->r); w->r = NULL; }
  if(w->charge) { free(w->charge); w->charge = NULL; }
  if(w->iheavyat) { free(w->iheavyat); w->iheavyat = NULL; }
  if(w->ihydrogen) { free(w->ihydrogen); w->ihydrogen = NULL; }
  if(w->idummy) { free(w->idummy); w->idummy = NULL; }
  if(w->isfrozen){ free(w->isfrozen) ; w->isfrozen = NULL; }
  if(w->gamma) { free(w->gamma); w->gamma = NULL; }
  if(w->alpha) { free(w->alpha); w->alpha = NULL; }
  if(w->delta) { free(w->delta); w->delta = NULL; }
  if(w->sgamma) { free(w->sgamma); w->sgamma = NULL; }
  if(w->salpha) { free(w->salpha); w->salpha = NULL; }
  if(w->sdelta) { free(w->sdelta); w->sdelta = NULL; }
  if(w->ab) { free(w->ab); w->ab = NULL; }
  if(w->hbtype) { free(w->hbtype); w->hbtype = NULL; }
  if(w->hbcorr) { free(w->hbcorr); w->hbcorr = NULL; }
  if(w->sp) { free(w->sp); w->sp = NULL; }
  if(w->br) { free(w->br); w->br = NULL; }
  if(w->surf_area) { free(w->surf_area); w->surf_area = NULL; }
  if(w->dgbdr) { free(w->dgbdr); w->dgbdr = NULL; }
  if(w->dvwdr) { free(w->dvwdr); w->dvwdr = NULL; }  
  if(w->decav) { free(w->decav); w->decav = NULL; }
  if(w->dehb) { free(w->dehb); w->dehb = NULL; }
  if(w->np) { agbnpf_delete_np(w->np); w->np = NULL; }
  if(w->typeid) { free(w->typeid); w->typeid = NULL; }
  if(w->conntbl){ 
    nblist_delete_neighbor_list(w->conntbl); 
    free(w->conntbl) ; 
    w->conntbl = NULL; 
  }
  if(w->occupancy) { free(w->occupancy); w->occupancy = NULL; }
  if(w->alt_site) { free(w->alt_site); w->alt_site = NULL; }
  if(w->alt_id) { free(w->alt_id); w->alt_id = NULL; }
  return;
}

void agbnpf_reset_w(AGBNPfw *w){
  if(!w){
    return;
  }
  w->natoms = 0;
  w->agbnp_tag = -1;
  w->mmnblist_tag = -1;
  w->mmnblistg_tag = -1;
  w->neigh_list = NULL;
  w->excl_neigh_list = NULL;
  w->nheavyat = w->nhydrogen = w->ndummy = 0;
  w->x = NULL;
  w->y = NULL;
  w->z = NULL;
  w->r = NULL;
  w->charge = NULL;
  w->iheavyat = NULL;
  w->ihydrogen = NULL;
  w->idummy = NULL;
  w->isfrozen = NULL;
  w->dopbc = 0;
  w->gamma = NULL;
  w->alpha = NULL;
  w->delta = NULL;
  w->sgamma = NULL;
  w->salpha = NULL;
  w->sdelta = NULL;
  w->ab = NULL;
  w->hbtype = NULL;
  w->hbcorr = NULL;
  w->sp = NULL;
  w->br = NULL;
  w->surf_area = NULL;
  w->dgbdr = NULL;
  w->dvwdr = NULL;  
  w->decav = NULL;
  w->dehb = NULL;
  w->np = NULL;
  w->conntbl = NULL;
  w->occupancy = NULL;
  w->alt_site = NULL;
  w->alt_id = NULL;
  w->typeid = NULL;
  return;
}

int agbnpf_allocate_w(AGBNPfw *w){
  int natoms = w->natoms;
  w->x = (double *)realloc(w->x, natoms*sizeof(double));
  w->y = (double *)realloc(w->y, natoms*sizeof(double));
  w->z = (double *)realloc(w->z, natoms*sizeof(double));
  w->r = (double *)realloc(w->r, natoms*sizeof(double));
  w->charge = (double *)realloc(w->charge, natoms*sizeof(double));
  w->iheavyat = (int *)realloc(w->iheavyat, natoms*sizeof(int));
  w->ihydrogen = (int *)realloc(w->ihydrogen, natoms*sizeof(int));
  w->idummy = (int *)realloc(w->idummy, natoms*sizeof(int));
  w->sp = (double *)realloc(w->sp, natoms*sizeof(double));
  w->br = (double *)realloc(w->br, natoms*sizeof(double));
  w->surf_area = (double *)realloc(w->surf_area, natoms*sizeof(double));
  w->dgbdr = (double (*)[3])realloc(w->dgbdr, natoms*sizeof(double [3]));
  w->gamma = (double *)realloc(w->gamma, natoms*sizeof(double));
  w->alpha = (double *)realloc(w->alpha, natoms*sizeof(double));
  w->delta = (double *)realloc(w->delta, natoms*sizeof(double));
  w->sgamma = (double *)realloc(w->sgamma, natoms*sizeof(double));
  w->salpha = (double *)realloc(w->salpha, natoms*sizeof(double));
  w->sdelta = (double *)realloc(w->sdelta, natoms*sizeof(double));
  w->ab = (double *)realloc(w->ab, natoms*sizeof(double));
  w->hbtype = (int *)realloc(w->hbtype, natoms*sizeof(int));
  w->hbcorr = (double *)realloc(w->hbcorr, natoms*sizeof(double));
  w->dvwdr = (double (*)[3])realloc(w->dvwdr, natoms*sizeof(double [3]));
  w->decav = (double (*)[3])realloc(w->decav, natoms*sizeof(double [3]));
  w->dehb = (double (*)[3])realloc(w->dehb, natoms*sizeof(double [3]));
  w->typeid = (int *)realloc(w->typeid, natoms*sizeof(int));
  if(!(w->x && w->y && w->z && w->r && w->charge && w->iheavyat &&
       w->ihydrogen && w->idummy && w->sp && w->br && 
       w->surf_area && w->dgbdr && w->gamma && w->alpha && w->delta && 
       w->sgamma && w->salpha && w->sdelta && w->ab &&
       w->dvwdr && w->decav && w->dehb && w->typeid)){
    fprintf(stderr,"agbnpf_allocate_data(): unable to allocate memory for one or more arrays.\n");
    return -1;
  }
  if(!(w->hbtype && w->hbcorr)){
    fprintf(stderr,"agbnpf_allocate_data(): unable to allocate memory for one or more arrays.\n");
    return -1;
  }

  memset(w->dgbdr,0,natoms*sizeof(double [3]));
  memset(w->dvwdr,0,natoms*sizeof(double [3]));
  memset(w->decav,0,natoms*sizeof(double [3]));
  memset(w->dehb, 0,natoms*sizeof(double [3]));

  return 1;
}

int agbnpf_assign_parameters(AGBNPfw *w, double * rx, double *ry, double *rz, 
			     int *typep, double *charge,
			     double *epsilonlj, double *sigmalj, 
			     int *ishydrogen){
  int iat, ispec, i, j, n;
  static const double sigmaw = 3.15365; /* LJ sigma of TIP4P water oxygen */
  static const double epsilonw = 0.155; /* LJ epsilon of TIP4P water oxygen */
  static const double rho = 0.033428;   /* water number density */
  double sigma, epsilon;
  double sigmaiw, epsiloniw;
  int type_found;
  double ws, s6;
  char atom_type[5];
  double pi = 3.14159265359;
  int nffield = 0;

  /*
  for(iat = 0; iat < w->natoms; iat++){
    printf("typeid: %d %d\n",iat+1,typep[iat]);
    } */

  for(iat = 0; iat < w->natoms; iat++){

      /* set coordinates */
      w->x[iat] = rx[iat];
      w->y[iat] = ry[iat];
      w->z[iat] = rz[iat];

      /* set charges */
      w->charge[iat] = charge[iat];

      /* retrieves LJ sigma and epsilon and combines with LJ water parameters*/
      sigma = sigmalj[iat];
      epsilon =  epsilonlj[iat];
      if(nffield > 1){ /* OPLS combining rules */
	sigma = sigma*sigma; /* impact stores sqrt(sigma) in sigtmp */
	epsilon = epsilon*epsilon/4.0; /* impact stores sqrt(4 epsilon) in eptmp */
	sigmaiw = sqrt(sigma*sigmaw);
	epsiloniw = sqrt(epsilon*epsilonw);
      }else{                          /* AMBER combining rules */
	sigma = 2.0*sigma;   /* impact stores sigma/2 in sigtmp */
	epsilon = epsilon*epsilon/4.0; /* impact stores sqrt(4 epsilon) in eptmp */
	sigmaiw = 0.5*(sigma+sigmaw);
	epsiloniw = sqrt(epsilon*epsilonw);
      }

      /* compare with known atom types and assign parameters */
      type_found = 0;
      for(n = 0 ; n < w->np->ntypes ; n++){
	if(w->np->typeid[n] == typep[iat]){ 
	  type_found = 1;
	  w->r[iat] = w->np->radius[n];
	  if(w->r[iat] <= 0.0){
	    fprintf(stderr,"agbnpf_assign_parameters(): fatal error: atoms with radius <= 0 are not allowed.\n");
	    return -1;
	  }
	  /* scaling factor for alpha based on LJ sigma and epsilon 
	     LJ weight: -16*pi*rho*epsilon*sigma^6/3 */
	  s6 = sigmaiw;
	  s6 = s6*s6 ; s6 = s6*s6*s6; /* sigmaiw^6 */
	  ws = - 16.0 * pi * rho * epsiloniw * s6 / 3.0;
	  /* ideal np parameters */
	  w->gamma[iat] =    w->np->igamma[n];
	  w->alpha[iat] = ws*w->np->ialpha[n];
	  w->delta[iat] =    w->np->idelta[n];
	  /* np correction parameters */
	  w->sgamma[iat] =    w->np->sgamma[n];
	  w->salpha[iat] = ws*w->np->salpha[n];
	  w->sdelta[iat] =    w->np->sdelta[n];
	  /* GB distance correction factor */
	  w->ab[iat] = w->np->ab[n];
	  if(w->version == 2){
	    /* HB corrections */
	    w->hbtype[iat] = w->np->hbtype[n];
	    w->hbcorr[iat] = w->np->hbcorr[n];
	  }
	  w->typeid[iat] = n;

	  break;
	}
      }
      if(!type_found){
	  fprintf(stderr,"agbnpf_assign_parameters(): atom type %s is an unknown agbnp atom type.\n",atom_type);
	  return -1;
      }



    /* set hydrogen flags */
      if(ishydrogen[iat]) {
	w->ihydrogen[w->nhydrogen++] = iat;
      }else{
	w->iheavyat[w->nheavyat++] = iat;
      }

    }



#ifdef AGBNP_WRITE_APARAMS
  iat = 0;
  for(iat= 0; iat<w->natoms; iat++){
      if(w->version == 2){
	printf("BPAR: %d %5s %f %f %f %f %f %f %f %f %d %f\n",
		iat+1, atom_type, w->r[iat],
		w->gamma[iat], w->alpha[iat], w->delta[iat],
		w->sgamma[iat], w->salpha[iat], w->sdelta[iat],
		w->ab[iat],
		w->hbtype[iat], w->hbcorr[iat]);
      }else{
	printf("BPAR: %d %d %25s %f %f %f %f %f %f %f %f\n",
	       iat+1, w->typeid[iat], w->np->symbol[w->typeid[iat]], w->r[iat],
		w->gamma[iat], w->alpha[iat], w->delta[iat],
		w->sgamma[iat], w->salpha[iat], w->sdelta[iat],
		w->ab[iat]);
      }
  }
#endif

  return 1;
}

static int agbnpf_set_frozen(AGBNPfw *w){
  int iat, ispec, i;
  int do_frozen = 0;

  /* set frozen atoms if at least one atom is frozen */
  /*
  for(ispec = 0, iat = 0; ispec<COMMON(impsize).nspec; ispec++){
    for(i=0; i<COMMON(impsize).natinsp[ispec]; i++, iat++){
      if(COMMON(freeze).ifrzat[ispec][i] == -1){
	do_frozen = 1;
      }
    }
  }
  */

  if(!do_frozen){
    if(w->isfrozen){
      free(w->isfrozen);
    }
    w->isfrozen = NULL;
    return 1;
  }

#ifdef NOTNOW

  w->isfrozen = (int *)realloc(w->isfrozen,w->natoms*sizeof(int));
  if(!w->isfrozen){
    fprintf(stderr,"agbnpf_set_frozen(): unable to allocate frozen flag array (%d ints)\n", w->natoms);
    return -1;
  }

  for(ispec = 0, iat = 0; ispec<COMMON(impsize).nspec; ispec++){
    for(i=0; i<COMMON(impsize).natinsp[ispec]; i++, iat++){
      w->isfrozen[iat] = 0;
      if(COMMON(freeze).ifrzat[ispec][i] == -1){
	w->isfrozen[iat] = 1;
      }
    }
  }

#endif /* NOTNOW */

  return 1;
}

#ifdef NOTNOW

static int agbnpf_set_conntbl(AGBNPfw *w){
  double nlsize_increment = 1.2;
  int nconnat = 4;
  int nlsize;
  int ispec, iat, i, j;
  int isize, jat, jatom;

  if(!w){
    return 1;
  }

  if(w->conntbl){
    nblist_delete_neighbor_list(w->conntbl); 
    w->conntbl = NULL;
  }

  w->conntbl = (NeighList *)calloc(1,sizeof(NeighList));
  nblist_reset_neighbor_list(w->conntbl);
  nlsize = nconnat * w->natoms;
  if(nblist_reallocate_neighbor_list(w->conntbl,w->natoms,nlsize) !=NBLIST_OK){
    fprintf(stderr, "agbnpf_set_conntbl(): unable to (re)allocate connection table (size = %d ints).\n", nlsize);
    return -1;
  }

  /* fill connection table */
  isize = 0;
  for(ispec = 0, iat = 0; ispec<COMMON(impsize).nspec; ispec++){
    for(i=0; i<COMMON(impsize).natinsp[ispec]; i++, iat++){
      w->conntbl->nne[iat] = 0;
      w->conntbl->neighl[iat] = &(w->conntbl->neighl1[isize]);
      while(isize + COMMON(strctp).nconn[ispec][i] >= w->conntbl->neighl_size){
	nlsize = nlsize_increment*w->conntbl->neighl_size;
	if(nblist_reallocate_neighbor_list(w->conntbl,w->natoms,nlsize) 
	   !=NBLIST_OK){
	  fprintf(stderr, "agbnpf_set_conntbl(): unable to (re)allocate connection table (size = %d ints).\n", nlsize);
	  return -1;
	}
      }
      for(j=0;j<COMMON(strctp).nconn[ispec][i];j++){
	jatom = COMMON(strctp).conntbl[j][ispec][i]; 
	jat = COMMON(impsize).natinsp[0]*ispec + jatom - 1;
	w->conntbl->nne[iat] += 1;
	w->conntbl->neighl1[isize] = jat;
	isize += 1;
      }
    }
  }

  return 1;
}

#endif

static AGBNPfw *agbnp_w = NULL;

void FORTRAN(agbnpfdelete)(){
  if(agbnp_w){
    agbnpf_delete_w(agbnp_w);
    free(agbnp_w);                         
    agbnp_w = NULL;
  }
}



void FORTRAN(agbnpf_init)(int *natoms, double *rx, double *ry, double * rz, 
			  int *typep, double *charge, 
			  double *epsilonlj, double *sigmalj,
			  int *ishydrogen,
			  int *rescode){
  /* initializes agbnp interface. Read parameters. Allocates work arrays. 
     Calls agbnp_initialize() and agbnp_new().
   */
  AGBNPfw *w;
  int ispec;
  int i,j, iat, iatom;
  static int libagbnp_initialized = 0;
  int isym,off, nsym = 1;
  int ssize;
  double *xs, *ys, *zs;
  double (*rot)[3][3] = NULL;

  if(agbnp_w){
    /* delete previous data */
    agbnpf_delete_w(agbnp_w);
    free(agbnp_w);
    agbnp_w = NULL;
  }

  agbnp_w = (AGBNPfw *)calloc(1,sizeof(AGBNPfw));
  w = agbnp_w;
  if(!agbnp_w){
    fprintf(stderr,"agbnpf_init(): unable to allocate agbnp_w.\n");
    *rescode = -1;
    return;
  }
  agbnpf_reset_w(agbnp_w);

  /* agbnp version */
  w->version = 1;

  /* sets up number of atoms */
  w->natoms = *natoms;
  if(w->natoms <= 0){
    fprintf(stderr,"agbnpf_init(): invalid number of atoms %d\n", w->natoms);
    *rescode = -1;
    return;
  }

  /* set PBC if needed */
  w->dopbc = 0;
  if(w->dopbc){
    nsym = 0;
  }

  /* allocate work arrays */
  if(agbnpf_allocate_w(agbnp_w) < 0){
    fprintf(stderr,"agbnpf_init(): error in agbnpf_allocate_w().\n");
    *rescode = -1;
    return;
  }

  /* read np parameters */
  if(w->version == 2){
    agbnp_w->np = agbnp_param( w, "agbnp2.param" );
  }else{
    agbnp_w->np = agbnp_param( w, "agbnp.param" );
  }
  if(!agbnp_w->np){
    fprintf(stderr,"agbnpf_init(): error in agbnp_param().\n");
    *rescode = -1;
    return;
  }

  /* assign np parameters and coordinates */
  if(agbnpf_assign_parameters(agbnp_w, rx, ry, rz, typep, charge, 
			      epsilonlj, sigmalj, ishydrogen) < 0){
    fprintf(stderr, "agbnpf_init(): error in agbnpf_assign_parameters().\n");
    *rescode = -1;
    return;
  }

  /* set frozen atoms */
  if(agbnpf_set_frozen(agbnp_w) < 0){
    fprintf(stderr, "agbnpf_init(): error in agbnpf_set_frozen().\n");
    *rescode = -1;
    return;
  }


     /* w->neigh_list = get_neighf_address() */
  w->neigh_list = NULL;
  w->excl_neigh_list = NULL;

  /* initialize libagbnp if necessary */
  if(!libagbnp_initialized){
    if(agbnp_initialize() != AGBNP_OK){
      fprintf(stderr, "agbnpf_init(): error in agbnp_initialize().\n");
      *rescode = -1;
      return;
    }
    libagbnp_initialized = 1;
  }

  /* set PBC's */
#ifdef NOTNOW
  if(w->dopbc) {
    if(w->mmnblist_tag < 0 || w->mmnblistg_tag < 0) {
      fprintf(stderr, "agbnpf_init(): error: a group-based neighbor list is needed with periodic boundary conditions.\n");
      *rescode = -1;
      return;
    }
    if(mmnblist_get_coord_ptr(w->mmnblist_tag, &ssize, &xs, &ys, &zs) 
       != AGBNP_OK){
      fprintf(stderr, "agbnpf_init(): error in mmnblist_get_coord_ptr().\n");
      *rescode = -1;
      return;
    }
    nsym = COMMON(boxcmn).nsym;
    if(nsym>1){

      /*
       for(isym=0;isym<nsym;isym++){
	for(iat=0;iat<w->natoms;iat++){
	  iatom = isym*(w->natoms+1)+iat+1;
	  printf("%d %d %lf %lf %lf\n",iatom,w->neigh_list->int2ext[isym*w->natoms+iat],xs[iatom],ys[iatom],zs[iatom]);
	}
      }
      */

      rot = (double (*)[3][3])calloc(nsym,sizeof(double [3][3]));
      if(!rot){
	fprintf(stderr, "agbnpf_init(): unable to allocate memory for rotation matrices.\n");
	*rescode = -1;
	return;
      }
      for(isym=0;isym<nsym;isym++){
        for(i=0;i<3;i++){
          for(j=0;j<3;j++){
            rot[isym][i][j] = COMMON(boxcmn).rot[isym][j][i];
          }
        }
      }

      /* print rot matrices */
      fprintf(stderr,"Rotational matrices:\n");
      for(isym=0;isym<nsym;isym++){
        fprintf(stderr,"Sym= %d\n",isym);
        for(i=0;i<3;i++){
          for(j=0;j<3;j++){
            fprintf(stderr,"%7.3f ", rot[isym][i][j]);
          }
          fprintf(stderr,"\n");
        }
      }



    }
  }else{
    ssize = 0;
    nsym = 0;
    xs = ys = zs = NULL;
    rot = NULL;
  }
#endif /* NOTNOW */


  /* create a libagbnp object */
  if(w->version == 2){
    /* version 2 */
#ifdef NOTNOW
    if(agbnp2_new(&(w->agbnp_tag), w->natoms, 
		  w->x, w->y, w->z, w->r, 
		  w->charge, w->np->dielectric_in, w->np->dielectric_out,
		  w->gamma, w->sgamma,
		  w->alpha, w->salpha,
		  w->delta, w->sdelta,
		  w->hbtype, w->hbcorr,
		  w->nhydrogen, w->ihydrogen, 
		  w->ndummy, w->idummy, w->isfrozen,
		  w->neigh_list, w->excl_neigh_list,
		  w->dopbc, nsym, ssize, xs, ys, zs, rot,
		  w->conntbl) != AGBNP_OK){

      fprintf(stderr, "agbnpf_init(): error in agbnp2_new().\n");
      *rescode = -1;
      return;
    }
#endif
  }else{
    /* version 1 */
    if(agbnp_new(&(w->agbnp_tag), w->natoms, 
		  w->x, w->y, w->z, w->r, 
		  w->charge, w->np->dielectric_in, w->np->dielectric_out,
		  w->gamma, w->sgamma,
		  w->alpha, w->salpha,
		  w->delta, w->sdelta,
		  w->ab, 
		  w->nhydrogen, w->ihydrogen, 
		  w->ndummy, w->idummy, w->isfrozen,
		  w->neigh_list, w->excl_neigh_list,
		  w->dopbc, nsym, ssize, xs, ys, zs, rot,
		  NULL) != AGBNP_OK){

      fprintf(stderr, "agbnpf_init(): error in agbnp_new().\n");
      *rescode = -1;
      return;
    }
  }

  /* turn on init flag for first energy evaluation */
  w->ener_init = 1;

  /* reset counter */
  w->e_count = 0;

  /* initialize frozen atoms constant terms */
  if(w->version == 2){
#ifdef NOTNOW
    if(agbnp2_init_frozen(w->agbnp_tag, w->x, w->y, w->z, w->isfrozen) 
       != AGBNP_OK){
      fprintf(stderr, "agbnpf_init(): error in agbnp2_init_frozen().\n");
      *rescode = -1;
      return;
    }
#endif
  }else{
    if(agbnp_init_frozen(w->agbnp_tag, w->x, w->y, w->z, w->isfrozen) 
       != AGBNP_OK){
      fprintf(stderr, "agbnpf_init(): error in agbnp_init_frozen().\n");
      *rescode = -1;
      return;
    }
  }

  if(rot){ free(rot);}
  
  *rescode = 1;
  return;
}

#ifdef NOTNOW

/* Fortran function to force re-initialization of energy/force calculation
   (triggers recreation of water site spheres)
*/
void FORTRAN(agbnpf_reinit_on)(){
  agbnp_w->ener_init = 1;
}

/* fortran interface to agbnp_init_frozen(). Re-initializes constant terms
   after a neighbor list update for example. */ 
void FORTRAN(agbnpfinitfrozen)(coord_t rx, coord_t ry, coord_t rz,
				 int *rescode){
  int iat, ispec, i;

  if(!agbnp_w){
    /* if agbnp is not initialized do nothing */
    *rescode = 1;
    return;
  }

  /* fills in coordinate arrays */
  iat = 0;
  for(ispec = 0; ispec<COMMON(impsize).nspec; ispec++){
    for(i=0; i<COMMON(impsize).natinsp[ispec]; i++, iat++){
      /* set coordinates */
      agbnp_w->x[iat] = rx[ispec][i];
      agbnp_w->y[iat] = ry[ispec][i];
      agbnp_w->z[iat] = rz[ispec][i];
    }
  }


  /* calls agbnp_init_frozen() with current list of frozen atoms */
  if(agbnp_w->version == 2){
    if(agbnp2_init_frozen(agbnp_w->agbnp_tag, 
			  agbnp_w->x, agbnp_w->y, agbnp_w->z, 
			  agbnp_w->isfrozen) != AGBNP_OK){
      fprintf(stderr, "agbnpf_init_frozen(): error in agbnp_init_frozen()\n");
      *rescode = -1;
      return;
    }
  }else{
    if(agbnp_init_frozen(agbnp_w->agbnp_tag, 
			 agbnp_w->x, agbnp_w->y, agbnp_w->z, 
			 agbnp_w->isfrozen) != AGBNP_OK){
      fprintf(stderr, "agbnpf_init_frozen(): error in agbnp_init_frozen()\n");
      *rescode = -1;
      return;
    }
  }

  /* turn on energy init flag */
  FORTRAN(agbnpf_reinit_on)();

  *rescode = 1;
  return;
}
#endif

/* Impact fortran interface with libagbnp */
void FORTRAN(agbnpf)(double *rx, double *ry, double *rz,
		     double *fx, double *fy, double *fz,
		     double *wtagbnp, 
		     double *agbnp_energy,
		     int *rescode){
  int ispec; /* species counter */
  int i, iat;   /* atom counters */
  double egb = 0.0;        /* GB energy */
  double evdw = 0.0;       /* van der Waals energy */
  double ecav = 0.0;       /* cavity energy */
  double egnp = 0.0;       /* NP energy (cavity + vdw) */
  double ecorr, ecorr_cav = 0.0, ecorr_vdw = 0.0; /* Correction energy */
  double ehb = 0.0; /* HB correction */
  AGBNPfw *w = agbnp_w;

#ifdef AGBNP_TESTDERS
  /* variables to test derivatives */
  int natoms = agbnp_w->natoms; /* total number of atoms */
  static double *x_old = NULL, *y_old=NULL, 
                *z_old=NULL;            /* previous atom coordinates */
  static double egb_old;                /* previous GB energy */
  static double evdwt_old, ecavt_old;   /* previous NP energy */
  static int iteration;
  double evdwt, ecavt, de, deder;
  static double ehbt_old;
  double ehbt;
#endif

#ifdef MMTIMER
  static MMtimer main_timer = -1, cav_timer = -1, agb_timer = -1; 
  double user_time, system_time, real_time;
#endif

  double wcav = *wtagbnp*1.0;
  double wvdw = *wtagbnp*1.0;
  double wagb = *wtagbnp*1.0;
  double whb  = *wtagbnp*1.0;

  if(!agbnp_w){
    fprintf(stderr, "agbnpf: agbnp module is not initialized. Call agbnpf_init() first.\n");
    *rescode = -1;
    return;
  }

  if(agbnp_w->version != 2){
    whb = 0.; 
  }
  
  /* fills in coordinate arrays */
  for(iat = 0; iat < w->natoms; iat++){
    /* set coordinates */
    agbnp_w->x[iat] = rx[iat];
    agbnp_w->y[iat] = ry[iat];
    agbnp_w->z[iat] = rz[iat];
  }


  if(agbnp_w->version == 2){

    /* if energy initialization has not been requested (ener_init = 0)
       decide to re-initialize based on counter */
    if(!agbnp_w->ener_init) {
      if(agbnp_w->e_count % agbnp_w->e_count_reinit == 0){
	agbnp_w->ener_init = 1;
	agbnp_w->e_count = 0;
      }
    }else{
      agbnp_w->e_count = 0;
    }
#ifdef NOTNOW
    /* in AGBNP2 there is a single energy call */
    if(agbnp2_ener(agbnp_w->agbnp_tag, agbnp_w->ener_init,
		  agbnp_w->x, agbnp_w->y, agbnp_w->z,
		  agbnp_w->sp,  agbnp_w->br, 
		  &(agbnp_w->mol_volume), agbnp_w->surf_area, 
		  &egb,              agbnp_w->dgbdr,
		  &evdw, &ecorr_vdw, agbnp_w->dvwdr, 
		  &ecav, &ecorr_cav, agbnp_w->decav,
		  &ehb,              agbnp_w->dehb) != AGBNP_OK){
      fprintf(stderr, "agbnpf: error in agbnp_ener().\n");
      *rescode = -1;
      return;
    }
#endif

    /* turn off init flag */
    agbnp_w->ener_init = 0;

  }else{
    /* in AGBNP1 there are separate calls for the electrostatic/vdw and 
       cavity components */

    /* cavity energy */
    if(agbnp_cavity_energy(agbnp_w->agbnp_tag, 
			   agbnp_w->x, agbnp_w->y, agbnp_w->z,
			   &(agbnp_w->mol_volume), agbnp_w->surf_area, 
			   &ecav, &ecorr_cav, 
			   agbnp_w->decav) != AGBNP_OK){
      fprintf(stderr, "agbnpf(): fatal error in agbnp_cavity_energy().\n");
      *rescode = -1;
      return;
    }
    
    /* agb/vdw energy */
    if(agbnp_agb_energy(agbnp_w->agbnp_tag, agbnp_w->x, agbnp_w->y, agbnp_w->z,
			agbnp_w->sp,  agbnp_w->br, &egb, agbnp_w->dgbdr,
			&evdw, agbnp_w->dvwdr, &ecorr_vdw) != AGBNP_OK){
      fprintf(stderr, "agbnpf: error in agbnp_agb_energy().\n");
      *rescode = -1;
      return;
    }

  }

#ifdef AGBNP_WRITE_STRUCT
  printf("egb = %f\n", egb);
  printf("ehb = %f\n", ehb);
  printf("evdw = %f\n",evdw);
  printf("ecorr_vdw = %f\n", ecorr_vdw);
  printf("ecav = %f\n",ecav);
  printf("ecorr_ecav = %f\n",ecorr_cav);
  printf("weights: %f %f %f %f\n",wcav,wvdw,wagb,whb);
#endif

  ecorr_cav = wcav*ecorr_cav;
  ecorr_vdw = wvdw*ecorr_vdw;
  ecav = wcav*ecav;
  evdw = wvdw*evdw;
  egb = wagb*egb;
  ehb = whb*ehb;  
  ecorr = ecorr_cav + ecorr_vdw;
  egnp = ecav + evdw;
  *agbnp_energy = egnp + egb;

    
  /* update forces */
  for(iat = 0; iat < w->natoms; iat++){
    fx[iat] -= wagb*w->dgbdr[iat][0] + wvdw*w->dvwdr[iat][0] + wcav*w->decav[iat][0] + whb*w->dehb[iat][0];
    fy[iat] -= wagb*w->dgbdr[iat][1] + wvdw*w->dvwdr[iat][1] + wcav*w->decav[iat][1] + whb*w->dehb[iat][1];
    fz[iat] -= wagb*w->dgbdr[iat][2] + wvdw*w->dvwdr[iat][2] + wcav*w->decav[iat][2] + whb*w->dehb[iat][2];
  }

#ifdef AGBNP_WRITE_STRUCT
  /* write coordinates etc and results to stderr */
  {
    int j;
    char atom_type[5];
    double all_area;
    

    all_area = 0.0;
    printf("        Atom    ntype stype         Area            x           y         z       Radius   Scaled_Volume_Factor  Born_r   Charge\n");  
    
    for(iat = 0; iat<w->natoms; iat++){
      printf("#agbnp_log %4d    %5d %4s  %12.3f A^2  %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	     iat+1,w->typeid[iat], w->np->symbol[w->typeid[iat]],agbnp_w->surf_area[iat],
	     agbnp_w->x[iat],agbnp_w->y[iat],agbnp_w->z[iat], 
	     agbnp_w->r[iat],    
	     agbnp_w->sp[iat], agbnp_w->br[iat], agbnp_w->charge[iat]);
      all_area += agbnp_w->surf_area[iat];
    }
    printf("Total surface area = %12.3f A^2\n", all_area);

  }
#endif

#ifdef AGBNP_TESTDERS
  if(x_old){
    iteration += 1;
    evdwt = evdw + ecorr_vdw;
    ecavt = ecav + ecorr_cav;
    ehbt = ehb; 
    de = (egb-egb_old) + (evdwt - evdwt_old) + 
      (ecavt - ecavt_old) + (ehbt - ehbt_old);
    deder = 0.0;
    for(iat=0;iat<natoms;iat++){
      deder += 
	(wagb*w->dgbdr[iat][0] + wvdw*w->dvwdr[iat][0] + 
	 wcav*w->decav[iat][0] + whb*w->dehb[iat][0])*
	  (w->x[iat]-x_old[iat]) +
	(wagb*w->dgbdr[iat][1] + wvdw*w->dvwdr[iat][1] + 
	 wcav*w->decav[iat][1] + whb*w->dehb[iat][1])*
	  (w->y[iat]-y_old[iat]) +
	(wagb*w->dgbdr[iat][2] + wvdw*w->dvwdr[iat][2] + 
	 wcav*w->decav[iat][2] + whb*w->dehb[iat][2])*
	  (w->z[iat]-z_old[iat]);
    }
    fprintf(stderr,"agbnp: %d %e %e %e\n", iteration, 
	  egb + evdwt + ecavt + ehbt, de, deder);
        printf("agbnp: %d %e %e %e\n", iteration, 
	  egb + evdwt + ecavt + ehbt, de,  deder); 
    egb_old = egb;
    evdwt_old = evdwt;
    ecavt_old = ecavt;
    ehbt_old = ehbt;
    for(iat=0;iat<natoms;iat++){
      x_old[iat] = w->x[iat];
      y_old[iat] = w->y[iat];
      z_old[iat] = w->z[iat];
    }
  }else{
    /* first time through */
    x_old = (double *)realloc(x_old, natoms*sizeof(double));
    y_old = (double *)realloc(y_old, natoms*sizeof(double));
    z_old = (double *)realloc(z_old, natoms*sizeof(double));
    iteration = 0;
    evdwt = evdw + ecorr_vdw;
    ecavt = ecav + ecorr_cav;
    ehbt = ehb;
    egb_old = egb;
    evdwt_old = evdwt;
    ecavt_old = ecavt;
    ehbt_old = ehbt;
    for(iat=0;iat<natoms;iat++){
      x_old[iat] = w->x[iat];
      y_old[iat] = w->y[iat];
      z_old[iat] = w->z[iat];
    }
  }
#endif

  /* update counter */
  w->e_count += 1;

  /* done */
  *rescode = 1;
  return;
}

