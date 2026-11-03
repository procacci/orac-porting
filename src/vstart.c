#include <math.h>
#define   TRUE    1
#define   FALSE   0
 
#define   OK      1
#define   FAIL    0
 
#include "config.H"
 
#define   X       0
#define   Y       1
#define   Z       2
#define   D       3
 
#define   EPS     1e-5
#define   MAX     1e9
#define   CUTTER  0.99
 
#define   det(A,B,C)  (idata[n1][A]*((idata[n3][B]*idata[i][C])-\
                                     (idata[n3][C]*idata[i][B])))+\
                      (idata[n1][B]*((idata[n3][C]*idata[i][A])-\
                                     (idata[n3][A]*idata[i][C])))+\
                      (idata[n1][C]*((idata[n3][A]*idata[i][B])-\
                                     (idata[n3][B]*idata[i][A])))
 
#define   abs(X)      (fabs((double)X)) 
#define   sq(X)       ((X)*(X))
#define   dist(A,B,C) (sq(A)+sq(B)+sq(C))
#define   deg(A)      (idata[MAXPLA-*ndeg-1][A])

 
double fabs();

/**********************************************************************
*
*  VORONOI computes a face of a Voronoi polyhedron.
*
*  Inputs to the subroutine are:
*
*  idata     set of planes on which all the faces of the polyhedron lay
*
*  x,y,z     position of one vertex of the face from which we start
*            the computation of all other vertices
*
*  n1,n2,n3  the indices of the three planes in the array idata
*            intersecting at point
*
*  ndeg      number of degenerate vertex found
*
*  Output of the subroutine are:
*
*  odata     array containing all the vertices computed sorted by face
*
*  vertex    is the number of vertices for each face
*
***********************************************************************/
 
int   voronc(x,y,z,n1,n2,n3,idata,odata,vertex,ndeg)
int   n1, n2, n3, *ndeg;
long int vertex[MAXPLA];
double x, y, z, idata[MAXPLA][4], odata[MAXPLA][MAXVER][3];
 
{
int    npla1[MAXVER], npla2[MAXVER], npla3[MAXVER];
int    i, j, km, lgo, code;
double  rmin, fmix, abc, dbc, adc, abd, r, regn;
 
/*----------------------------------------------------------------------
*       check if the face has been already computed
*---------------------------------------------------------------------*/

if (!vertex[n1])
 
/*----------------------------------------------------------------------
*       find all vertices over face N1
*---------------------------------------------------------------------*/
 
    {
    km = -1;
    npla1[vertex[n1]] = n1;
    npla2[vertex[n1]] = n2;
    npla3[vertex[n1]] = n3;
    odata[n1][vertex[n1]][X] = x;
    odata[n1][vertex[n1]][Y] = y;
    odata[n1][vertex[n1]][Z] = z;
 
    while (km != npla2[0])
         {
         rmin = MAX;
         vertex[n1] ++;
         for (i=0 ; i<MAXPLA ; i++)
              {
              if (i == n1 || (i == npla2[0] && vertex[n1] < 2))
                   continue;
              lgo = TRUE;
              for (j=0 ; j<vertex[n1] ; j++)
                   {
                   if (i == npla3[j])
                        {
                        lgo = FALSE;
                        break;
                        }
                   }
              if (lgo)
                   {
                   fmix = det(X,Y,Z);
                   if (abs(fmix) >= EPS)
                        {
                        abc  = 1/fmix;
                        dbc  = det(D,Y,Z);
                        adc  = det(X,D,Z);
                        abd  = det(X,Y,D);
 
                        x = dbc * abc;
                        y = adc * abc;
                        z = abd * abc;
 
                        r = dist(x-odata[n1][vertex[n1]-1][X],
                                 y-odata[n1][vertex[n1]-1][Y],
                                 z-odata[n1][vertex[n1]-1][Z]);
 
                        regn = x * idata[n2][X] +
                               y * idata[n2][Y] +
                               z * idata[n2][Z] -
                                   idata[n2][D];
 
                        if (abs(r) < EPS && abs(regn) < EPS)
                             {
                             deg(X) = x*CUTTER;
                             deg(Y) = y*CUTTER;
                             deg(Z) = z*CUTTER;
                             deg(D) = dist(deg(X),deg(Y),deg(Z));
                             (*ndeg)++;
                             for (j=0 ; j<MAXPLA ; j++)
                                  vertex[j] = 0;
                             return(FAIL);
                             }
                        if (r < rmin && regn < 0)
                             {
                             rmin = r;
                             odata[n1][vertex[n1]][X] = x;
                             odata[n1][vertex[n1]][Y] = y;
                             odata[n1][vertex[n1]][Z] = z;
                             km = i;
                             }
                        }
                   }
              }
 
         n2 = n3;
         n3 = km;
         npla1[vertex[n1]] = n1;
         npla2[vertex[n1]] = n2;
         npla3[vertex[n1]] = n3;
         }
 
    vertex[n1] ++;
 
    for (i=0 ; i<vertex[n1] ; i++)
         if ((code = voronc(odata[n1][i][X],
                             odata[n1][i][Y],
                             odata[n1][i][Z],
                             npla2[i],
                             npla1[i],
                             npla3[i],
                             idata,
                             odata,
                             vertex,
                             ndeg)) == FAIL) return(FAIL);
    }
return(OK);
}
 
/**********************************************************************
*
*  VSTART computes the starting point of a face of a Voronoi polyhedron.
*
*  Inputs to the subroutine are:
*
*  idata     set of planes on which all the faces of the polyhedron lay
*
*  Output of the subroutine are:
*
*  odata     array containing all the vertices computed sorted by face
*
*  vertex    is the number of vertices for each face
*
*  The subroutine calls VORONOI after computing the starting vertex.
*
***********************************************************************/



/* compile this file with -DFC_UNDERSCORES if the fortran compiler
   adds a trailing underscore to symbol names
   (e.g. gfortran, ifort )
*/

#ifdef FC_UNDERSCORES
void  vstart_(idata,odata,vertex)
#else
void  vstart(idata,odata,vertex)
#endif

long int   vertex[MAXPLA];
double idata[MAXPLA][4], odata[MAXPLA][MAXVER][3];
 
{
int   i, j, k, h, kk, k3, code, ndegen, sum;
double rmin, ab, aabb, bbab, aaab, ccc, dete, alfa, beta, c1o, c2o, c3o,
      cco, x2, y2, z2, rx, ry, rz, r, cr1, cr2, cr3, ccr, xr2, yr2,
      zr2, x11, x22, xlambda, p1n, p2n, p3n, restore;
struct
    {
    int count;
    int normver[2];
    } connect[MAXPLA];
 
/*----------------------------------------------------------------------
*   Procedure for constructing voronoi poliedra starts here
*   for I=1 is the atom is at the minimum distance from the
*   central atom; the starting informations consists in finding
*   The first vertex and the first two edges of the poliedrum
*---------------------------------------------------------------------*/
 
for (i=0 ; i<MAXPLA ; i++)
    vertex[i] = 0;
 
code = FAIL;
ndegen = 0;
 
/*   -----  find intersection line ------    */
 
while (code == FAIL)
    {
     rmin = MAX;
     for (i=1 ; i<MAXPLA ; i++)
         {
         ab = idata[0][X]*idata[i][X] +
              idata[0][Y]*idata[i][Y] +
              idata[0][Z]*idata[i][Z];
         aabb = idata[0][D]*idata[i][D];
         bbab = idata[i][D]*ab;
         aaab = idata[0][D]*ab;
 
         ccc = aabb - ab*ab;
 
         if (ccc > EPS)
              {
              dete = 1/ccc;
              alfa = (aabb-bbab)*dete;
              beta = (aabb-aaab)*dete;
              c1o  = idata[0][Y]*idata[i][Z] -
                     idata[0][Z]*idata[i][Y];
              c2o  = idata[0][Z]*idata[i][X] -
                     idata[0][X]*idata[i][Z];
              c3o  = idata[0][X]*idata[i][Y] -
                     idata[0][Y]*idata[i][X];
              cco  = dist(c1o,c2o,c3o);
              x2   = alfa*idata[0][X]+beta*idata[i][X];
              y2   = alfa*idata[0][Y]+beta*idata[i][Y];
              z2   = alfa*idata[0][Z]+beta*idata[i][Z];
              rx   = x2-idata[0][X];
              ry   = y2-idata[0][Y];
              rz   = z2-idata[0][Z];
              r    = dist(rx,ry,rz);
              if (r < rmin)
                   {
                   rmin = r;
                   kk   = i;
                   cr1  = c1o;
                   cr2  = c2o;
                   cr3  = c3o;
                   ccr  = cco;
                   xr2  = x2;
                   yr2  = y2;
                   zr2  = z2;
                   }
              }
         }
 
/*   -----  find first vertex ------    */
 
     rmin = MAX;
     for (i=1 ; i<MAXPLA ; i++)
         {
         if (i != kk)
              {
              x11 = idata[i][D] - (xr2*idata[i][X] +
                                    yr2*idata[i][Y] +
                                    zr2*idata[i][Z]);
              x22 = cr1*idata[i][X] +
                    cr2*idata[i][Y] +
                    cr3*idata[i][Z];
              if (abs(x22) < EPS) continue;
              xlambda = x11/x22;
              r      = xlambda*xlambda*ccr;
              if (r < rmin)
                   {
                   k3   = i;
                   rmin = r;
                   p1n  = xr2 + xlambda*cr1;
                   p2n  = yr2 + xlambda*cr2;
                   p3n  = zr2 + xlambda*cr3;
                   }
              }
         }
         code = voronc(p1n,p2n,p3n,0,kk,k3,idata,odata,vertex,&ndegen);
    }
if (ndegen)
    {
    restore = 1/CUTTER;
    for (i=MAXPLA-1 ; i>MAXPLA-1-ndegen ; i--)
         {
         x2 = idata[i][X]*restore;
         y2 = idata[i][Y]*restore;
         z2 = idata[i][Z]*restore;
         for (j=0 ; j<MAXPLA ; j++)
              connect[j].count = 0;
         for (j=0 ; j<vertex[i] ; j++)
              {
              for (k=0 ; k<MAXPLA-ndegen ; k++)
                   {
                   for (h=0 ; h<vertex[k] ; h++)
                        {
                        if ((abs(odata[k][h][X] - odata[i][j][X]) <= EPS) &&
                            (abs(odata[k][h][Y] - odata[i][j][Y]) <= EPS) &&
                            (abs(odata[k][h][Z] - odata[i][j][Z]) <= EPS))
                             {
                             connect[k].normver[connect[k].count] = h;
                             connect[k].count++;
                             }
                        }
                   }
              }
         for (j=0 ; j<MAXPLA ; j++)
              if (connect[j].count != 0)
                   {
                   if(abs(connect[j].normver[0]-connect[j].normver[1]) ==
                         (vertex[j]-1))
                        {
                        odata[j][0][X] = x2;
                        odata[j][0][Y] = y2;
                        odata[j][0][Z] = z2;
                        }
                   else
                        {
                        sum = connect[j].normver[0] + connect[j].normver[1];
                        odata[j][(sum-1)/2][X] = x2;
                        odata[j][(sum-1)/2][Y] = y2;
                        odata[j][(sum-1)/2][Z] = z2;
                        for (k=(sum+1)/2 ; k<vertex[j]-1 ; k++)
                             {
                             odata[j][k][X] = odata[j][k+1][X];
                             odata[j][k][Y] = odata[j][k+1][Y];
                             odata[j][k][Z] = odata[j][k+1][Z];
                             }
                        }
                   vertex[j]--;
                   }
         vertex[i] = 0;
         }
    }
}


