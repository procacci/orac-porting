      SUBROUTINE brot(read_co,a,b,c,alfa,beta,gamma,i,j,k,co,oc,volume)

************************************************************************
*                                                                      *
*                                                                      *
*     This subroutine computes the rotation matrix that transforms     *
*     simulation box coordinates into orthogonalized coordinates       *
*     and its inverse.                                                 *
*                                                                      *
*     A       :  Unit cell length along the x direction.               *
*     B       :  Unit cell length along the y direction.               *
*     C       :  Unit cell length along the z direction.               *
*     ALFA    :  Angle between the y and z axis.                       *
*     BETA    :  Angle between the z and x axis.                       *
*     GAMMA   :  Angle between the x and y axis.                       *
*     I       :  Number of cell along x.                               *
*     J       :  Number of cell along y.                               *
*     K       :  Number of cell along z.                               *
*     CO      :  Rotation matrix from simulation box frame to          *
*                orthogonalized frame.                                 *
*                >> real*8 CO(3,3) <<                                  *
*     OC      :  Inverse of CO.                                        *
*                >> real*8 OC(3,3) <<                                  *
*     VOLUME  :  Volume of the simulation box cell in                  *
*                Angstroem**3.                                         *
*                                                                      *
*----- Last update 04/13/89 -------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*                                                                      *
*     EXTERNALS  MATINV, NEAR0.                                        *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER i,j,k
      REAL*8 a,b,c,alfa,beta,gamma,volume,co(3,3),oc(3,3)
      LOGICAL read_co

*-------------------- EXTERNAL FUNCTION --------------------------------

      EXTERNAL near0
      LOGICAL near0

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER m,n
      REAL*8  degrad,qt,alf,bet,gam,ax,bx,by,cx,cy,cz

*==================== EXECUTABLE STATEMENTS ============================

      IF(.NOT. read_co) THEN
         co(2,1)=0.0D0
         co(3,1)=0.0D0
         co(3,2)=0.0D0
         degrad=pi/180.0d0
         ax=a
         alf=DCOS(degrad*alfa)
         bet=DCOS(degrad*beta)
         qt=DSIN(degrad*gamma)
         gam=DCOS(degrad*gamma)
         bx=b*gam
         by=b*qt
         cx=c*bet
         cy=c*(alf-bet*gam)/qt
         cz=dsqrt(c*c-cx*cx-cy*cy)
         co(1,1)=ax
         co(1,2)=bx
         co(1,3)=cx
         co(2,2)=by
         co(2,3)=cy
         co(3,3)=cz
      END IF
      DO 10 m=1,3
          co(m,1)=DFLOAT(i)*co(m,1)/boxl
          co(m,2)=DFLOAT(j)*co(m,2)/boxl
          co(m,3)=DFLOAT(k)*co(m,3)/boxl
10    CONTINUE
      CALL matinv(3,3,co,oc,volume)
      DO 20 m=1,3
          DO 30 n=1,3
              IF(near0(co(m,n))) co(m,n)=0.0D+0
              IF(near0(oc(m,n))) oc(m,n)=0.0D+0
30        CONTINUE
20    CONTINUE
      volume=volume*boxl**3

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
