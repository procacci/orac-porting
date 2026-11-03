      SUBROUTINE zeroa(ax,ay,az,n,m)

************************************************************************
*                                                                      *
*                                                                      *
*     This subroutine set the elements of 3 arrays to zero.            *
*                                                                      *
*     AX      :  Arrays whose element must be set to zero.             *
*     AY         >> real*8 AX(N*M), AY(N*M), AZ(N*M) <<                *
*     AZ                                                               *
*                                                                      *
*     N       :  Physical row dimension of A.                          *
*     M       :  Physical column dimension of A.                       *
*                                                                      *
*---- Last update 04/23/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS NONE                                                   *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER n,m
      REAL*8 ax(n*m),ay(n*m),az(n*m)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,l

*==================== EXECUTABLE STATEMENTS ============================

      l=0
      DO 10 i=1,m
          DO 20 j=1,n
              ax(j+l)=0.0D+0
              ay(j+l)=0.0D+0
              az(j+l)=0.0D+0
20        CONTINUE
          l=l+n
10    CONTINUE

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
