      SUBROUTINE print_cm(x,y,z,mass,ntap,massinfty)

************************************************************************
*   Time-stamp: <98/02/10 16:09:11 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Feb 10 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  x(*),y(*),z(*),mass(*),massinfty
      INTEGER ntap

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i
      REAL*8  totmass,totmassm,xtot,ytot,ztot,xm,ym,zm,xx,yy,zz

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      xtot=0.0D0
      ytot=0.0D0
      ztot=0.0D0
      xm=0.0D0
      ym=0.0D0
      zm=0.0D0
      totmass=0.0D0
      totmassm=0.0D0
      DO i=1,ntap
        if(mass(i).lt.massinfty) THEN
          xtot=xtot+x(i)*mass(i)
          ytot=ytot+y(i)*mass(i)
          ztot=ztot+z(i)*mass(i)
          totmass=totmass+mass(i)
        ELSE
          xm=xm+x(i)*mass(i)
          ym=ym+y(i)*mass(i)
          zm=zm+z(i)*mass(i)
          totmassm=totmassm+mass(i)
        END IF
      END DO
      xx=xtot+xm
      yy=ytot+ym
      zz=ztot+zm
      write(6,10) xtot/totmass,ytot/totmass,ztot/totmass
      if(totmassm.gt.1.D-15) THEN 
        write(6,11) xm/totmassm,ym/totmassm,zm
     &       /totmassm
        totmass=totmass+totmassm
        write(6,12) xx/totmass,yy/totmass,zz/totmass
      END IF
10    FORMAT(' COM of light atoms ', 3f15.8)
11    FORMAT(' COM of heavy atoms ', 3f15.8)
12    FORMAT(' COM of all atoms   ', 3f15.8)

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
