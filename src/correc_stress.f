      SUBROUTINE correc_stress(cpress,kinetic,nprot,co,oc,vco,vcax,vcay
     &     ,vcaz,fcax,fcay,fcaz,stressd,stressr,volume,press,press_kin
     &     ,pext,tmass,masspp,tm,tm2)

************************************************************************
*   Time-stamp: <98/02/13 21:48:51 marchi>                             *
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

      use unit

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      LOGICAL cpress,kinetic
      INTEGER nprot
      REAL*8  tm,tm2,volume,press,press_kin,pext,masspp(*),tmass(*),co(3
     &     ,3),oc(3,3),stressd(3,3),stressr(3,3),vco(3,3),vcax(*),vcay(
     &     *),vcaz(*),fcax(*),fcay(*),fcaz(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  prt(3,3),st(3,3)
      INTEGER i,j

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      CALL correc(vcax,vcay,vcaz,fcax,fcay,fcaz,tmass,nprot,tm2)

      IF(cpress) THEN
c-------------- Compute total stress -----------------------------------
         
         CALL comp_stress_conf(stressd,stressr,prt,oc,volume
     &        ,unitp,press)

         IF(kinetic) THEN
            CALL comp_stress_kinetic(vcax,vcay,vcaz,tmass,co,nprot
     &           ,volume,unitp,st,press_kin)
            
            CALL comp_forcep(prt,st,oc,volume,pext)
         ELSE
            CALL DGEMM('N','T',3,3,3,1.0D0,prt,3,oc,3,0.0D0,st,3)
            CALL dcopy(9,st,1,prt,1)
         END IF

c-------------- Correct velocities of the cofm by the total stress -----
         
         CALL correc3x3(vco,prt,masspp,tm)
      END IF
      CALL correc(vcax,vcay,vcaz,fcax,fcay,fcaz,tmass,nprot,tm2)
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
