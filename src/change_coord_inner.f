      SUBROUTINE change_coord_inner(ntap,nprot,protl,ngrp,grppt,co,oc
     &     ,pmass,lx,ly,lz,xp0,yp0,zp0,xpa,ypa,zpa,xpcm,ypcm,zpcm,xpcma
     &     ,ypcma,zpcma)

************************************************************************
*   Time-stamp: <98/03/03 12:21:40 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Wed Feb 11 1998 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER ntap,nprot,ngrp,protl(*),grppt(2,*)
      REAL*8  co(3,3),oc(3,3),pmass(*),lx(*),ly(*),lz(*),xp0(*),yp0(*)
     &     ,zp0(*),xpa(*),ypa(*),zpa(*),xpcm(*),ypcm(*),zpcm(*),xpcma(*)
     &     ,ypcma(*),zpcma(*)

*------------------------- LOCAL VARIABLES ----------------------------*


*----------------------- EXECUTABLE STATEMENTS ------------------------*


*--- After advancing lx,ly,lz recompute xp0,yp0,zp0

      CALL change_origin(-1,nprot,protl,xp0,yp0,zp0,lx,ly
     &     ,lz,xpcma,ypcma,zpcma,co)
                  
*---        go from real particle coordinates to scaled ----------------
      
      CALL change_frame(co,oc,-1,ntap,xp0,yp0,zp0,xpa,ypa
     &     ,zpa)
      CALL change_frame(co,oc,1,nprot,xpcma,ypcma,zpcma,xpcm
     &     ,ypcm,zpcm)
      
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
