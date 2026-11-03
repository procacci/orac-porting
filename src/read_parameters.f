      SUBROUTINE read_parameters(iret,errmsg)

************************************************************************
*   Time-stamp: <97/07/14 11:52:04 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Nov 19 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use parst

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER iret
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN COMMON --------------------------*

      INCLUDE 'parameters.h'
c-
C$$ Awful patchwork to avoid nasty bugs in rdaprm 
      INTEGER nmax_bond
      INTEGER nmax_bend
      INTEGER nmax_ptor
      INTEGER nmax_itor
      COMMON /dim_for_prm/ nmax_bond,nmax_bend,nmax_ptor,nmax_itor
      nmax_bond = bondp 
      nmax_bend = bendp 
      nmax_ptor = torsp
      nmax_itor = itorp
C$$ END of Awful patchwork to avoid nasty bugs in rdaprm 
c-

*----------------------- EXECUTABLE STATEMENTS ------------------------*

!=======================================================================
!------------ Read in the macromolecular topology and the --------------
!------------ model potential parameters -------------------------------
!=======================================================================

      xnbond=0
      xntor=0
      xnitor=0
      CALL crdrtf(alpha,natype,qchge,natop,jngrp,jgrppt,jbnd,nmbo,jbend
     &     ,nbend,jtor,ntor,jitor,nitor,NOATO,cnat,alphb,BACKP,nback
     &     ,jacc,jdon,nacc,ndon,HYDP,albon,pbon1,pbon2,lpbon,albnd,pbnd1
     &     ,pbnd2,pbnd3,pbnd4,pbnd5,pbnd6,pbnd7,lpbnd,altor,ptor1,ntor2
     &     ,potj,TORSP,lptor,alito,pito1,pito2,pito3,lpito,massu,alhyd
     &     ,xbond,xnbond,rbond,xtor,xntor,rtor,xitor,xnitor,ritor,stpr
     &     ,stsc,nstpr,nstsc,nstrg,NORSRE,NOATRE,iret,errmsg,qmss)

!----------------- If iret .EQ. 1 STOP !! ------------------------------
!                                 ====

      IF(lpbon .GT. bondp) THEN
          WRITE(errmsg,'(''In JOIN: Dimensions of bond parameters'',
     x   '' exceeded. '',i6,''>'',i6,'' Abort.'')') lpbon,bondp
          iret=1
      END IF

      IF(lpbnd .GT. bendp) THEN
          errmsg='In JOIN: Dimensions of angle parameters exceeded.'//
     x            ' Abort.'
          iret=1
      END IF
      IF(lptor .GT. torsp) THEN
          errmsg='In JOIN: Dimensions of torsional parameters '//
     x            'exceeded. Abort.'
          iret=1
      END IF
      IF(lpito .GT. itorp) THEN
          WRITE(errmsg,'(''In JOIN: Dimensions of i-torsional'',
     x    '' parameters exceeded. LPITO ='',i6)') lpito
          iret=1
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END

