      SUBROUTINE readco(co,oc,volume)

************************************************************************
*   Time-stamp: <95/07/02 11:24:21 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Sun Jul  2 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit
      use rem
      use sge, only: sge_update_rst, sim_tempering

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  co(3,3),oc(3,3),volume

*------------------------- LOCAL VARIABLES ----------------------------*
      
      INTEGER i,j
      REAL*8  oc1(3,3)
      
*----------------------- EXECUTABLE STATEMENTS ------------------------*
      
      REWIND kdump
#ifdef _MPI_
      if(rem_run) then
        read(kdump,*)
      elseif(sim_tempering) then
#else
      if(sim_tempering) then
#endif
        if( sge_update_rst > 0 ) then
          do i = 1, 6 + 2*sge_update_rst
            read(kdump,*)
          enddo
        else
          do i = 1, 5
            read(kdump,*)
          enddo
        endif
      endif

      READ(kdump,*) ((co(i,j),j=1,3),i=1,3),((oc(i,j),j=1,3),i=1,3)
      CALL matinv(3,3,co,oc1,volume)
      volume=volume*boxl**3

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
