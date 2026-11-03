      SUBROUTINE plotc ( fstep,beta,co,xp0,yp0,zp0,ntap,nres,m1,prsymb
     &     ,para0)

************************************************************************
*                                                                      *
*     Print a snapshot of the solute molecules on unit kplot.          *
*                                                                      *
*     CO      :  Transform the coordinates to the orthogonalized  (I)  *
*                crystallographic frame.                               *
*                >> real*8  CO(3,3) <<                                 *
*     BETA    :  List of labels for each site of the solute.      (I)  *
*                >> character*7 BETA(NTAP) <<                          *
*     XP0     :  Coordinates of the solute molecules, packed by   (I)  *
*     YP0        site.                                                 *
*     ZP0        >> real*8  XP0(NTAP), ... <<                          *
*                                                                      *
*     NTAP    :  Number of sites on the solute molecules.         (I)  *
*                                                                      *
*---- Last update 05/03/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS NONE.                                                  *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit
      use rem, only: para_index, rem_run
      use sge, only: sim_tempering

      implicit none

*----------------------- ARGUMENTS -------------------------------------

      integer :: ntap,m1
      integer :: nres(m1,*)
      integer, intent(in)  :: para0
      real(8), intent(in) :: xp0(*), yp0(*), zp0(*), co(3,3)
      real(8), intent(in) :: fstep
      character(7) :: beta(*)
      character(8) :: prsymb(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      integer :: i, j, l
      real(8) :: xb, yb, zb
      character(5) :: bet, bet2
      character(3) :: rsd

*==================== EXECUTABLE STATEMENTS ============================

#ifdef _MPI_
      if ( rem_run .or. sim_tempering ) then
#else
      if ( sim_tempering ) then
#endif
        write(kout,2) fstep,
     &       2.d0*co(1,1), 2.d0*co(1,2), 2.d0*co(1,3),
     &       2.d0*co(2,1), 2.d0*co(2,2), 2.d0*co(2,3),
     &       2.d0*co(3,1), 2.d0*co(3,2), 2.d0*co(3,3),
     &       para_index,para0
        
2       format('REMARK   Time:', 1x, f15.2, 1x, 'CO matrix:', 1x, 
     &       9f10.5, 1x,'Ensemble:', 1x, 2i4)
      else
        write (kout,3) fstep,
     &       2.d0*co(1,1), 2.d0*co(1,2), 2.d0*co(1,3),
     &       2.d0*co(2,1), 2.d0*co(2,2), 2.d0*co(2,3),
     &       2.d0*co(3,1), 2.d0*co(3,2), 2.d0*co(3,3)
3       format('REMARK   Time:', 1x, f15.2, 1x, 'CO matrix:', 1x,
     &       9f10.5)
      endif

      DO i = 1, ntap
        j = nres(i,2)
        xb = xp0(i)
        yb = yp0(i)
        zb = zp0(i)
        bet(1:5) = beta(i)(1:5)
        CALL low_up(bet,5)
        DO l = 1, 5
          bet2(l:l)=bet(l:l)
        ENDDO
        rsd(1:3)=prsymb(j)(1:3)
        CALL low_up(rsd,3)
        WRITE(kout,1)'ATOM ',i,bet2(1:5),rsd,nres(i,1),xb,yb,zb
      ENDDO

      WRITE(kout,'(a)') 'END'
1     FORMAT(a5,i6,1x,a5,a3,1x,i5,4x,3f8.3)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
