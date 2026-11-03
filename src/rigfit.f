      SUBROUTINE rigfit(puretr,n,xyz,xyz0,wt,wr,q,dcm,xyzfit,
     .                  error,iret,msg)
 
**********************************************************************
*                                                                    *
*     SR RIGFIT fits a rigid structure of N sites into a given       *
*     structure of N sites.                                          *
*                                                                    *
*     AUTHOR: Gerald Kneller                                         *
*             IBM France                                             *
*             94-96, rue Reaumur                                     *
*             F-75002 Paris, FRANCE                                  *
*                                                                    *
*---- LAST UPDATE: 26.07.1990 ---------------------------------------*
*                                                                    *
*     ARGUMENTS:                                                     *
*                                                                    *
*     PURETR  : Parameter to determine type of the fit:          (I) *
*               PURETR = 0: The fit includes a translation and       *
*                           a rotation.                              *
*               PURETR = 1: A pure translational fit is done.        *
*               PURETR = 2: A pure rotational fit is done.           *
*               >> INTEGER PURETR <<                                 *
*     N      : Number of sites.                                 (I)  *
*              >> INTEGER N <<                                       *
*     XYZ    : Input coordinates.                               (I)  *
*              >> REAL*8 XYZ(3,*) <<                                 *
*     XYZ0   : Definition of the rigid structure.               (I)  *
*              >> REAL*8 XYZ0(3,*) <<                                *
*     WT     : Relative weights of the sites for the CM/CR.     (I)  *
*              >> REAL*8 WT(*) <<                                    *
*     WR     : Relative weights of the sites for the rot. fit.  (I)  *
*              >> REAL*8 WR(*) <<                                    *
*     Q      : Fitted quaternions.                              (O)  *
*              >> REAL*8 Q(0:7) <<                                   *
*     DCM    : CM/CR-difference between input and reference          *
*              structure.                                       (O)  *
*              >> REAL*8 DCM(3) <<                                   *
*     XYZFIT : Fitted coordinates.                              (O)  *
*              >> REAL*8 XYZFIT(3,*) <<                              *
*     ERROR  : RMS error / site.                                (O)  *
*              >> REAL*8 ERROR <<                                    *
*     IRET   : Return-code identical with LSQQTN.               (O)  *
*              >> INTEGER iret <<                                    *
*     MSG    : Message corresponding to IRET.                   (O)  *
*              >> CHARACTER*80 msg <<                                *
*                                                                    *
*     EXTERNALS:                                                     *
*                                                                    *
*     - subroutine LSQQTN -> sr. EIGRS                               *
*                           (or equiv.)                              *
*                                                                    *
**********************************************************************
 
*==== DECLARATIONS: =================================================*
 
      IMPLICIT NONE
 
*---- ARGUMENTS: ----------------------------------------------------*
 
      integer      puretr,n,iopt,iret
      real*8       xyz(3,*),xyz0(3,*),xyzfit(3,*),wt(*),wr(*),
     .             error,q(0:7),dcm(3)
      CHARACTER*80 MSG
 
*---- LOCAL CONSTANTS: ----------------------------------------------*
 
      integer   nmax
      parameter (nmax=10000)
 
*---- LOCAL VARIABLES: ----------------------------------------------*
 
      integer          iret2,i
      REAL*8 xs,ys,zs,x0s,y0s,z0s,
     .                 r(3,nmax),r0(3,nmax),rfit(3,nmax)
      character*80     msg2
 
      equivalence (rfit,r)
 
*==== EXECUTABLE STATEMENTS: ========================================*
 
*---- CHECKS: -------------------------------------------------------*
 
      if (puretr .ne. 0 .and. puretr .ne. 1 .and. puretr .ne. 2)
     .   goto 8900
      if (n .gt. nmax) goto 9000
      if (n .le. 1)    goto 9010
      iret = 0
 
*---- REFER STRUCTURES TO THE CENTER OF MASS/ROTATION ---------------*
 
*.....INPUT STRUCTURE:
      xs = 0.d0
      ys = 0.d0
      zs = 0.d0
      do 200 i = 1,n
         if (wt(i) .eq. 0.d0) goto 200
         xs = xs + wt(i)*xyz(1,i)
         ys = ys + wt(i)*xyz(2,i)
         zs = zs + wt(i)*xyz(3,i)
200   continue
      do 205 i = 1,n
         r(1,i) = xyz(1,i) - xs
         r(2,i) = xyz(2,i) - ys
         r(3,i) = xyz(3,i) - zs
205   continue
 
*.....REFERENCE STRUCTURE:
      x0s = 0.d0
      y0s = 0.d0
      z0s = 0.d0
      do 210 i = 1,n
         if (wt(i) .eq. 0.d0) goto 210
         x0s = x0s + wt(i)*xyz0(1,i)
         y0s = y0s + wt(i)*xyz0(2,i)
         z0s = z0s + wt(i)*xyz0(3,i)
210   continue
      do 215 i = 1,n
         r0(1,i) = xyz0(1,i) - x0s
         r0(2,i) = xyz0(2,i) - y0s
         r0(3,i) = xyz0(3,i) - z0s
215   continue
 
*---- DO ROTATIONAL FIT AND SHIFTING OF THE CONTROIDS ACCORDING -----*
*---- TO THE VALUE OF "PURETR": -------------------------------------*
 
*.....A TRANSLATIONAL AND ROTATIONAL FIT HAS TO BE PERFORMED:
 
      if (puretr .eq. 0) then
*........DO THE ROTATIONAL FIT:
         iopt = 1
         call LSQQTN(n,r,r0,wr,q,rfit,error,iopt,iret2,msg2)
         if (iret2 .eq. 99) goto 9020
         iret = iret2
         write(msg,'(a11,a69)')' RIGFIT -> ',msg2
*........CALCULATE THE POSITION DIFFERENCE OF THE CMs/CRs:
         dcm(1) = xs - x0s
         dcm(2) = ys - y0s
         dcm(3) = zs - z0s
*....... SHIFT FITTED CM/CR-COORDINATES:
         do 300 i = 1,n
            xyzfit(1,i) = rfit(1,i) + xs
            xyzfit(2,i) = rfit(2,i) + ys
            xyzfit(3,i) = rfit(3,i) + zs
300      continue
 
      endif
 
*.....A PURE TRANSLATIONAL FIT HAS TO BE PERFORMED:
 
      if (puretr .eq. 1) then
*........SET QUATERNION(S) FOR ZERO ANGLE ROTATION:
         q(0) = 1.0d0
         q(1) = 0.0d0
         q(2) = 0.0d0
         q(3) = 0.0d0
         q(4) = 0.0d0
         q(5) = 0.0d0
         q(6) = 0.0d0
         q(7) = 0.0d0
*........CALCULATE THE POSITION DIFFERENCE OF THE CMs/CRs:
         dcm(1) = xs - x0s
         dcm(2) = ys - y0s
         dcm(3) = zs - z0s
*....... SHIFT UNROTATED CM/CR-COORDINATES OF THE REFERENCE STRUCTURE:
         do 400 i = 1,n
            xyzfit(1,i) = r0(1,i) + xs
            xyzfit(2,i) = r0(2,i) + ys
            xyzfit(3,i) = r0(3,i) + zs
400      continue
*........CALCULATE RMS-DEVIATION "BY HAND" SINCE IT IS NOT AVAILABLE
*........FROM SR. "LSQQTN" AS SMALLEST EIGENVALUE OF THE MATRIX "M":
         error = 0.0d0
         do 410 i = 1,n
            error = error + (xyzfit(1,i)-xyz(1,i))**2
     .                    + (xyzfit(2,i)-xyz(2,i))**2
     .                    + (xyzfit(3,i)-xyz(3,i))**2
410      continue
         error = error/DFLOAT(n)
 
      endif
 
*.....A PURE ROTATIONAL FIT HAS TO BE PERFORMED:
 
      if (puretr .eq. 2) then
*........DO THE ROTATIONAL FIT:
         iopt = 1
         call LSQQTN(n,r,r0,wr,q,rfit,error,iopt,iret2,msg2)
         if (iret2 .eq. 99) goto 9020
         iret = iret2
         write(msg,'(a11,a69)')' RIGFIT -> ',msg2
*........SET THE POSITION DIFFERENCE OF THE CMs/CRs TO ZERO:
         dcm(1) = 0.0d0
         dcm(2) = 0.0d0
         dcm(3) = 0.0d0
*....... NO SHIFT FOR THE FITTED CM/CR-COORDINATES:
         do 500 i = 1,n
            xyzfit(1,i) = rfit(1,i)
            xyzfit(2,i) = rfit(2,i)
            xyzfit(3,i) = rfit(3,i)
500      continue
 
      endif
 
*---- STOP: ---------------------------------------------------------*
 
      return
 
8900  iret = 99
      msg = ' sr RIGFIT: PURETR must be 0 or 1 or 2'
      return
 
9000  iret = 99
      write(msg,'(a25,i5)')' RIGFIT: NMAX must be >= ',n
      return
 
9010  iret = 99
      write(msg,'(a36,i5)')' RIGFIT: Nonsense input for N, N = ',N
      return
 
9020  iret = 99
      write(msg,'(a11,a69)')' RIGFIT -> ',msg2
      return
 
      end
