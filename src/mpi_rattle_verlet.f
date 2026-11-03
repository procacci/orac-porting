      SUBROUTINE rattle_verlet(dt,xp1,yp1,zp1,xp0,yp0,zp0,vpx,vpy,vpz
     &     ,nato,cnst,dss,coeff,ncnstr,mass,dnit,nprot,cnst_protl
     &     ,mim_lim,iret,errmsg)

************************************************************************
*   Time-stamp: <98/03/03 12:46:09 marchi>                             *
*                                                                      *
*     It will supply a new set of coordinates for which                *
*     the set of constraints are satisfied.                            *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Feb 13 1997 -                                     *
*                                                                      *
************************************************************************

!  Parallel version only use SHAKE algorithm nprot is the number of
!  "constraints networks" (CN). If in a molecule all bonds are
!  constrained then the CN coincides with the actual molecule.  If in a
!  molecule one bond (not belonging to a cycle) is not constrained, then
!  the molecule is splitted in two CNs. Constraints need to be satisfied
!  in these two sub units independently due to the intercalated
!  stretching separating the molecule in the two distinct CNs
!  Parallelization may hence be done on the nprot independent CNs using
!  the stride mechanism as usual.  If the CN distribution is
!  pathologically disparate, e.g. a single solvated protein making up
!  one single CN out of nprot and nprot-1 small solvent CN-molecules,
!  then the program insert bonds in the large CN so as to divide it in a
!  set of M smaller identical CNs such that M is exactly equal to
!  NPROCf. In this each thread is assigned exactly one a "big" CN
!  sub-molecole. This optimizes the balance among threads. The inserted
!  bond is always between heavy atoms with the corresponding FF force
!  constant.  TODO: to achieve the above it suffices to erase the
!  intercalated constraints in the array cnstp(:,:) and to add the
!  corresponding stretching at the end of the bond list (as done in
!  steer dynamics) augmenting potbo lstretch etc etc. The get_prot_cnstr
!  (SHAKE version) that finds out the number of independent constraint
!  networks can be called safely after the cnstp update.
*======================= DECLARATIONS ==================================

      use parst

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER ncnstr,nato,iret,nprot,mim_lim
      INTEGER cnst(2,*),cnst_protl(*)
      REAL*8  xp1(*),yp1(*),zp1(*), xp0(*),yp0(*),zp0(*),vpx(*)
     &     ,vpy(*),vpz(*),dss(*),coeff(*),mass(*),dnit,dt
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN COMMON --------------------------*

      INTEGER nmat
      PARAMETER (nmat = 4)
      REAL*8  x0(m9),y0(m9),z0(m9),mass0(m1)
      REAL*8  mat(nmat,nmat),gam(nmat),gamo(nmat),matx(2*nmat*nmat)
     &     ,maty(2*nmat*nmat),matz(2*nmat*nmat),xc(nmat),yc(nmat)
     &     ,zc(nmat),dd(nmat),aux,xd,yd,zd
      LOGICAL mask(m1)
      INTEGER ipiv(2*nmat),info
      COMMON /rag1/ x0,y0,z0,mass0,mat,gam,gamo,matx,maty,matz,ipiv,mask

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER la,lb,k,iter,iox,i,cnstp,iter1,count,ka,lab,lbb,k1,k2,kb
     &     ,count1
      REAL*8 tol,xab,yab,zab,dpx,dpy,dpz,dpp,dps,gg,amsla,amslb,dpax
     &     ,dpay,dpaz,dpbx,dpby,dpbz
      REAL*8  gg1,gg2,dcnt,dti,xk,yk,zk,tol_mim,two,aux1,aux2,aux3,det
     &     ,zero,four,a11,a12,a13,a21,a22,a23,a31,a32,a33,b11,b12,b13
     &     ,b21,b22,b23,b31,b32,b33
      DATA tol/1.0D-8/tol_mim/1.0D-7/two/2.0D0/zero/0.0D0/four/4.0D0/

*==================== EXECUTABLE STATEMENTS ============================

      IF(ncnstr .EQ. 0) RETURN
      IF(nmat .LT. mim_lim) THEN
         errmsg=' While SHAKEing dimension of the cnstr matrix exceeds'
     &        //' physical dimanesions.'
         iret=1
         RETURN
      END IF

      IF(m9 .LT. ncnstr) THEN
          errmsg=' In CNSTRP. Physical dimensions of working arrays'//
     x ' are insufficient. Abort!'
          iret=1
          WRITE(6,'('' M9 = '',i6,'' MCNSTR = '',i6)') m9,ncnstr
          RETURN
      END IF

*=======================================================================
*---- Copy the coordinates of the molecule IA to a temporary array -----
*=======================================================================


      if(dt.ne.0.D0) THEN 
        dti=1.0D0/dt
      ELSE
        dti=0.d0
      END IF

      DO i=1,nato
         mass0(i)=1.0D0/mass(i)
         mask(i)=.TRUE.
      END DO

*=======================================================================
*---- Compute the vectors associated with each constraint and store ----
*---- them in a temporary array ----------------------------------------
*=======================================================================

      DO k=1,ncnstr
         la=cnst(1,k)
         lb=cnst(2,k)
         xab=xp0(la)-xp0(lb)
         yab=yp0(la)-yp0(lb)
         zab=zp0(la)-zp0(lb)
         x0(k)=xab
         y0(k)=yab
         z0(k)=zab
      END DO

*=======================================================================
*---- SHAKE loop -------------------------------------------------------
*=======================================================================

      iter1=0
      count=0
      DO i=1,nprot
        cnstp=cnst_protl(1+count)
        iter=0
1000    CONTINUE
        iox=0
        DO ka=1,cnstp
          k=cnst_protl(1+count+ka)
          la=cnst(1,k)
          lb=cnst(2,k)
          dpx=xp1(la)-xp1(lb)
          dpy=yp1(la)-yp1(lb)
          dpz=zp1(la)-zp1(lb)
          dpp=dpx**2+dpy**2+dpz**2
          dcnt=dpp-dss(k)
          IF(DABS(dcnt).GE.tol) THEN
            iox=1
            xk=x0(k)
            yk=y0(k)
            zk=z0(k)
            amsla=mass0(la)
            amslb=mass0(lb)
            dps=dpx*xk+dpy*yk+dpz*zk
            gg=dcnt/(dps*coeff(k))
            gg1=-gg*amsla
            gg2=gg*amslb
            xp1(la)=xp1(la)+xk*gg1
            xp1(lb)=xp1(lb)+xk*gg2
            
            yp1(la)=yp1(la)+yk*gg1
            yp1(lb)=yp1(lb)+yk*gg2
            
            zp1(la)=zp1(la)+zk*gg1
            zp1(lb)=zp1(lb)+zk*gg2
            
            gg1=gg1*dti
            gg2=gg2*dti
            
            vpx(la)=vpx(la)+xk*gg1
            vpx(lb)=vpx(lb)+xk*gg2
            
            vpy(la)=vpy(la)+yk*gg1
            vpy(lb)=vpy(lb)+yk*gg2
            
            vpz(la)=vpz(la)+zk*gg1
            vpz(lb)=vpz(lb)+zk*gg2
          END IF
        END DO
        IF(iox.ne.0) THEN
          iter1=iter1+1
          iter=iter+1
          IF(iter.GT.5000)THEN
            iret=1
            errmsg=
     &           ' While SHAKEing : The iteration procedure did not'
     x           //' converge.'
            RETURN
          END IF
          GOTO 1000
        END IF
        count=count+cnstp+1
      END DO

      dnit=dnit+DFLOAT(iter1)/DFLOAT(nprot)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
