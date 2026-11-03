      SUBROUTINE rattle_correc(dt,xp0,yp0,zp0,vpx,vpy,vpz,nato
     &     ,cnst,dss,coeff,ncnstr,mass,dnit,nprot,cnst_protl,mim_lim
     &     ,iret,errmsg)

************************************************************************
*   Time-stamp: <98/03/03 12:43:51 marchi>                             *
*                                                                      *
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

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use parst
      use omp_integr, only:omp_timing,time_cn,tcpu_cn,nthr,m1t,nthr2
     &     ,array2_omp,omp_dynamic
#ifdef _OMP_
#ifdef _BGQ_
      use omp_lib !omp layer 
#endif
#endif

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      REAL*8  dt,xp0(nato),yp0(nato),zp0(nato),vpx(nato)
     &     ,vpy(nato),vpz(nato),dss(*),coeff(*),mass(*),dnit
#ifdef _OMP_
      INTEGER nato,ncnstr,iret,cnst(2,*),nprot,cnst_protl(m1t,nthr2)
     &     ,mim_lim
#ifndef _BGQ_
      include 'omp_lib.h'       ! OMP layer
#endif
#else
      INTEGER nato,ncnstr,iret,cnst(2,*),nprot,cnst_protl(*),mim_lim
#endif
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN COMMON --------------------------*

      REAL*8  xp1(m9),yp1(m9),zp1(m9),vpx1(m1),vpy1(m1),vpz1(m1)
     &     ,mass0(m1),dssi(m9),coeffi(m9)
      LOGICAL mask(m1)

      INTEGER nmat
      PARAMETER (nmat = 20)
      REAL*8  mat(nmat,nmat),gam(nmat),matx(2*nmat*nmat)
     &     ,maty(2*nmat*nmat),matz(2*nmat*nmat),xc(nmat),yc(nmat)
     &     ,zc(nmat),dd(nmat),aux,xd,yd,zd,aux1,aux2,aux3,xd1,yd1,zd1
      INTEGER ipiv(2*nmat),info

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8  tela1,tela2,tcpu1,tcpu2
      REAL*8  treal
      integer t1,t2
      integer clock_rate, clock_max 
      INTEGER i,j,k,iter,iox,la,lb,cnstp,iter1,count,ka,kb,count1,k1,k2
     &     ,lab,lbb
      REAL*8  xpa,ypa,zpa,gab,dpx,dpy,dpz,vdpx,vdpy,vdpz,dpp,gg1,gg2,xab
     &     ,yab,zab,tol,tolh,two,zero,dpax,dpay,dpaz,dpbx,dpby,dpbz,det
     &     ,a11,a12,a13,a21,a22,a23,a31,a32,a33,b11,b12,b13,b21,b22,b23
     &     ,b31,b32,b33
      DATA  tol/1.0D-6/two/2.0D0/zero/0.0D0/

*----------------------- EXECUTABLE STATEMENTS ------------------------*


*------------  OMP related variable -----------------------------------

      INTEGER itask ! these are local ifndef _OMP_

      if(omp_timing)  call timer(treal,tcpu1,tela1)
      IF(nmat .LT. mim_lim) THEN
         errmsg=' While RATTLEing dimension of the cnstr matrix exceeds'
     &        //' physical dimanesions.'
         iret=1
         RETURN
      END IF
      tolh=tol
      DO i=1,nato
         mass0(i)=1.0D0/mass(i)
         mask(i)=.TRUE.
      END DO
      DO k=1,ncnstr
         la=cnst(1,k)
         lb=cnst(2,k)
         xpa=xp0(la)-xp0(lb)
         ypa=yp0(la)-yp0(lb)
         zpa=zp0(la)-zp0(lb)
         xp1(k)=xpa
         yp1(k)=ypa
         zp1(k)=zpa
         coeffi(k)=1.0D0/coeff(k)
         dssi(k)=1.0D0/dss(k)
      END DO
#ifdef _OMP_
      if(omp_dynamic) call omp_set_dynamic(.TRUE.)
!$OMP  PARALLEL DEFAULT(PRIVATE) NUM_THREADS(nthr2)
!$OMP& SHARED(nprot,cnst_protl,mim_lim,cnst,mask)
!$OMP& SHARED(xp1,yp1,zp1,coeffi,tolh,dssi,mass0,two,zero)
!$OMP& SHARED(vpx,vpy,vpz,array2_omp)
      itask=1+OMP_GET_THREAD_NUM()
#endif
      iter1=0
      count=0
      DO i=1,nprot
#ifdef _OMP_
        if(array2_omp(i,itask)) THEN
         cnstp=cnst_protl(1+count,itask)
#else
         cnstp=cnst_protl(1+count)
#endif
         IF(cnstp .GT. mim_lim) THEN
            iter=0
1000        CONTINUE
            iox=0
            DO ka=1,cnstp
#ifdef _OMP_
               k=cnst_protl(1+count+ka,itask)
#else
               k=cnst_protl(1+count+ka)
#endif
               la=cnst(1,k)
               lb=cnst(2,k)
               dpx=xp1(k)
               dpy=yp1(k)
               dpz=zp1(k)
               vdpx=vpx(la)-vpx(lb)
               vdpy=vpy(la)-vpy(lb)
               vdpz=vpz(la)-vpz(lb)
               dpp=dpx*vdpx+dpy*vdpy+dpz*vdpz
               IF(DABS(dpp) .GT. tolh) THEN
                  gab=-dpp*dssi(k)*coeffi(k)
                  iox=1
                  gg1=mass0(la)*gab
                  gg2=-mass0(lb)*gab
                  vpx(la)=vpx(la)+gg1*dpx
                  vpx(lb)=vpx(lb)+gg2*dpx
                  
                  vpy(la)=vpy(la)+gg1*dpy
                  vpy(lb)=vpy(lb)+gg2*dpy
                  
                  vpz(la)=vpz(la)+gg1*dpz
                  vpz(lb)=vpz(lb)+gg2*dpz
               END IF
            END DO
            IF(iox.ne.0) THEN
               iter1=iter1+1
               iter=iter+1
               IF(iter.GT.5000)THEN
                 iret=1
                 errmsg='While RATTLEing: The iteration' 
     &               // ' procedure did not converge.'
                 CALL xerror(errmsg,80,1,2)
#ifdef _OMP_
!$OMP             BARRIER 
#endif
                 STOP
               END IF
               GOTO 1000
            END IF
         ELSE

*=======================================================================
*---  Use Matrix Inversion Method for Constraints ---------------------- 
*=======================================================================

*=========== Build constraint matrix ====================================

            DO ka=1,cnstp
#ifdef _OMP_
               k1=cnst_protl(1+count+ka,itask)
#else
               k1=cnst_protl(1+count+ka)
#endif
               la=cnst(1,k1)
               lb=cnst(2,k1)
               xc(ka)=vpx(la)-vpx(lb)
               yc(ka)=vpy(la)-vpy(lb)
               zc(ka)=vpz(la)-vpz(lb)
            END DO
               
*======== First loop on constraints k1 =================================


            count1=0
            DO ka=1,cnstp
#ifdef _OMP_
               k1=cnst_protl(1+count+ka,itask)
#else
               k1=cnst_protl(1+count+ka)
#endif
               la=cnst(1,k1)
               lb=cnst(2,k1)
               dd(ka)=xc(ka)*xp1(k1)+yc(ka)*yp1(k1)+zc(ka)*zp1(k1)
               xd1=xp1(k1)
               yd1=yp1(k1)
               zd1=zp1(k1)
               
*============= Second loop on constraints k2 ===========================
                     
               DO kb=1,cnstp
#ifdef _OMP_
                  k2=cnst_protl(1+count+kb,itask)
#else
                  k2=cnst_protl(1+count+kb)
#endif
                  lab=cnst(1,k2)
                  lbb=cnst(2,k2)
                  
                  xd=xp1(k2)
                  yd=yp1(k2)
                  zd=zp1(k2)
                  
                  aux1=mass0(la)
                  aux2=mass0(lb)
                  
                  IF(la .NE. lab .AND. la .NE. lbb .AND. lb .NE.
     &                 lab.AND. lb .NE. lbb) THEN
                     mat(ka,kb)=zero
                  ELSE
                     dpax=zero
                     dpay=zero
                     dpaz=zero
                     dpbx=zero
                     dpby=zero
                     dpbz=zero
                     IF(la .EQ. lab) THEN
                        dpax=-xd
                        dpay=-yd
                        dpaz=-zd
                     ELSE IF(la .EQ. lbb) THEN
                        dpax=xd
                        dpay=yd
                        dpaz=zd
                     END IF
                     
                     IF(lb .EQ. lab) THEN
                        dpbx=-xd
                        dpby=-yd
                        dpbz=-zd
                     ELSE IF(lb .EQ. lbb) THEN
                        dpbx=xd
                        dpby=yd
                        dpbz=zd
                     END IF
                     dpx=-dpax
                     dpy=-dpay
                     dpz=-dpaz

                     count1=count1+1
                     matx(count1)=dpax
                     maty(count1)=dpay
                     matz(count1)=dpaz
                     aux3=(xd1*dpx+yd1*dpy+zd1*dpz)*aux1

                     dpx=dpbx
                     dpy=dpby
                     dpz=dpbz
                     count1=count1+1
                     matx(count1)=dpbx
                     maty(count1)=dpby
                     matz(count1)=dpbz

                     aux3=aux3+(xd1*dpx+yd1*dpy+zd1*dpz)*aux2
                     aux3=aux3*two
                     mat(ka,kb)=-aux3
                  END IF
               END DO
            END DO

*=======================================================================
*---- Solve the linear system ------------------------------------------
*=======================================================================

            info=0
            DO ka=1,cnstp
               gam(ka)=-dd(ka)
            END DO
            IF(cnstp .EQ. 1) THEN
               gam(1)=gam(1)/mat(1,1)
            ELSE IF(cnstp .EQ. 2) THEN
               a11=mat(1,1)
               a12=mat(1,2)
               a21=mat(2,1)
               a22=mat(2,2)
               b11=a22
               b12=-a21
               b21=-a12
               b22=a11
               det=a11*a22-a21*a12
               det=1.0D0/det
               xd=b11*gam(1)+b21*gam(2)
               yd=b12*gam(1)+b22*gam(2)
               gam(1)=xd*det
               gam(2)=yd*det

            ELSEIF(cnstp .EQ. 3) THEN
               a11=mat(1,1)
               a12=mat(1,2)
               a13=mat(1,3)
               a21=mat(2,1)
               a22=mat(2,2)
               a23=mat(2,3)
               a31=mat(3,1)
               a32=mat(3,2)
               a33=mat(3,3)
               b11= a22*a33-a32*a23
               b12=-a21*a33+a31*a23
               b13= a21*a32-a22*a31
               b21=-a12*a33+a32*a13
               b22= a11*a33-a31*a13
               b23=-a11*a32+a31*a12
               b31= a12*a23-a22*a13
               b32=-a11*a23+a21*a13
               b33= a11*a22-a21*a12
               
               det=a11*b11+a21*b21+a31*b31
               det=1.0D0/det
               xd=b11*gam(1)+b21*gam(2)+b31*gam(3)
               yd=b12*gam(1)+b22*gam(2)+b32*gam(3)
               zd=b13*gam(1)+b23*gam(2)+b33*gam(3)
               gam(1)=xd*det
               gam(2)=yd*det
               gam(3)=zd*det
            ELSE
               CALL dgesv(cnstp,1,mat,nmat,ipiv,gam,nmat,info)
               IF(info .NE. 0) THEN
                  iret=1
                  errmsg=
     & ' While constraining with MIM: matrix inversion failed. '
                  CALL xerror(errmsg,80,1,2)
#ifdef _OMP_
!$OMP             BARRIER 
#endif
                  STOP
               END IF
            END IF

*=======================================================================
*---- Compute corrected positions and velocities -----------------------
*=======================================================================

            count1=0
            DO ka=1,cnstp
#ifdef _OMP_
               k1=cnst_protl(1+count+ka,itask)
#else
               k1=cnst_protl(1+count+ka)
#endif
               la=cnst(1,k1)
               lb=cnst(2,k1)
               xd=mass0(la)*two
               yd=mass0(lb)*two
               DO kb=1,cnstp
                  count1=count1+1
                  IF(mask(la)) THEN
                     aux=gam(kb)
                     aux1=matx(count1)*aux
                     aux2=maty(count1)*aux
                     aux3=matz(count1)*aux
                     vpx(la)=vpx(la)+aux1*xd
                     vpy(la)=vpy(la)+aux2*xd
                     vpz(la)=vpz(la)+aux3*xd
                  END IF
                  count1=count1+1
                  IF(mask(lb)) THEN
                     aux=gam(kb)
                     aux1=matx(count1)*aux
                     aux2=maty(count1)*aux
                     aux3=matz(count1)*aux
                     vpx(lb)=vpx(lb)+aux1*yd
                     vpy(lb)=vpy(lb)+aux2*yd
                     vpz(lb)=vpz(lb)+aux3*yd
                  END IF
               END DO
               mask(la)=.FALSE.
               mask(lb)=.FALSE.
            END DO
         END IF
         count=count+cnstp+1
#ifdef _OMP_
         ENDIF  ! close the OpenMP task condition
#endif
      END DO
#ifdef _OMP_
!$OMP END PARALLEL
      if(omp_dynamic)  THEN 
        call omp_set_dynamic(.FALSE.)
        call OMP_SET_NUM_THREADS(nthr) !restore thread numebr
      END IF
#endif 
      dnit=dnit+DFLOAT(iter1)/DFLOAT(nprot)
      if(omp_timing) THEN 
        call timer(treal,tcpu2,tela2)
        time_cn=time_cn+tela2-tela1
        tcpu_cn=tcpu_cn+tcpu2-tcpu1
      END IF
*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
