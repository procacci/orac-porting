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


*======================= DECLARATIONS ==================================

      use parst
      use unit, only:kprint
      use omp_integr, only:omp_timing,time_cn,tcpu_cn,nthr,m1t
     &     ,array2_omp,nthr2,omp_dynamic
#ifdef _OMP_
#ifdef _BGQ_
      use omp_lib !omp layer 
#endif
#endif

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER ncnstr,nato,iret,nprot,mim_lim
#ifdef _OMP_
      INTEGER cnst(2,*),cnst_protl(m1t,nthr2)
#ifndef _BGQ_
      include 'omp_lib.h'       ! OMP layer
#endif
#else
      INTEGER cnst(2,*),cnst_protl(*)
#endif
      REAL*8  xp1(nato),yp1(nato),zp1(nato), xp0(nato),yp0(nato)
     &     ,zp0(nato),vpx(nato),vpy(nato),vpz(nato),dss(*),coeff(*)
     &     ,mass(*),dnit,dt
      CHARACTER*80 errmsg

*----------------------- VARIABLES IN COMMON --------------------------*

      INTEGER nmat
      PARAMETER (nmat = 20)
      REAL*8  x0(m9),y0(m9),z0(m9),mass0(m1)
      REAL*8  mat(nmat,nmat),gam(nmat),gamo(nmat),matx(2*nmat*nmat)
     &     ,maty(2*nmat*nmat),matz(2*nmat*nmat),xc(nmat),yc(nmat)
     &     ,zc(nmat),dd(nmat),aux,xd,yd,zd
      LOGICAL mask(m1)
      INTEGER ipiv(2*nmat),info

*-------------------- LOCAL VARIABLES ----------------------------------

      LOGICAL redo
      REAL*8  tela1,tela2,tcpu1,tcpu2
      REAL*8  treal
      integer t1,t2
      INTEGER la,lb,k,iter,iox,i,cnstp,iter1,count,ka,lab,lbb,k1,k2,kb
     &     ,count1,iredo
      REAL*8 tol,xab,yab,zab,dpx,dpy,dpz,dpp,dps,gg,amsla,amslb,dpax
     &     ,dpay,dpaz,dpbx,dpby,dpbz
      REAL*8  gg1,gg2,dcnt,dti,xk,yk,zk,tol_mim,two,aux1,aux2,aux3,det
     &     ,zero,four,a11,a12,a13,a21,a22,a23,a31,a32,a33,b11,b12,b13
     &     ,b21,b22,b23,b31,b32,b33
      DATA tol/1.0D-8/tol_mim/1.0D-7/two/2.0D0/zero/0.0D0/four/4.0D0/

*------------  OMP related variable -----------------------------------

      INTEGER itask   ! these are local ifndef _OMP_

*==================== EXECUTABLE STATEMENTS ============================


      if(omp_timing)  call timer(treal,tcpu1,tela1)
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

#ifdef _OMP_
      if(omp_dynamic) call omp_set_dynamic(.TRUE.)
!$OMP  PARALLEL DEFAULT(PRIVATE) NUM_THREADS(nthr2)
!$OMP& SHARED(nprot,cnst_protl,mim_lim,cnst,nthr,dss)
!$OMP& SHARED(x0,y0,z0,coeff,tol,tol_mim,mass0,mask,two,zero)
!$OMP& SHARED(xp1,yp1,zp1,vpx,vpy,vpz,four,dti,array2_omp,kprint)
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
     & ' While SHAKEing : The iteration procedure did not'
     x              //' converge.'
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

5000        redo=.false.
            iredo=0
            DO ka=1,cnstp
               gam(ka)=0.0D0
            END DO

*=========== Build constraint matrix ====================================

            DO ka=1,cnstp
#ifdef _OMP_
               k1=cnst_protl(1+count+ka,itask)
#else
               k1=cnst_protl(1+count+ka)
#endif
               la=cnst(1,k1)
               lb=cnst(2,k1)
               xc(ka)=xp1(la)-xp1(lb)
               yc(ka)=yp1(la)-yp1(lb)
               zc(ka)=zp1(la)-zp1(lb)
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
               dd(ka)=dss(k1)-(xc(ka)*xc(ka)+yc(ka)*yc(ka)+zc(ka)
     &              *zc(ka))
               
*============= Second loop on constraints k2 ===========================
                     
               DO kb=1,cnstp
#ifdef _OMP_
                  k2=cnst_protl(1+count+kb,itask)
#else
                  k2=cnst_protl(1+count+kb)
#endif
                  lab=cnst(1,k2)
                  lbb=cnst(2,k2)
                  
                  xd=x0(k2)
                  yd=y0(k2)
                  zd=z0(k2)
                  
                  aux1=mass0(la)
                  aux2=mass0(lb)
                  
                  IF(la .NE. lab .AND. la .NE. lbb .AND. lb .NE.
     &                 lab.AND. lb .NE. lbb) THEN
                     mat(ka,kb)=0.0D0
                  ELSE
                     dpax=0.0D0
                     dpay=0.0D0
                     dpaz=0.0D0
                     dpbx=0.0D0
                     dpby=0.0D0
                     dpbz=0.0D0
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

                     aux3=(dpx*xc(ka)+dpy*yc(ka)+dpz*zc(ka))*aux1

                     dpx=dpbx
                     dpy=dpby
                     dpz=dpbz
                     count1=count1+1
                     matx(count1)=dpbx
                     maty(count1)=dpby
                     matz(count1)=dpbz

                     aux3=aux3+(dpx*xc(ka)+dpy*yc(ka)+dpz*zc(ka))*aux2
                     aux3=aux3*four

                     mat(ka,kb)=-aux3


                  END IF
               END DO
            END DO

*=======================================================================
*---- Iteration to solve the non linear system -------------------------
*=======================================================================

            iox=1
            iter=0
            info=0
            DO WHILE(iox .EQ. 1)
               DO ka=1,cnstp
                  gamo(ka)=gam(ka)
               END DO

               IF(iter .EQ. 0) THEN
                  DO ka=1,cnstp
                     gam(ka)=dd(ka)
                  END DO
                  IF(cnstp .EQ. 1) THEN
                     gam(1)=gam(1)/mat(1,1)
                  ELSEIF(cnstp .EQ. 2) THEN
                     b11=mat(2,2)
                     b12=-mat(2,1)
                     b21=-mat(1,2)
                     b22=mat(1,1)
                     det=b11*b22-b21*b12
                     det=1.0D0/det
                     xd=b11*gam(1)+b21*gam(2)
                     yd=b12*gam(1)+b22*gam(2)
                     gam(1)=xd*det
                     gam(2)=yd*det
                  ELSE IF(cnstp .EQ. 3) THEN
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
                  END IF
               ELSE
                  count1=0
                  DO ka=1,cnstp
                     xd=0.0D0
                     yd=0.0D0
                     zd=0.0D0
                     gg1=mass0(la)*two
                     gg2=mass0(lb)*two
                     DO kb=1,cnstp
                        count1=count1+1
                        aux=gam(kb)*gg1
                        xd=xd+matx(count1)*aux
                        yd=yd+maty(count1)*aux
                        zd=zd+matz(count1)*aux
                        count1=count1+1
                        aux1=-matx(count1)
                        aux2=-maty(count1)
                        aux3=-matz(count1)
                        aux=gam(kb)*gg2
                        xd=xd+aux1*aux
                        yd=yd+aux2*aux
                        zd=zd+aux3*aux
                     END DO
                     xc(ka)=dd(ka)-(xd*xd+yd*yd+zd*zd)
                  END DO

                  DO ka=1,cnstp
                     gam(ka)=xc(ka)
                  END DO

                  IF(cnstp .EQ. 1) THEN

                     gam(1)=gam(1)/mat(1,1)
                  ELSE IF(cnstp .EQ. 2) THEN
                     xd=b11*gam(1)+b21*gam(2)
                     yd=b12*gam(1)+b22*gam(2)
                     gam(1)=xd*det
                     gam(2)=yd*det
                  ELSE IF(cnstp .EQ. 3) THEN
                     xd=b11*gam(1)+b21*gam(2)+b31*gam(3)
                     yd=b12*gam(1)+b22*gam(2)+b32*gam(3)
                     zd=b13*gam(1)+b23*gam(2)+b33*gam(3)
                     gam(1)=xd*det
                     gam(2)=yd*det
                     gam(3)=zd*det
                  ELSE
                     CALL dgetrs('N',cnstp,1,mat,nmat,ipiv,gam,nmat,info
     &                    )
                  END IF
               END IF
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
               iter=iter+1
               IF(iter.GT.120)THEN
                 if(iredo.le.4) THEN
                   iret=1
                   errmsg=
     & ' While constraining with MIM: The iteration procedure did not'
     x              //' converge.'
                   CALL xerror(errmsg,80,1,20)
                   DO ka=1,cnstp
#ifdef _OMP_
                     k1=cnst_protl(1+count+ka,itask)
#else
                     
                     k1=cnst_protl(1+count+ka)
#endif
                     la=cnst(1,k1)
                     lb=cnst(2,k1)
#ifdef _OMP_
                     write(kprint,132) la,lb,itask
132                  format(" -- constraint betweem atoms:",2i7,
     &                    "   thread =",i3)
#else
                     write(kprint,132) la,lb
132                  format(" -- constraint betweem atoms:",2i7)
#endif
                   end do
                   redo=.true.
                   iredo=iredo+1
                 ELSE
                   errmsg="MIM failed four times. Exiting.."  
                   write(kprint,*) " constraint ->",la,lb
                   do ka=1,nato
                     write(kprint,223) ka,1/mass0(ka),xp0(ka),yp0(ka)
     &                    ,zp0(ka),xp1(ka),yp1(ka)
     &                    ,zp1(ka)
223                  format(i7,f8.2,3f15.6,5x,3f15.6)
                   end do
                   CALL xerror(errmsg,80,1,2)
                 END IF
               END IF
               iox=0
               DO ka=1,cnstp
                  IF(gam(ka) .NE. 0) THEN
                     aux=(gam(ka)-gamo(ka))/gam(ka)
                     aux=DABS(aux)
                     IF(aux .GT. tol_mim) iox=1
                  END IF
               END DO
               if(redo) iox=0
            END DO

            iter1=iter1+iter

*=======================================================================
*---- Compute corrected position and velocities ------------------------
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
                     aux=gam(kb)*xd
                     aux1=matx(count1)
                     aux2=maty(count1)
                     aux3=matz(count1)
                     xp1(la)=xp1(la)+aux1*aux
                     yp1(la)=yp1(la)+aux2*aux
                     zp1(la)=zp1(la)+aux3*aux
                     aux=aux*dti
                     vpx(la)=vpx(la)+aux1*aux
                     vpy(la)=vpy(la)+aux2*aux
                     vpz(la)=vpz(la)+aux3*aux
                  END IF
                  count1=count1+1
                  IF(mask(lb)) THEN
                     aux=gam(kb)*yd
                     aux1=matx(count1)
                     aux2=maty(count1)
                     aux3=matz(count1)
                     xp1(lb)=xp1(lb)+aux1*aux
                     yp1(lb)=yp1(lb)+aux2*aux
                     zp1(lb)=zp1(lb)+aux3*aux
                     aux=aux*dti
                     vpx(lb)=vpx(lb)+aux1*aux
                     vpy(lb)=vpy(lb)+aux2*aux
                     vpz(lb)=vpz(lb)+aux3*aux
                  END IF
               END DO
               mask(la)=.FALSE.
               mask(lb)=.FALSE.
            END DO
            if(.not.redo) THEN 
              count=count+cnstp+1
            ELSE
              DO ka=1,cnstp
#ifdef _OMP_
                k1=cnst_protl(1+count+ka,itask)
#else
                k1=cnst_protl(1+count+ka)
#endif
                la=cnst(1,k1)
                lb=cnst(2,k1)
                mask(la)=.true.
                mask(lb)=.true.
              end do
              goto 5000
            end if
         END IF
#ifdef _OMP_ 
         endif  ! close OpenMP task condition
#endif
      END DO
#ifdef _OMP_
!$OMP END PARALLEL
      if(omp_dynamic) THEN
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
*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
