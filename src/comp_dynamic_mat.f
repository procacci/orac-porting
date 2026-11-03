      SUBROUTINE comp_dynamic_mat(mapnl,mapdn,nmapdn,nnlpp0,nnlpp
     &     ,fudgec,xp0,yp0,zp0,ma,mb,fpx1,fpy1,fpz1,d_mat,wk
     &     ,eigvl,eigvc)

************************************************************************
*   Time-stamp: <98/03/07 15:14:40 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Feb 26 1998 -                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*    comp_dynamic_mat externals:       	       	       	       	       *
*          eigrs f1dim_der interpol_dyna low_up xerror zero0	       *
*       							       *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit
      use parst
      use cpropar

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER ma,mb,mapnl(*),mapdn(2,*),nmapdn(*)
      INTEGER nnlpp0(*),nnlpp(*)
      REAL*8  fudgec
      REAL*8  xp0(*),yp0(*),zp0(*),fpx1(ma,*),fpy1(ma,*),fpz1(ma,*)
     &     ,d_mat(*),wk(*),eigvl(*),eigvc(mb,*)

*------------------------- LOCAL VARIABLES ----------------------------*

      REAL*8   hst
      INTEGER nstep,i,j,ia,ib,k,l,kf,ier
      CHARACTER*80 errmsg
      REAL*8  fc(cheb_order),ccoeff(cheb_order),cder(cheb_order),dyna
     &     ,f1dim_der
      REAL*8  hstep,aux1,aux2,utotal,a,b,bma,bpa,y,w1
     &     ,w2,w,unit_freq
      CHARACTER*1 type
      EXTERNAL f1dim_der

*----------------------- EXECUTABLE STATEMENTS ------------------------*
      
      IF(n_mat .LT. 3*ntap) THEN
         errmsg=' Dimensions the dynamical matrix too small. Increase '
     &        //' _DYNAMIC_DIM_ in config.H'
         CALL xerror(errmsg,80,1,2)
         STOP
      END IF
      IF(cheb_order .LT. nstep_freq) THEN
         errmsg=' Input Chebyshev order too large. Increase '
     &        //' _CHEB_ORDER_ in config.H'
         CALL xerror(errmsg,80,1,2)
         STOP
      END IF

*=======================================================================
*---- Pick an initial step size ----------------------------------------
*=======================================================================

      hstep=hstep_freq
      nstep=nstep_freq
      b=-hstep
      a=hstep
      bma=0.5D0*(b-a)
      bpa=0.5D0*(b+a)

*=======================================================================
*---- Zero dynamical matrix --------------------------------------------
*=======================================================================

      CALL zero0(eigvc,n_mat*n_mat)

*=======================================================================
*---- Loop to compute dynamical matrix numerically ---------------------
*=======================================================================


      DO i=1,ntap
         DO j=1,3

*---- Compute Forces at the Chebyshev points and build up a polynomial

            DO k=1,nstep
               y=DCOS(PI*(k-0.5D0)/nstep)
               hst=y*bma+bpa
               utotal=f1dim_der(i,j,mapnl,mapdn,nmapdn,nnlpp0,nnlpp
     &              ,fudgec,hst,xp0,yp0,zp0,fpx1(1,k)
     &              ,fpy1(1,k),fpz1(1,k))
            END DO
            ia=(i-1)*3+j
            DO k=1,ntap
               ib=(k-1)*3

*---- Interpolate the derivative of the polynomial for each particle 
*---- coordinate

               DO l=1,nstep
                  fc(l)=-fpx1(k,l)
               END DO
               CALL interpol_dyna(a,b,ccoeff,cder,fc,nstep,dyna)
               eigvc(ia,ib+1)=eigvc(ia,ib+1)+dyna

               DO l=1,nstep
                  fc(l)=-fpy1(k,l)
               END DO
               CALL interpol_dyna(a,b,ccoeff,cder,fc,nstep,dyna)
               eigvc(ia,ib+2)=eigvc(ia,ib+2)+dyna

               DO l=1,nstep
                  fc(l)=-fpz1(k,l)
               END DO
               CALL interpol_dyna(a,b,ccoeff,cder,fc,nstep,dyna)
               eigvc(ia,ib+3)=eigvc(ia,ib+3)+dyna
            END DO
         END DO
      END DO

*=======================================================================
*---- Symmetrize the dynamical matrix ----------------------------------
*=======================================================================

      DO i=1,ntap*3-1
         DO j=i+1,ntap*3
            aux1=eigvc(i,j)
            aux2=eigvc(j,i)
            aux1=(aux1+aux2)*0.5D0
            eigvc(i,j)=aux1
            eigvc(j,i)=aux1
         END DO
      END DO

c---  trasform cartesian to mass weighted cartesian 
 
      DO i=1,ntap 
         DO k=1,3
            ia = 3*(i-1)+k
            DO j=1,ntap
               DO l=1,3
                  ib=3*(j-1)+l
                  eigvc(ia,ib)=eigvc(ia,ib)/DSQRT(mass(i)*mass(j))
               END DO
            END DO
         END DO
      END DO

      k=0
      DO i=1,ntap*3
         DO j=1,i
            k=k+1
            d_mat(k)=eigvc(i,j)
         END DO
      END DO

c--------------------------------------------------------------------------
c    Diagonalization of dynamical matrix and print out
c--------------------------------------------------------------------------

      WRITE(kprint,100)
      CALL eigrs(d_mat,3*ntap,2,eigvl,eigvc,n_mat,wk,ier)

      kf = 0
      unit_freq=1.0D0/2.997925D10/unitt/(2.0D0*pi)
      DO i=1,3*ntap
         w2 = eigvl(i) 
         w1 = DSQRT(DABS(w2))
         w = 0.0D0
c--     write for xmol format 
         IF(w1 .gt. 0.0D0) THEN 
            kf = kf + 1
            w = w2/w1*unit_freq
            WRITE(kfreq,'(i7)') ntap 
            WRITE(kfreq,'(i5,''-th Frequency (cm-1):'',f12.4)') kf,w
            DO j=1,ntap
               type(1:1)=beta(j)(1:1)
               CALL low_up(type,1)
               WRITE(kfreq,'(a1,3f12.4,3f10.5)') type,xp0(j),yp0(j)
     &              ,zp0(j),(eigvc((j-1)*3+k,i),k=1,3)
            END DO 
         END IF
      END DO

      WRITE(kprint,200) 3*ntap
      WRITE(kprint,300) (i,eigvl(i)/DSQRT(DABS(eigvl(i)))*unit_freq,i=1
     &     ,3*ntap)
      WRITE(kprint,400) 

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

100   FORMAT(80('.'))
200   FORMAT(//
     &     21x,'***************************************'/
     &     21x,'*     ',i6,'  Frequencies in cm-1     *'/
     &     21x,'***************************************'//)
300   FORMAT(4(1x,i4,1x,f10.3,2x))
400   FORMAT(//
     &     21x,'***************************************'/
     &     21x,'*         Eigenvectors Written        *'/
     &     21x,'*            to Output File           *'/
     &     21x,'***************************************'//)
      RETURN
      END
