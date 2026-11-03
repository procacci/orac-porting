      SUBROUTINE write_gofrp(rewind,fstep,igr,maxint,offset,wca,whe,ntap
     &     ,delrg,gofr_cut,ngrdon,iret,errmsg)

************************************************************************
*   Time-stamp: <97/08/21 21:31:20 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Dec  7 1995 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      use unit

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER iret,maxint,offset,ngrdon,ntap,igr(maxint,*)
      CHARACTER*80  errmsg
      REAL*8  gofr_cut,delrg,wca(*),whe(*),fstep
      LOGICAL rewind

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,k,max,numa,numb,nca,nhe,kmmax,kmax
      PARAMETER (kmmax = 10000, kmax=1500)
      REAL*8  delrgi,facrdf,rcon,qr,gint,fns,gr,sumca,sumhe,xcut
     &     ,sk(kmmax),gor(kmmax),ssk,rk,rp,vol,aux,ssk0
      COMMON /rag1/ sk,gor

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      IF(kmmax .LT. maxint) THEN
         iret=1
         errmsg=
     &'In WRITE_GOFRP: Physical dimensions of G of R exceed limit of'
     &//' 4000. Abort.'
         RETURN
      END IF
      IF(ngrdon.LE.0) RETURN
      IF(rewind) REWIND kgofr_sk
      WRITE(kgofr_sk,'(/''Time Step ='',f15.5)') fstep
      xcut=gofr_cut**2
      delrgi=1.0D+0/delrg
      facrdf=1.0D0/(4.0d0*pi*delrg)
      max=MIN0(IDINT(delrgi*DSQRT(xcut)+0.5d0),(maxint-1))
      vol=(4.0D0*pi/3.0D0)*(DFLOAT(max)*delrg)**3
      sumca=0.0D0
      sumhe=0.0D0
      DO i=1,ntap
         sumca=sumca+wca(i)
         sumhe=sumhe+whe(i)
      END DO
      vol=0.0D0
      DO j=1,max
         rp=j*delrg
         vol=vol+rp*rp*delrg
      END DO
      vol=vol*4.0D0*pi

      nca=NINT(sumca)
      nhe=NINT(sumhe)
      numa=nca
      numb=nca
      IF(nca .NE. 0) THEN
         rcon=2.0D0/DFLOAT(ngrdon*numa)
         WRITE(kgofr_sk,1000) 
         qr=0.0d0
         gint=0.0d0
         DO i=1,max
            qr=qr+delrg
            fns=rcon*DFLOAT(igr(i,1+offset))
            gint=gint+fns
            gr=facrdf*fns/(qr*qr)
            gor(i)=gr
            WRITE(kgofr_sk,500) qr,gr,gint
         END DO
         WRITE(kgofr_sk,1500) 
         DO k=1,kmax
            rk=2.0D0*pi*DFLOAT(k)/kmax
            ssk=0.0D0
            aux=0.0D0
            DO j=1,max
               rp=j*delrg
               ssk=ssk+rp*rp*(gor(j)-(gint+1.0)/vol)*DSIN(rk*rp)*delrg
     &              /(rk*rp)
            END DO
            ssk=1.0D0+4.0D0*pi*ssk
            WRITE(kgofr_sk,700) rk,ssk
         END DO
      END IF

      numa=nhe
      numb=nhe
      IF( nhe .NE. 0) THEN
         rcon=2.0D0/DFLOAT(ngrdon*numa)
         WRITE(kgofr_sk,2000) 
         qr=0.0d0
         gint=0.0d0
         DO i=1,max
            qr=qr+delrg
            fns=rcon*DFLOAT(igr(i,2+offset))
            gint=gint+fns
            gr=facrdf*fns/(qr*qr)
            gor(i)=gr
            WRITE(kgofr_sk,500) qr,gr,gint
         END DO
         WRITE(kgofr_sk,2500) 
         DO k=1,kmax
            rk=2.0D0*pi*DFLOAT(k)/kmax
            ssk=0.0D0
            DO j=1,max
               rp=j*delrg
               ssk=ssk+rp*rp*(gor(j)-(gint+1.0)/vol)*DSIN(rk*rp)*delrg
     &              /(rk*rp)
            END DO
            ssk=1.0D0+4.0D0*pi*ssk
            WRITE(kgofr_sk,700) rk,ssk 
         END DO
      END IF

      numa=ntap
      numb=ntap
      rcon=2.0D0/DFLOAT(ngrdon*numa)
      WRITE(kgofr_sk,3000) 
      qr=0.0d0
      gint=0.0d0
      DO i=1,max
         qr=qr+delrg
         fns=rcon*DFLOAT(igr(i,3+offset))
         gint=gint+fns
         gr=facrdf*fns/(qr*qr)
         gor(i)=gr
         WRITE(kgofr_sk,500) qr,gr,gint
      END DO
      WRITE(kgofr_sk,3500) 
      DO k=1,kmax
         rk=2.0D0*pi*DFLOAT(k)/kmax
         ssk=0.0D0
         DO j=1,max
            rp=j*delrg
            ssk=ssk+rp*rp*(gor(j)-(gint+1.0)/vol)*DSIN(rk*rp)*delrg
     &           /(rk*rp)
         END DO
         ssk=1.0D0+4.0D0*pi*ssk
         WRITE(kgofr_sk,700) rk,ssk 
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

500   FORMAT(f15.6,e15.6,f15.6)
700   FORMAT(2f15.8)
1000  FORMAT(/'----------------  CA - CA  G ( r ) ---------------'/)
2000  FORMAT(/'----------------  HE - HE  G ( r ) ---------------'/)
3000  FORMAT(/'----------------  AL - AL  G ( r ) ---------------'/)
1500  FORMAT(/'----------------  CA - CA  S ( K ) ---------------'/)
2500  FORMAT(/'----------------  HE - HE  S ( K ) ---------------'/)
3500  FORMAT(/'----------------  AL - AL  S ( K ) ---------------'/)
      RETURN
      END
