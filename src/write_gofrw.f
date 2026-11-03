      SUBROUTINE write_gofrw(rewind,fstep,igr,maxint,g1,delrg,gofr_cut
     &     ,itypes,beta_slv,volume,nspec,nmol,ngrdon,offset,iret,errmsg)

************************************************************************
*   Time-stamp: <97/08/21 11:26:58 marchi>                             *
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

      INTEGER nmol,maxint,nspec,g1,ngrdon,igr(maxint,*),offset,itypes(*)
      CHARACTER*1 beta_slv(*)
      REAL*8  vol,volume,delrg,gofr_cut,fstep
      LOGICAL rewind

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,j,isp,jsp,itype,max,lpoint,numa,numb,kmmax,iret,kmax,k
     &     ,nts
      CHARACTER*80  errmsg
      PARAMETER (kmmax = 10000, kmax=1500)
      REAL*8  delrgi,boxvol,facrdf,rcon,qr,gint,ginr,fns,gr,xcut
     &     ,sk(kmmax),gor(kmmax),gor1(kmmax),ssk,rk,rp,fact
      COMMON /rag1/ sk,gor,gor1

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      IF(kmmax .LT. maxint) THEN
         iret=1
         errmsg=
     &'In WRITE_GOFRW: Physical dimensions of G of R exceed limit of'
     &//' 4000. Abort.'
         RETURN
      END IF
      IF(ngrdon.LE.0) RETURN
      IF(rewind) REWIND kgofr_sk
      WRITE(kgofr_sk,'(''Time Step ='',f15.5)') fstep
      nts=0
      DO jsp=1,nspec
         nts=nts+itypes(jsp)
      END DO
      nts=nts*nmol
      xcut=gofr_cut**2
      delrgi=1.0D+0/delrg
      boxvol=volume
      facrdf=1.0D0/(4.0d0*pi*delrg)
      max=MIN0(IDINT(delrgi*DSQRT(xcut)+0.5d0),(maxint-1))
      vol=0.0D0
      DO j=1,max 
         rp=j*delrg
         vol=vol+rp*rp*delrg
      END DO
      vol=vol*4.0D0*pi
      
      DO 20 jsp=1,nspec
         numa=nmol*itypes(jsp)
         numb=numa
         rcon=2.0D0/DFLOAT(ngrdon*numa)
         qr=0.0d0
         gint=0.0d0
         lpoint=jsp*(jsp-1)/2+jsp
         WRITE(kgofr_sk,1000) beta_slv(jsp),beta_slv(jsp)
         DO 30 i=1,max
            qr=qr+delrg
            fns=rcon*DFLOAT(igr(i,lpoint+offset))
            gint=gint+fns
            gr=boxvol*facrdf*fns/(qr*qr*DFLOAT(numb))
            gor(i)=gr
            WRITE(kgofr_sk,500) qr,gr,gint
30       CONTINUE

         WRITE(kgofr_sk,1500) beta_slv(jsp),beta_slv(jsp)
         DO k=1,kmax
            rk=2.0D0*pi*DFLOAT(k)/kmax
            ssk=0.0D0
            DO j=1,max
               rp=j*delrg
               ssk=ssk+rp*rp*(gor(j)*DFLOAT(numb)/boxvol-(gint+1.0D0)
     &              /vol)*DSIN(rk*rp)*delrg/(rk*rp)
            END DO
            ssk=(1.0D0+4.0D0*pi*ssk)*DFLOAT(numa)/DFLOAT(nts)
            WRITE(kgofr_sk,700) rk,ssk
         END DO
20    CONTINUE
      itype=nspec
      IF(nspec.ne.1) THEN
         DO 40 isp=1,nspec-1
            numa=nmol*itypes(isp)
            rcon=1.0d0/DFLOAT(ngrdon*numa)
            DO 50 jsp=isp+1,nspec
               itype=itype+1
               numb=itypes(jsp)*nmol 
               WRITE(kgofr_sk,1000) beta_slv(isp),beta_slv(jsp)
               qr=0.0d0
               gint=0.0d0
               ginr=0.0d0
               lpoint=jsp*(jsp-1)/2+isp
               DO 60 i=1,max
                  qr=qr+delrg
                  fns=rcon*DFLOAT(igr(i,lpoint+offset))
                  gint=gint+fns
                  gr=boxvol*facrdf*fns/(DFLOAT(numb)*qr*qr)
                  gor(i)=gr
                  ginr=ginr+fns*DFLOAT(numa)/DFLOAT(numb)
                  WRITE(kgofr_sk,600) qr,gr,gint,ginr
60             CONTINUE
               WRITE(kgofr_sk,1500)beta_slv(isp),beta_slv(jsp)
               DO k=1,kmax
                  rk=2.0D0*pi*DFLOAT(k)/kmax
                  ssk=0.0D0
                  DO j=1,max
                     rp=j*delrg
                     ssk=ssk+rp*rp*(gor(j)*DFLOAT(numb)/boxvol-gint/vol)
     &                    *DSIN(rk*rp)*delrg/(rk*rp)
                  END DO
                  fact=DFLOAT(numa)/DFLOAT(nts)
                  ssk=4.0D0*pi*ssk*fact
                  WRITE(kgofr_sk,700) rk,ssk
               END DO
50          CONTINUE
40       CONTINUE
      END IF

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

500   FORMAT(3f15.6)
600   FORMAT(4f15.6)
700   FORMAT(2f15.8)
1000  FORMAT(/16x,a1,2x,1h-,2x,a1,5x,19hRadial Distribution,
     x     1x,8hFunction)
1500  FORMAT(/16x,a1,2x,1h-,2x,a1,5x,19h      S ( K )      ,
     x     1x,8hFunction)
      RETURN
      END
