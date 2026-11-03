c=======================================================================
      SUBROUTINE field_chromo(xp,yp,zp,vi,cut,tkvi)
c=======================================================================
c   Transit Routine: calls field_intra field_pme field_direct
C   Drive routine for field calculation. Computes field at chromophores
c   specified in input READ_SOLUTE by DEF_COFACTOR directive (see ORAC
c   manual). Written by P PROCACCI CECAM SEPT 96
c======================================================================
c    FIELD routines: compute the potential at each atom 
c    of selected chromphores. Intrachromophores contribution 
c    are *not* included in the potential.
c           =============================   
c    routine FIELD_PME:     compute PME contribution
c    routine FIELD_DIRECT:  compute direct lattice contribution
c    routine FIELD_INTRA:   compute intrachromophores contrib. 
c    routine FIELD_DUMP:    dump potential and coordinates onto disk
c    routine FIELD_INIT:    print out segments  
c    routine FIELD_CHROM:   transit routine called by mtsmd/drvmda
c
c=======================================================================

c======================= DECLARATIONS ==================================

      use unit
      use parst
      use cpropar

      IMPLICIT NONE

c----------------------- ARGUMENTS -------------------------------------

      REAL*8  vi(*),xp(*),yp(*),zp(*),cut,tkvi,pir

c-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,k,ibeg,iend,n1t,n2t
      REAL*8  x0(m1),y0(m1),z0(m1),xp0(m1),yp0(m1),zp0(m1)

!==================== EXECUTABLE STATEMENTS ============================

c--------------find out solute and solvent coordinates ----------------
      j=0 
      k=0 
      do i=1,ntap
        if(ss_index(i).eq.1) then 
          j = j+1
          xp0(j)=xp(i)
          yp0(j)=yp(i)
          zp0(j)=zp(i)
        else if (ss_index(i).eq.2) THEN 
          k = k+1 
          x0(k)= xp(i) 
          y0(k)= yp(i) 
          z0(k)= zp(i) 
        end if
      end do
      n1t = j
      n2t = k
      if(ncofactor.eq.ntap) ncofactor=n1t
        
c--   computes reciprocal lattice contribution to field 
      If(clewld) 
     &CALL field_pme(vi,n2t,n1t,x0,y0,z0,xp0,yp0,zp0,chrge,co
     &     ,oc,volume,alphaf,pme_order,nfft1,nfft2,nfft3,slt_exist
     &     ,slv_exist,ncofactor,nsegs)

c--   computes direct lattice contribution to field  and
c---  subtract intrachromophore electrostatic energy. 
      pir = 2.0/dsqrt(3.14159265358979323846d0) 
      DO i =1,ncofactor
         ibeg = nsegs(1,i)
         iend = nsegs(2,i)
         CALL field_direct(xp0,yp0,zp0,x0,y0,z0,vi,n1t,n2t,chrge,co
     &        ,cut,alphaf,clewld,protei,slv_exist,ibeg,iend)
         If(clewld) THEN 
            CALL field_intra(ibeg,iend,chrge,co,xp0,yp0,zp0,vi,alphaf)
            do j=ibeg,iend
               vi(j) = vi(j)-alphaf*chrge(j)*pir
            end do
         END IF
      END DO

c-   dump potential onto disk

      CALL field_dump(tkvi,ncofactor,nsegs,co,xp0,yp0,zp0,chrge,vi
     &     ,lform_field)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END 

c======================================================================
      subroutine field_pme(vi,nts,ntap,xs,ys,zs,xp,yp,zp,pmechg
     &     ,co,oc,volume,alphal,order,nfft1,nfft2,nfft3,pflag,sflag
     &     ,ncofactor,nsegs)
c=======================================================================
c   Computes vi PME reciprocal potential at positions xi yi zi.   
c=======================================================================

      use spme
      IMPLICIT none

#ifdef FFTW
      include "fftw3.f"
#endif

c---  stuff for b-spline interpolation and FFT 

      INTEGER   NFFT1,NFFT2,NFFT3,NUMATOMS,ORDER,i,j,NCOFACTOR

#ifdef FFTW
      real*8 Q(nfft1,nfft2,nfft3),theta1(nts*ntap*order), theta2(nts
     &     *ntap*order),theta3(nts*ntap*order),dtheta1(nts*ntap*order)
     &     ,dtheta2(nts*ntap*order),dtheta3(nts*ntap*order),fr1(nts
     &     *ntap), fr2(nts*ntap),fr3(nts*ntap),recip(3,3),x(max_atm)
     &     ,y(max_atm),z(max_atm),fx(max_atm),fy(max_atm),fz(max_atm)
      INTEGER nsegs(2,*)
      complex*16 FQ(nfft1,nfft2,(nfft3/2+1))
#else
      REAL*8  theta1(mth),dtheta1(mth),theta2(mth),dtheta2(mth)
     &     ,theta3(mth),dtheta3(mth)
      REAL*8  Q(maxt),fr1(max_atm),fr2(max_atm),fr3(max_atm),recip(3,3)
     &     ,x(max_atm),y(max_atm),z(max_atm),fx(max_atm),fy(max_atm)
     &     ,fz(max_atm)
      REAL*8 sx,sy,sz,stc1,stc2,stc3,stc4,stc5,stc6,stc7,stc8,stc9
      INTEGER k,nta,nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft
     &     ,sffw,nsegs(2,*)
#endif      

c---  non local stuff 
      
      INTEGER  nts,ntap
      REAL*8   xs(*),ys(*),zs(*),xp(*),yp(*),zp(*),pmechg(*),vi(*)
      REAL*8   co(3,3),oc(3,3),volume,alphal,stressc(3,3),eer
      LOGICAL pflag,sflag

c---  transform to cartesian coordinates

      IF(pflag) THEN
         do i=1,ntap
            vi(i)=0.d0
            x(i)=co(1,1)*xp(i)+co(1,2)*yp(i)+co(1,3)*zp(i)
            y(i)=co(2,1)*xp(i)+co(2,2)*yp(i)+co(2,3)*zp(i)
            z(i)=co(3,1)*xp(i)+co(3,2)*yp(i)+co(3,3)*zp(i)
         end do
      END IF
      IF(sflag) THEN
         do i=1,nts
            x(i+ntap)=co(1,1)*xs(i)+co(1,2)*ys(i)+co(1,3)*zs(i)
            y(i+ntap)=co(2,1)*xs(i)+co(2,2)*ys(i)+co(2,3)*zs(i)
            z(i+ntap)=co(3,1)*xs(i)+co(3,2)*ys(i)+co(3,3)*zs(i)
         end do
      END IF
      numatoms=nts+ntap

c---  PME wants the transpose of our oc matrix

      do i=1,3
         do j=1,3
            recip(i,j) = 0.5*oc(j,i)
         end do
      end do

c---  call all PME routines 
#ifdef FFTW
      siztheta=numatoms*order
      call get_scaled_fractionals(numatoms,x,y,z,recip,nfft1,nfft2,nfft3
     &     ,fr1,fr2,fr3)
      call get_bspline_coeffs(numatoms,fr1,fr2,fr3,order,theta1,theta2
     &     ,theta3,dtheta1,dtheta2,dtheta3)
      call fill_charge_grid(numatoms,pmechg,theta1,theta2,theta3,fr1,fr2
     &     ,fr3,order,nfft1,nfft2,nfft3,Q)
      call dfftw_execute_dft_r2c(planf, Q, FQ)
      call field_scalar_sum(FQ,alphal,volume,recip,bsp_mod1
     &     ,bsp_mod2,bsp_mod3,nfft1,nfft2,nfft3)
      call dfftw_execute_dft_c2r(planr, FQ, Q)
      call field_sum(numatoms,nts,vi,nsegs,ncofactor,theta1,theta2
     &     ,theta3,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,Q)
#else
      call get_fftdims(nfft1,nfft2,nfft3,
     $       nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)
      call get_scaled_fractionals(numatoms,x,y,z,recip,nfft1,nfft2,nfft3
     &     ,fr1,fr2,fr3)
      call get_bspline_coeffs(numatoms,fr1,fr2,fr3,order,theta1,theta2
     &     ,theta3,dtheta1,dtheta2,dtheta3)
      call fill_charge_grid(numatoms,pmechg,theta1,theta2,theta3,fr1,fr2
     &     ,fr3,order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q)
      call fft_back(Q,fftable,ffwork,nfft1,nfft2,nfft3,nfftdim1,nfftdim2
     &     ,nfftdim3,nfftable,nffwork)
      call field_scalar_sum(Q,alphal,volume,recip,bsp_mod1,bsp_mod2
     &     ,bsp_mod3,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3)
      call fft_forward(Q,fftable,ffwork,nfft1,nfft2,nfft3,nfftdim1
     &     ,nfftdim2,nfftdim3,nfftable,nffwork)
      call field_sum(numatoms,nts,vi,nsegs,ncofactor,theta1,theta2
     &     ,theta3,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2
     &     ,nfftdim3,Q)
#endif

      return
      end 

c======================================================================
#ifdef FFTW
      SUBROUTINE field_sum(numatoms,nts,vi,nsegs,ncofactor,theta1,theta2
     &     ,theta3,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,Q)
#else
      SUBROUTINE field_sum(numatoms,nts,vi,nsegs,ncofactor,theta1,theta2
     &     ,theta3,fr1,fr2,fr3,order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2
     &     ,nfftdim3,Q)
#endif
c======================= DECLARATIONS ==================================

      IMPLICIT NONE
      integer order,nfft1,nfft2,nfft3,nsegs(2,*),ncofactor,numatoms
      REAL*8 fr1(numatoms),fr2(numatoms),fr3(numatoms)
      REAL*8 theta1(order,numatoms),theta2(order,numatoms),
     $     theta3(order,numatoms),vi(numatoms)
#ifdef FFTW
      REAL*8 Q(nfft1,nfft2,nfft3)
#else
      REAL*8 Q(2,nfftdim1,nfftdim2,nfftdim3)
      integer nfftdim1,nfftdim2,nfftdim3
#endif

      integer n1,n,ith1,ith2,ith3,i0,j0,k0,i,j,k,ic,nts
      REAL*8 f1,f2,f3,term

!==================== EXECUTABLE STATEMENTS ============================

      do ic=1,ncofactor
         do n1 = nsegs(1,ic),nsegs(2,ic)
            n = n1
            f1 = 0.d0
            k0 = int(fr3(n)) - order
            do 200 ith3 = 1,order
               k0 = k0 + 1
               k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
               j0 = int(fr2(n)) - order
               do 150 ith2 = 1,order
                  j0 = j0 + 1
                  j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
                  i0 = int(fr1(n)) - order
                  do 100 ith1 = 1,order
                     i0 = i0 + 1
                     i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
#ifdef FFTW
                     term = Q(i,j,k)
#else
                     term = Q(1,i,j,k)
#endif
                     f1 = f1+term*theta1(ith1,n)*theta2(ith2,n)
     &                    *theta3(ith3,n)
 100              continue
 150           continue
 200        continue
            vi(n1) = f1
 400        continue
         end do
      end do
      
*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END

c======================================================================
      SUBROUTINE field_intra(n1,n2,charge,co,x0,y0,z0,vi,alphal)
c======================================================================
c     This subroutine sutract the Ewald intra-term from the Vi potents 
c----------------------------------------------------------------------

c======================= DECLARATIONS ==================================

      IMPLICIT NONE

c----------------------- ARGUMENTS -------------------------------------

      INTEGER n1,n2
      REAL*8  charge(*),x0(*),y0(*),z0(*),vi(*),co(3,3),t1

c-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER k,ia,ib
      REAL*8  xab,yab,zab,dx,dy,dz,rsq,rsp,alphal
      REAL*8  a1,a2,a3,a4,a5,qp,qt,expcst,erfcst,erfst,alphar,pir
      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/

!==================== EXECUTABLE STATEMENTS ============================

      if(n1.eq.n2) return
         
      DO ia=n1,n2-1
         DO ib=ia+1,n2
            xab=x0(ia)-x0(ib)
            yab=y0(ia)-y0(ib)
            zab=z0(ia)-z0(ib)
            dx=co(1,1)*xab+co(1,2)*yab+co(1,3)*zab
            dy=co(2,1)*xab+co(2,2)*yab+co(2,3)*zab
            dz=co(3,1)*xab+co(3,2)*yab+co(3,3)*zab
            rsq=dx*dx+dy*dy+dz*dz
            rsp=DSQRT(rsq)
            alphar=alphal*rsp
            qt=1.0D0/(1.0e0+qp*alphar)
            expcst=exp(-alphar*alphar)
            erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*
     x           expcst
            erfst=1.0D0-erfcst
            vi(ia)=vi(ia)-erfst*charge(ib)/rsp
            vi(ib)=vi(ib)-erfst*charge(ia)/rsp
         END DO
      END DO 
      
*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END

c======================================================================
      SUBROUTINE field_init(nsegs,segs,kprint,j) 
c======================================================================
c     This subroutine initializes coordinates and intra list arrays    
c     only atoms defined as solute are considered. 
c----------------------------------------------------------------------

c======================= DECLARATIONS ==================================

      IMPLICIT NONE

c----------------------- ARGUMENTS -------------------------------------

      INTEGER   nsegs,segs(2,*),kprint
  
c-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER   i,k,j,ibeg,iend

!==================== EXECUTABLE STATEMENTS ============================


c---  fill array of i-atoms using input segments 

      j=0
      write(kprint,50) 
 50   format(
     &//'  ====================================================='/
     & ,'  = V computed for the following atoms of the solute: ='/
     & ,'  ====================================================='/
     & ,'     segment       from atom             to atom       '/)

      DO i =1,nsegs
         ibeg = segs(1,i)
         iend = segs(2,i)
         write(kprint,10) i,ibeg,iend
 10      format(2x,i7,4x,i10,13x,i10)
         DO k=ibeg,iend
            j=j+1
         END DO
      END DO
      write(kprint,70) j
 70   format(/' Total number of atoms =',i10/)
      
*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END

c======================================================================
      SUBROUTINE field_dump(time,nsegs,segs,co,x,y,z,qq,vi,lf)
c======================================================================
c     This subroutine dumps onto disk Vi potential and coordinates
c----------------------------------------------------------------------

c======================= DECLARATIONS ==================================

      use unit

      IMPLICIT NONE

c----------------------- ARGUMENTS -------------------------------------

      INTEGER   nsegs,segs(2,*)
      REAL*8    x(*),y(*),z(*),qq(*),co(3,3),vi(*),time
      logical   lf

c-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER   i,l,nat_seg,ibeg,iend,k
      real*4    vi4,x4,y4,z4,q4

!==================== EXECUTABLE STATEMENTS ============================

      if(.not.lf) THEN 
c---    write unformatted file 
        write(kvi) time
        DO i =1,nsegs
          ibeg = segs(1,i)
          iend = segs(2,i)
          DO k=ibeg,iend
            q4=dsqrt(unitc)*qq(k)
            vi4= efact*vi(k)/1000.d0/dsqrt(unitc)
            x4=co(1,1)*x(k)+co(1,2)*y(k)+co(1,3)*z(k)
            y4=co(2,1)*x(k)+co(2,2)*y(k)+co(2,3)*z(k)
            z4=co(3,1)*x(k)+co(3,2)*y(k)+co(3,3)*z(k)
            write(kvi) vi4,x4,y4,z4,q4
          END DO
        END DO
c---    write formatted file 
      ELSE
        write(kvi,10) time
10      FORMAT(f20.4) 
        DO i =1,nsegs
          ibeg = segs(1,i)
          iend = segs(2,i)
          DO k=ibeg,iend
            vi4= efact*vi(k)/1000.d0/dsqrt(unitc)
            q4=dsqrt(unitc)*qq(k)
            x4=co(1,1)*x(k)+co(1,2)*y(k)+co(1,3)*z(k)
            y4=co(2,1)*x(k)+co(2,2)*y(k)+co(2,3)*z(k)
            z4=co(3,1)*x(k)+co(3,2)*y(k)+co(3,3)*z(k)
            write(kvi,20) vi4,x4,y4,z4,q4
20          FORMAT(F20.5,3f12.4,f9.5) 
          END DO
        END DO
      ENDIF
*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
      


c======================================================================
      SUBROUTINE field_direct(xp,yp,zp,xs,ys,zs,vi,ntap,nts,pmechg,co
     &     ,cut,alphal,ewald,pflag,sflag,ibeg,iend)
c======================================================================
c     This subroutine computes the direct part contr. to Vi
c----------------------------------------------------------------------

c======================= DECLARATIONS ==================================

      IMPLICIT NONE

      INTEGER  ntap,nts,i,j,k,ibeg,iend
      REAL*8   xp(*),yp(*),zp(*),xs(*),ys(*),zs(*),pmechg(*),vi(*),cut
     &     ,rsp
      REAL*8   co(3,3),alphal,xc,yc,zc,rc2,xab,yab,zab
      LOGICAL pflag,sflag,ewald
      REAL*8  a1,a2,a3,a4,a5,qp,qt,expcst,erfcst,erfst,alphar
      DATA a1,a2,a3/0.2548296d0,-0.28449674d0,1.4214137d0/
      DATA a4,a5/-1.453152d0,1.0614054d0/
      DATA qp/0.3275911d0/

c----------------------- ARGUMENTS -------------------------------------

      INCLUDE 'pbc.h'

!==================== EXECUTABLE STATEMENTS ============================


c--   computes cofactor-rest_of_the_world direct space interaction 

      do i = ibeg,iend
         do j=1,nts
            xab=xp(i)-xs(j)
            yab=yp(i)-ys(j)
            zab=zp(i)-zs(j)
            xab=xab-2.0*PBC(xab)
            yab=yab-2.0*PBC(yab)
            zab=zab-2.0*PBC(zab)
            xc=co(1,1)*xab+co(1,2)*yab+co(1,3)*zab
            yc=co(2,1)*xab+co(2,2)*yab+co(2,3)*zab
            zc=co(3,1)*xab+co(3,2)*yab+co(3,3)*zab
            rc2=xc*xc+yc*yc+zc*zc
            if(rc2.gt.cut) goto 1001
            rsp = dsqrt(rc2)
            if(ewald) then
               alphar=alphal*rsp
               qt=1.0d0/(1.0d0+qp*alphar)
               expcst=DEXP(-alphar*alphar)
               erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)
     x              *qt+a1)*qt*expcst
            else
               erfcst=1.d0
            END IF
            vi(i) = vi(i) + pmechg(j+ntap)*erfcst/rsp
 1001       continue
         end do
         do j=1,ntap
            if(j.ge.ibeg.and.j.le.iend) go to 2001
            xab=xp(i)-xp(j)
            yab=yp(i)-yp(j)
            zab=zp(i)-zp(j)
            xab=xab-2.0*PBC(xab)
            yab=yab-2.0*PBC(yab)
            zab=zab-2.0*PBC(zab)
            xc=co(1,1)*xab+co(1,2)*yab+co(1,3)*zab
            yc=co(2,1)*xab+co(2,2)*yab+co(2,3)*zab
            zc=co(3,1)*xab+co(3,2)*yab+co(3,3)*zab
            rc2=xc*xc+yc*yc+zc*zc
            if(rc2.gt.cut) goto 2001
            rsp = dsqrt(rc2)
            if(ewald) then
               alphar=alphal*rsp
               qt=1.0d0/(1.0d0+qp*alphar)
               expcst=DEXP(-alphar*alphar)
               erfcst=((((a5*qt+a4)*qt+a3)*qt+a2)
     x              *qt+a1)*qt*expcst
            else
               erfcst=1.d0
            END IF
            vi(i) = vi(i) + pmechg(j)*erfcst/rsp
 2001       continue
         end do
      end do

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END 

c======================================================================
#ifdef FFTW
      subroutine field_scalar_sum(
     $     Q,ewaldcof,volume,recip,bsp_mod1,bsp_mod2,bsp_mod3,nfft1
     &     ,nfft2,nfft3)
      implicit none
      integer nfft1,nfft2,nfft3
      complex*16 Q(1+nfft1/2,nfft2,nfft3)
      REAL*8 bsp_mod1(nfft1),bsp_mod2(nfft2),
     +   bsp_mod3(nfft3),ewaldcof,volume,rkcut
      REAL*8 eer,vir(3,3)
      REAL*8 recip(3,3)
      REAL*8 pi,fac,denom,eterm
      integer k,k1,k2,k3,m1,m2,m3,nff,ind,jnd,indtop
      integer nf1,nf2,nf3,ndim1
      REAL*8 mhat1,mhat2,mhat3,msq,struc2,rkcut2,fac12
      ndim1=1+nfft1/2
      indtop = ndim1*nfft2*nfft3
      pi = 3.14159265358979323846
      rkcut2=0.25*rkcut*rkcut/pi**2
      fac = pi**2/ewaldcof**2
      nff = ndim1*nfft2
      nf1 = nfft1/2
      if ( 2*nf1 .lt. nfft1 )nf1 = nf1+1
      nf2 = nfft2/2
      if ( 2*nf2 .lt. nfft2 )nf2 = nf2+1
      nf3 = nfft3/2
      if ( 2*nf3 .lt. nfft3 )nf3 = nf3+1
      do 100 ind = 1,indtop-1
c get k1,k2,k3 from the relationship
c           ind = (k1-1) + (k2-1)*nfft1 + (k3-1)*nfft2*nfft1
         k3 = ind/nff + 1
         jnd = ind - (k3-1)*nff
         k2 = jnd/ndim1 + 1
         k1 = jnd - (k2-1)*ndim1 +1
         m1 = k1 - 1
         if ( k1 .gt. nf1 )m1 = k1 - 1 - nfft1
         m2 = k2 - 1
         if ( k2 .gt. nf2 )m2 = k2 - 1 - nfft2
         m3 = k3 - 1
         if ( k3 .gt. nf3 )m3 = k3 - 1 - nfft3
         mhat1 = recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3
         mhat2 = recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3
         mhat3 = recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3
         msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
         IF(msq .GT. rkcut2) THEN 
            eterm = 0.d0 
         ELSE
            eterm = dexp(-fac*msq)/denom
         END IF
!        B,C and F[Q] are in reciprocal space; H_rec and Q are in direct space.
!        this BC*F[Q] = F[H_rec]*F[Q] where H_rec=F^-1[BC]  
!        so that F^-1(F[H]*F(Q)) =H_rec x Q (where "x" indicate convolution) 
1098     Q(k1,k2,k3) = eterm * Q(k1,k2,k3) ! array is overwritten by the product BC*F[Q]
100   continue
      return
      end
#else
      subroutine field_scalar_sum(
     $         Q,ewaldcof,volume,recip,bsp_mod1,bsp_mod2,bsp_mod3,
     $     nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3)
      implicit none
      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
      REAL*8 Q(2,nfftdim1,nfftdim2,nfftdim3)
      REAL*8 bsp_mod1(nfft1),bsp_mod2(nfft2),
     +   bsp_mod3(nfft3),ewaldcof,volume,rkcut
      REAL*8 recip(3,3)
      REAL*8 pi,fac,denom,eterm
      integer k,k1,k2,k3,m1,m2,m3,nff,ind,jnd,indtop
      integer nf1,nf2,nf3
      REAL*8 mhat1,mhat2,mhat3,msq,struc2,rkcut2
      indtop = nfft1*nfft2*nfft3
      pi = 3.14159265358979323846
      rkcut2=0.25*rkcut*rkcut/pi**2
      fac = pi**2/ewaldcof**2
      nff = nfft1*nfft2
      nf1 = nfft1/2
      if ( 2*nf1 .lt. nfft1 )nf1 = nf1+1
      nf2 = nfft2/2
      if ( 2*nf2 .lt. nfft2 )nf2 = nf2+1
      nf3 = nfft3/2
      if ( 2*nf3 .lt. nfft3 )nf3 = nf3+1
      do 100 ind = 1,indtop-1
c get k1,k2,k3 from the relationship
c           ind = (k1-1) + (k2-1)*nfft1 + (k3-1)*nfft2*nfft1
         k3 = ind/nff + 1
         jnd = ind - (k3-1)*nff
         k2 = jnd/nfft1 + 1
         k1 = jnd - (k2-1)*nfft1 +1
         m1 = k1 - 1
         if ( k1 .gt. nf1 )m1 = k1 - 1 - nfft1
         m2 = k2 - 1
         if ( k2 .gt. nf2 )m2 = k2 - 1 - nfft2
         m3 = k3 - 1
         if ( k3 .gt. nf3 )m3 = k3 - 1 - nfft3
         mhat1 = recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3
         mhat2 = recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3
         mhat3 = recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3
         msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
         IF(msq .GT. rkcut2) THEN 
            eterm = 0.d0 
         ELSE
            denom = pi*volume*bsp_mod1(k1)*bsp_mod2(k2)*bsp_mod3(k3)*msq
            eterm = dexp(-fac*msq)/denom
         END IF
1098     Q(1,k1,k2,k3) = eterm * Q(1,k1,k2,k3)
         Q(2,k1,k2,k3) = eterm * Q(2,k1,k2,k3)
100   continue
      return
      end
#endif

c======================================================================
      subroutine field_restart(kprint,kvi,mode,time,nsegs,segs,lf)
c======================================================================
*     Rudimental restart routine for field_chromo
c======================================================================

      implicit none

      integer  kprint,kvi,nsegs,segs(2,*),nsegs_old,seg1,seg2,i
      character*5 mode
      real*8  time,time_old
      logical lf

c--------------------------------------------------------------------------
c      close kvi unit and put a EOF
c--------------------------------------------------------------------------

      if(mode.eq.'write') then
c---     write end-of-file and close all units
         endfile(kvi)
         rewind(kvi)
         if(.not.lf) THEN 
5          read(kvi,end=10) 
           go to 5
         ELSE
6          read(kvi,*,end=10) 
           go to 6
         END IF
 10      continue
      end if

c--------------------------------------------------------------------------
c      check restart file and read kvi up to EOF
c--------------------------------------------------------------------------

      if(mode.eq.'read ') then
         rewind(kvi)
         if(lf) goto 900 
         read(kvi) time_old,nsegs_old 
         if(time_old.ne.time) then
            write(kprint,100) time_old,time
 100        format(
     &      //'ERROR: Time interval for cofactor does not match:'/
     &           ' old time = ',f10.4, 'current time = ',f10.4//)  
            STOP
         END IF
         if(nsegs_old.ne.nsegs) then
            write(kprint,200) nsegs_old,nsegs
 200        format(
     &      //'ERROR: Number of def. cofactors does not match:'/
     &           ' old N  = ',I5, 'current N =  ',I5//)  
            STOP
         END IF
         do i=1,nsegs
            read(kvi) seg1,seg2
            if ((seg1.ne.segs(1,i)).and.(seg2.ne.segs(2,i))) THEN
               write(kprint,300) i,seg1,seg2,segs(1,i),segs(2,i)
 300           format(
     &              //'ERROR: Number of atoms in',i3, 
     &              'cofactor does not match'/
     &              ' old end atoms  = ',2I5, 'current end atoms = ',2I5
     &              //)  
               STOP
            END IF
         end do
c--      all tests are passed: now can reposition the tape of to the EOF
900      if(.not.lf) THEN 
15         read(kvi,end=110) 
           go to 15
         ELSE
16         read(kvi,*,end=110) 
           go to 16
         END IF
110      CONTINUE
      ENDIF
      RETURN
      END
