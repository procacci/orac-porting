      subroutine lc_list(ncx,ncy,ncz,nind,indxi,indxj,indxk,rcut,co,xp
     &     ,yp,zp,natp,nnlpp,kprint)

************************************************************************
*  Written by Marc Souaille, CECAM 96. Compute linked cell list.       *
*  Cell indexing is done in lc_index.                                  *
************************************************************************

*======================= DECLARATIONS ==================================
      
      use parst
#ifdef _OMP_
      use omp_integr, only: nthr,mpp8,time_u,tcpu_u,omp_timing,array_omp
#ifdef _BGQ_  
      use omp_lib       ! OMP layer
#endif
#else
      use omp_integr, only: time_u,tcpu_u,omp_timing
#endif
      implicit none

      include 'lc_list.h' 
  
*----------------------- ARGUMENTS -------------------------------------

      integer ncx,ncy,ncz
      integer nind(2)
      integer indxi(2,*),indxj(2,*),indxk(2,*)
      integer natp

#ifdef _OMP_
      INTEGER  nnlpp(mpp8,nthr)
#ifndef _BGQ_  
      include 'omp_lib.h'       ! OMP layer
#endif
#else
      integer nnlpp(*)
#endif
      integer   kprint

      real*8 xp(*),yp(*),zp(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8  tela1,tela2,tcpu1,tcpu2
      integer treal
      character*80 errmsg 
      integer headw(CELLMAX),headp(CELLMAX)
      integer chainp(tgroup)
      integer cellpi(tgroup),cellpj(tgroup),cellpk(tgroup)
      integer nx,ny,nz
      integer i,j,k,l,m,n,nv,numcell,nppp
      integer iv,jv,kv,nmin
      integer itask,istride,dimmax,ndims
#ifdef _OMP_
      integer nvtot(nthr)
      real*8 nvtot_av,nvtot_av2
#else
      integer nvtot
#endif
      real*8 dx,dy,dz,co(3,3),rcut
      real*8 x1,y1,z1,x2,y2,z2,xx,yy,zz 
      real*8 sqcut,d

c--   local stuff is dumped in a scratch common 

      COMMON /rag1/ headw,headp,chainp,cellpi,cellpj,cellpk

 
      INCLUDE 'pbc.h'
   
*=======================================================================
*                      Initialization
*=======================================================================

      sqcut = rcut**2
      dx=2.d0/ncx
      dy=2.d0/ncy
      dz=2.d0/ncz

      do 10 n=1,ncx*ncy*ncz
         headp(n)=0
10    continue

*=======================================================================
*     Compute chain list for protein
*=======================================================================

        if(omp_timing) call timer(treal,tcpu1,tela1) 
#ifdef _OMP_
!!$OMP PARALLEL DO SCHEDULE(STATIC) 
#endif
      do n=1,natp
         x1=xp(n)/dx
         y1=yp(n)/dy
         z1=zp(n)/dz
         nx=int(x1)+(sign(1.d0,x1-int(x1))-1.)/2
         ny=int(y1)+(sign(1.d0,y1-int(y1))-1.)/2
         nz=int(z1)+(sign(1.d0,z1-int(z1))-1.)/2
         nx=mod(mod(nx,ncx)+ncx,ncx)
         ny=mod(mod(ny,ncy)+ncy,ncy)
         nz=mod(mod(nz,ncz)+ncz,ncz)
         cellpi(n)=nx
         cellpj(n)=ny
         cellpk(n)=nz
         numcell=nz+ncz*(ny+ncy*nx)+1
         chainp(n)=headp(numcell)
         headp(numcell)=n
       END DO
#ifdef _OMP_
!!$OMP END PARALLEL DO 
#endif

*=======================================================================
*     Compute neighbor list nnlpp protein protein
*=======================================================================

      nppp=0
#ifdef _OMP_
!$OMP  PARALLEL DEFAULT(PRIVATE) 
!$OMP& SHARED(xp,yp,zp,co,natp,indxi,indxj,indxk,ncx,ncy,ncz,mpp8)
!$OMP& SHARED(array_omp)
!$OMP& SHARED(cellpi,cellpj,cellpk,chainp,headp,nind,sqcut,nnlpp,nvtot)
!$OMP& REDUCTION(+:nppp)
      dimmax=mpp8
      itask=1+OMP_GET_THREAD_NUM()
      nvtot(itask)=1
#else
      dimmax=mpp
      nppp=0
      nvtot=1
#endif
      do 60 n=1,natp
#ifdef _OMP_
         if(array_omp(n,itask)) THEN
#endif
         x1=xp(n)
         y1=yp(n)
         z1=zp(n)
         nv=0         
         i=cellpi(n)
         j=cellpj(n)
         k=cellpk(n)
         do 70 m=1,nind(2)
            iv=indxi(2,m)
            jv=indxj(2,m)
            kv=indxk(2,m)
            nx=mod(mod(i+iv,ncx)+ncx,ncx)
            ny=mod(mod(j+jv,ncy)+ncy,ncy)
            nz=mod(mod(k+kv,ncz)+ncz,ncz)
            numcell=nz+ncz*(ny+ncy*nx)+1
            l=headp(numcell)
            nmin=0
            if(m.eq.1) nmin=n
            do while(l.gt.nmin)
               if(l.gt.n) then
                  x2=x1-xp(l)
                  y2=y1-yp(l)
                  z2=z1-zp(l)
                  x2=x2-2.0*pbc(x2)
                  y2=y2-2.0*pbc(y2)
                  z2=z2-2.0*pbc(z2)
                  xx=co(1,1)*x2+co(1,2)*y2+co(1,3)*z2
                  yy=co(2,1)*x2+co(2,2)*y2+co(2,3)*z2
                  zz=co(3,1)*x2+co(3,2)*y2+co(3,3)*z2
                  d=xx**2+yy**2+zz**2
                  if(d.lt.sqcut) then
                     nppp=nppp+1
                     nv=nv+1
#ifdef _OMP_
                     nnlpp(nvtot(itask)+nv,itask)=l
#else
                     nnlpp(nvtot+nv)=l
#endif
                  endif
               endif
               l=chainp(l)
            end do
 70      continue         


#ifdef _OMP_
         nnlpp(nvtot(itask),itask)=nv
         nvtot(itask)=nvtot(itask)+nv+1
         ndims=nvtot(itask)
#else
         nnlpp(nvtot)=nv
         nvtot=nvtot+nv+1
         ndims=nvtot
#endif
         if(ndims.gt.dimmax) then
            errmsg=' IN LC_LIST : Dimensions of the NNLPP are 
     &           insufficient.ABORT ** '
            WRITE(kprint,'('' MAX ='',i7,'' NCOUNT = '',32i7)') dimmax
     &           ,nvtot
            CALL xerror(errmsg,80,1,2)
         endif
#ifdef _OMP_
         ENDIF
#endif
60    continue
#ifdef _OMP_
!$OMP END PARALLEL 
#endif

#ifdef _OMP_
      WRITE(kprint,*)
      nvtot_av=0.d0
      nvtot_av2=1.D-8
      DO I=1,nthr
        nvtot_av=nvtot_av+nvtot(i)
        nvtot_av2=nvtot_av2+nvtot(i)**2.
      end do
      nvtot_av=nvtot_av/float(nthr)
      nvtot_av2= nvtot_av2/float(nthr)
      IF(nthr.GT.1) THEN
        WRITE(kprint,10001) nvtot_av,dsqrt(nvtot_av2-nvtot_av**2)
      ELSE
        WRITE(kprint,10000) nint(nvtot_av)
      END IF
#else
      WRITE(kprint,*)
      WRITE(kprint,10000) nppp+natp
#endif      
*================= END OF EXECUTABLE STATEMENTS ========================

10000 FORMAT(5x,'Neighbor Lists Dimensions *neighbor(',i8,')* ')
10001 FORMAT(5x,'Average Neighbor Lists Dimension -> ',F15.1,' +/-'
     &     ,F15.1)

      if(omp_timing) THEN 
        call timer(treal,tcpu2,tela2)
        time_u=time_u+tela2-tela1
        tcpu_u=tcpu_u+tcpu2-tcpu1
      end if
      RETURN
      END
