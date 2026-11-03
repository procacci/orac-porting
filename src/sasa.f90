subroutine force_sasa(n,nbtype,betb,x,y,z,co,r2cut,rsolv,rvdw,sp,p,enersolv,fpx,fpy,fpz,sasa) 
  implicit none
  

! arguments (In)
  real*8 x(*),y(*),z(*),co(3,3)     ! coordinates in the crystal frame and co matrix
  real*8 rvdw(*)                    ! van Der Waals radii
  real*8 sp(*)                      ! free energy per unit area (kcal/mole/angs) 
  real*8 p(*)                       ! atomic SASA parameters (reported in table  1  of Hasel, Tetr. Comput  Methodol 1,103, 1988
  real*8 r2cut                      ! cut off (2*max(rvdv)+rsolv)
  real*8 rsolv                      ! solvent radius
  integer n                         ! number of atoms
  integer nbtype(*)                 ! atomic types index
  character*7  betb(*)              ! atomic types, e.g. ct, c3, hc ca, etc               
! arguments (Out)              
  real*8 fpx(*),fpy(*),fpz(*)  ! atomic forces with SASA contr added
  real*8 enersolv              ! SASA energy  
  real*8 sasa(n)               ! atomic SASA 

! local variables 
  real*8           pi,r2,factbi,factbj,factj,si,sj,xc,yc,zc,xij,yij,zij,pbcx,pbcy,pbcz,r,rth,r1
  real*8           bij,bji,dbij,dbji,bij1,dbij1,bij2,dbij2,bji2,dbji2,pij,hij,p_i,p_j,rvdwi,rvdwj,sasatot
  character*1      atoi,atoj
  integer          i,j,k,index0(n*50),map,ncount
  real*8           rr(n*50),ppij(n*50),xijk(n*50),yijk(n*50),zijk(n*50),qforce(n),qforcej,spj
  real*8           qforcetot
  logical cnoi,cnoj

  sasatot=0.d0
  pi=dacos(-1.d0) 
  enersolv=0.d0
  ncount=0

  do i=1,n
     rvdwi=rvdw(nbtype(i))
     p_i=p(nbtype(i))
     si=4.d0*pi*(rvdwi+rsolv)**2. 
     sasa(i)=si
     atoi=betb(i)
     cnoi=.false.
     cnoi=atoi.eq."c".or.atoi.eq."n".or.atoi.eq."o".or.atoi.eq."h".or.atoi.eq."f"
     factbi=pi*(rvdwi+rsolv)
     map=0
     do j=1,n
        if(i.eq.j) go to 100
        atoj=betb(j)
        cnoj=.false.
        cnoj=atoj.eq."c".or.atoj.eq."n".or.atoj.eq."o".or.atoj.eq."h".or.atoj.eq."f"
        xc=x(i)-x(j)
        yc=y(i)-y(j)
        zc=z(i)-z(j)
        pbcx=2.0D0*DNINT(0.5D0*xc)
        pbcy=2.0D0*DNINT(0.5D0*yc)
        pbcz=2.0D0*DNINT(0.5D0*zc)
        xc=xc-pbcx        
        yc=yc-pbcy
        zc=zc-pbcz
        xij=co(1,1)*xc+co(1,2)*yc+co(1,3)*zc
        yij=co(2,1)*xc+co(2,2)*yc+co(2,3)*zc
        zij=co(3,1)*xc+co(3,2)*yc+co(3,3)*zc
        r2=xij**2+yij**2+zij**2
        if(r2.lt.r2cut) THEN 
           r=dsqrt(xij**2+yij**2+zij**2)
           rvdwj=rvdw(nbtype(j))
           if(cnoj.and.cnoi) THEN 
              if(r.lt.1.6) THEN 
                 pij=0.8875     ! Hasel, Tetr. Comput  Methodol 1,103, 1988
              else 
                 pij=0.3516
              end if
           else 
              if(r.lt.2.0) THEN 
                 pij=0.8875    ! Hasel, Tetr. Comput  Methodol 1,103, 1988
              else 
                 pij=0.3516
              end if
           end if
           rth=rvdwi+rvdwj+rsolv
           if(r.gt.rth) THEN 
              bij=0.d0 
              dbij=0.d0
           else 
              map=map+1
              index0(1+ncount+map)=j
              rr(1+ncount+map)=r
              xijk(1+ncount+map)=xij
              yijk(1+ncount+map)=yij
              zijk(1+ncount+map)=zij
              ppij(1+ncount+map)=pij
              bij1=(rvdwi+rvdwj+2*rsolv-r)
              bij2=(1+(rvdwj-rvdwi)/r)
              bij=factbi*bij1*bij2
              factj=1-p_i*pij*bij/si
              sasa(i)=sasa(i)*factj
           end if
        end if
100     continue
     end do
     index0(1+ncount)=map
     ncount=ncount + (map+1)
     enersolv=enersolv+sasa(i)*sp(nbtype(i))
     sasatot=sasatot+sasa(i)
!    compute forces
     qforce(i)=sp(nbtype(i))*sasa(i)*p_i/si
  end do


  ncount=0
  do i=1,n
     rvdwi=rvdw(nbtype(i))
     map=index0(ncount+1)
     do k=1,map
        j=index0(ncount+1+k)
        spj=sp(nbtype(j))
        p_j=p(nbtype(j))
        rvdwj=rvdw(nbtype(j))
        r=rr(ncount+1+k)
        r1=1/r
        xij=xijk(ncount+1+k)
        yij=yijk(ncount+1+k)
        zij=zijk(ncount+1+k)
        pij=ppij(ncount+1+k)
        factbj=pi*(rvdwj+rsolv)
        sj=4.d0*pi*(rvdwj+rsolv)**2. 
        qforcej=spj*sasa(j)*p_j*pij/sj
        bij1=(rvdwi+rvdwj+2*rsolv-r)
        hij=(rvdwj-rvdwi)*r1
        bij2 = 1+hij
        bji2 = 1-hij
        dbij1= -1.d0
        dbij2 =-hij*r1
        dbji2 = hij*r1
        bij=factbi*bij1*bij2
        bji=factbj*bij1*bji2
        dbij=factbi*(dbij1*bij2 + bij1*dbij2)
        dbji=factbj*(dbij1*bji2 + bij1*dbji2)
        qforcetot=(qforce(i)*pij*dbij + qforcej*dbji)
        fpx(i)=fpx(i)+qforcetot*xij*r1
        fpy(i)=fpy(i)+qforcetot*yij*r1
        fpz(i)=fpz(i)+qforcetot*zij*r1
     end do
     ncount=ncount + (map+1)
  end do
end subroutine force_sasa
