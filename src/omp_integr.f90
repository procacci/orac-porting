module omp_integr

  use parst, only:m1
  implicit none
  integer mcor1,mcor2,mcor3
  parameter(mcor1=16,mcor2=16,mcor3=16)
  REAL*8  time_n0,time_n1,time_l,time_m,time_h,time_u,time_f,time_cn,time_in
  REAL*8  tcpu_n0,tcpu_n1,tcpu_l,tcpu_m,tcpu_h,tcpu_u,tcpu_f,tcpu_cn,tcpu_in
  integer nthr,nthr1,nthr2,m8t,mpp8,m1t,cache_line_size,cache_line_size1,cache_line_size2
  logical*1 array_omp(4*m1,mcor3),array1_omp(4*m1,mcor1),array2_omp(4*m1,mcor2)
  integer level_integr,level_cnst,level_n0,level_n1
  integer level_shell_m,level_shell_l,level_shell_h,level_update,level_fftw,iprot(m1),typei(m1)
  LOGICAL omp_timing,omp_dynamic

contains  

  !=======================================================================
  subroutine copyomp(n,x1,y1,z1,x2,y2,z2)
    !=======================================================================
    ! OMP copy version 
    !=======================================================================
    real*8 x1(n),y1(n),z1(n),x2(n),y2(n),z2(n)
    integer n,i

    !$OMP PARALLEL DO 
    do i=1,n
       x2(i)=x1(i)
       y2(i)=y1(i)
       z2(i)=z1(i)
    end do
    !$OMP END PARALLEL DO 
    return 
  end subroutine copyomp

  !=======================================================================
  subroutine zeroomp(x1,y1,z1,n)
    !=======================================================================
    ! OMP zero version 
    !=======================================================================
    implicit NONE
    real*8 x1(n),y1(n),z1(n)
    INTEGER n,i

    !$OMP PARALLEL DO 
    do i=1,n
       x1(i)=0.d0
       y1(i)=0.d0
       z1(i)=0.d0
    end do
    !$OMP END PARALLEL DO 
    return 
  end subroutine zeroomp

  !=======================================================================
  subroutine  coord_inner_omp(ntap,nprot,iprot,co,oc,lx,ly,lz,xp0,yp0,zp0,xpa,ypa,zpa,xpcm,ypcm,zpcm, & 
       xpcma,ypcma,zpcma) 
    !=======================================================================
    !    OMP recompute coord from lx and cm coords and change frame
    !=======================================================================
    implicit none 
    ! args
    real*8 lx(ntap),ly(ntap),lz(ntap) ! input: atomic ORT mol frame
    real*8 xpcma(nprot),ypcma(nprot),zpcma(nprot) !input: com cryst frame
    integer iprot(ntap)! input: vector labeling the i mol for j atom index
    integer ntap,nprot

    real*8 xp0(ntap),yp0(ntap),zp0(ntap) ! output: atomic ORT fixed frame
    real*8 xpcm(nprot),ypcm(nprot),zpcm(nprot) !output com ORT frame
    real*8 xpa(ntap),ypa(ntap),zpa(ntap) ! output: atomic  com cryst frame
    real*8 co(3,3),oc(3,3)
    !local 
    integer i,j
    real*8 xc(nprot),yc(nprot),zc(nprot)
    REAL*8  tela1,tela2,tcpu1,tcpu2
    REAL*8  treal

    if(omp_timing)  call timer(treal,tcpu1,tela1)
    !$OMP PARALLEL DO
    DO i=1,nprot
       !       change origin generate new xp0 coord in Orthog  syst.
       xc(i)=co(1,1)*xpcma(i)+co(1,2)*ypcma(i)+co(1,3)*zpcma(i)
       yc(i)=co(2,1)*xpcma(i)+co(2,2)*ypcma(i)+co(2,3)*zpcma(i)
       zc(i)=co(3,1)*xpcma(i)+co(3,2)*ypcma(i)+co(3,3)*zpcma(i)
       !       stor ORT xpcm coordinate
       xpcm(i)=xc(i)
       ypcm(i)=yc(i)
       zpcm(i)=zc(i)
    END DO
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO
    do i=1,ntap
       xp0(i)=lx(i)+xc(iprot(i))
       yp0(i)=ly(i)+yc(iprot(i))
       zp0(i)=lz(i)+zc(iprot(i))
       xpa(i)=oc(1,1)*xp0(i)+oc(1,2)*yp0(i)+oc(1,3)*zp0(i)
       ypa(i)=oc(2,1)*xp0(i)+oc(2,2)*yp0(i)+oc(2,3)*zp0(i)
       zpa(i)=oc(3,1)*xp0(i)+oc(3,2)*yp0(i)+oc(3,3)*zp0(i)
    END DO
    !$OMP END PARALLEL DO        
    if(omp_timing) THEN 
       call timer(treal,tcpu2,tela2)
       time_in=time_in+tela2-tela1
       tcpu_in=tcpu_in+tcpu2-tcpu1
    END IF
    return
  end subroutine coord_inner_omp

  !=======================================================================
  SUBROUTINE comp_fcm_omp(nprot,ntap,iprot,fpx,fpy,fpz,fcx,fcy,fcz, &     
       mass,tmass1,oc)
    !=======================================================================
    !     OMP version
    !======================== DECLARATIONS ================================*

    IMPLICIT none

    REAL*8  tela1,tela2,tcpu1,tcpu2
    REAL*8  treal

    !----------------------------- ARGUMENTS ------------------------------*

    INTEGER nprot,ntap,iprot(*)
    REAL*8  fpx(ntap),fpy(ntap),fpz(ntap),fcx(nprot),fcy(nprot),fcz(nprot),mass(nprot), & 
         tmass1(*),oc(3,3)

    !------------------------- LOCAL VARIABLES ----------------------------*

    INTEGER i
    REAL*8  xc(nprot),yc(nprot),zc(nprot)

    !----------------------- EXECUTABLE STATEMENTS ------------------------*

    if(omp_timing)  call timer(treal,tcpu1,tela1)

    call zeroomp(xc,yc,zc,nprot) 

    !     computes total (COM) forces on each of the NPROT "molecule"  
    !     This do OMP loop does not work with pgforgtran;  
    !!$OMP PARALLEL DO REDUCTION(+:xc,yc,zc) SHARED(iprot)
    do i=1,ntap 
       xc(iprot(i))=xc(iprot(i))+fpx(i)
       yc(iprot(i))=yc(iprot(i))+fpy(i)
       zc(iprot(i))=zc(iprot(i))+fpz(i)
    end do
    !!$OMP END PARALLEL DO       
    !------------------------------------------------------
    !     Reduces atomic foorce 
    !$OMP PARALLEL DO REDUCTION(+:fpx,fpy,fpz)
    do i=1,ntap 
       fpx(i)=fpx(i)-mass(i)*xc(iprot(i))*tmass1(iprot(i))
       fpy(i)=fpy(i)-mass(i)*yc(iprot(i))*tmass1(iprot(i))
       fpz(i)=fpz(i)-mass(i)*zc(iprot(i))*tmass1(iprot(i))
    end do
    !$OMP END PARALLEL DO       
    !------------------------------------------------------

    !     computes crystal coord COM forces
    !$OMP PARALLEL DO
    do i=1,nprot
       fcx(i)=oc(1,1)*xc(i)+oc(1,2)*yc(i)+oc(1,3)*zc(i)
       fcy(i)=oc(2,1)*xc(i)+oc(2,2)*yc(i)+oc(2,3)*zc(i)
       fcz(i)=oc(3,1)*xc(i)+oc(3,2)*yc(I)+oc(3,3)*zc(I)
    END DO
    !$OMP END PARALLEL DO       

    if(omp_timing) THEN 
       call timer(treal,tcpu2,tela2)
       time_in=time_in+tela2-tela1
       tcpu_in=tcpu_in+tcpu2-tcpu1
    END IF
    RETURN
  END SUBROUTINE comp_fcm_omp

  !=======================================================================
  subroutine assign_prot(nprot,protl,ss_index,iprot,typei,tmass,tmass1)
    !=======================================================================
    !    NON-OMP subroutine to compute ntap-vector label belonging mol 
    !=======================================================================
    implicit NONE 

    integer protl(*),ss_index(*)  ! input: molecule/type map for atomic indices
    integer nprot    ! input: number of molecules
    integer iprot(*),typei(*)! output map vectors
    real*8 tmass(*),tmass1(*)

    integer i,m ,count,i1,j,offset   !local 
    count=0
    do i=1,nprot
       offset=protl(count+1)
       typei(i)=ss_index(protl(count+1+1)) ! assign type to i-th mol
       count=count+offset+1
    end do
    count=0
    do i=1,nprot 
       m=protl(1+count)
       do i1=1,m
          j=protl(count+1+i1)
          iprot(j)=i  ! assign atom j to i-th mol
       end do
       count=count+m+1
       tmass1(i)=1.d0/tmass(i) ! needed in frequently computed do loop in comp_fcm_omp
    end do
    return
  end subroutine assign_prot

  !=======================================================================
  subroutine cache_line(n,nthr,lcl,cl_size)
    !=======================================================================
    !    NON-OMP subroutine to compute task splitting vector[s] according to
    !    cache-line size to avoid false sharing and cache misses
    !=======================================================================
    implicit NONE
    integer n,nthr,cl_size,nchunks,i,j,k,l,nbalance(nthr),m
    logical*1 lcl(n,nthr)
    nchunks=1+n/cl_size
    lcl=.false.
    do j=1,nthr
       nbalance(j)=0
       i=0
       do k=j,nchunks,nthr
          if(i.ge.n) exit
          do l=1,cl_size
             i=l+(k-1)*cl_size
             if(i.le.n) THEN 
                lcl(i,j)=.true.
             END IF
             nbalance(j)=nbalance(j)+1
          end do
       end do
    end do
    return
  end subroutine cache_line


  !------------------------------------------------------------------------
  subroutine split_map(nato,mapnl,grppt,ngrp,kprint,mapnlt,itask,istride) 
    !------------------------------------------------------------------------
    !  Only for PARALLEL execution.
    !  Receive the full 1-2 1-3 1-4 contact map and build the iproc
    !  contact map based on particle decomposition. N.B. On return mapnl is overwritten
    !  with the new contact map  of the iproc thread. 

    implicit none
    ! ARGUMENTS 
    INTEGER, intent(in)  ::  nato,grppt(2,*),ngrp,kprint,itask,istride
    INTEGER, intent(in)  ::  mapnl(*) 
    INTEGER, intent(out) :: mapnlt(*)     
    ! LOCAL VARIABLES 
    INTEGER :: itaski,na,nb,la,noff,tmp(nato),i,i1,j


    !     split 1-2 1-3 1-4  contact map 
    na=0
    nb=0
100 format(6x,'Rebuilding 1-2 1-3 1-4 map        Proc ---->',i4)
    do i=1,ngrp
       DO i1=grppt(1,i),grppt(2,i)
          la=mapnl(1+na)
          !       while unrolling the full map, built the map for the iproc thread 
          DO j=1,la
             tmp(j)=mapnl(j+1+na)
          END DO
          if(array_omp(i,itask)) THEN
             mapnlt(nb+1)=la
             do j=1,la
                mapnlt(nb+1+j)=tmp(j)
             end do
             noff=mapnlt(nb+1)+1
             nb=nb+noff
          endif
          noff=mapnl(1+na)+1
          na=na+noff
       end do
    end do
    if(itask.eq.1) write(kprint,101) na
101 format(16x, ' Old MAPNL dimension =',i10)

    WRITE(kprint,200) itask-1,nb 
200 format(16x, ' New map done for IPROC =',i4, / &  
         , 16x, ' New MAPNL dimension =',i10)  
  end subroutine split_map


  !---------------------------------------------------------------------------------------    
  subroutine split_cnstr(nprot,cnst_protl,cnst_protlt,cnstp_threads,itask,nthr)
    !---------------------------------------------------------------------------------------    
    !  Only for PARALLEL execution.  This routine enters with a core identifier ITASK. 
    !  Receive the full constraint map and build the ITASK contact map based on particle 
    !  decomposition. Old full constraint map is flushed in mtsmd on return.
    !  ITASK maps are then used with same PARALLEL construct in rattle_correc and rattel_verlet.
    !  NB: OMP or MPI  execution occurs only when STRETCHING HEAVY is used.
    !-----------------------------------written by P. Procacci 2015-------------------------    
    use unit, only:kprint
    implicit none

    !   arguments----------------
    integer nprot,cnst_protl(*),cnst_protlt(*),itask,cnstp_threads(nthr),nthr

    !   local variables 
    integer i,j,count,countt,tmp(nprot),itaski,cnstp
    count=0
    countt=0
    if(itask.eq.1) write(kprint,'(20x,"Constraint load balance",/)')
    Do i=1,nprot
       cnstp=cnst_protl(1+count)
       DO j=1,cnstp                   !while unrolling the full map, built the map for the iproc thread 
          tmp(j)=cnst_protl(j+1+count)
       END DO
       if(array2_omp(i,itask)) THEN
          cnst_protlt(countt+1)=cnstp
          do j=1,cnstp
             cnst_protlt(countt+1+j)=tmp(j)
          end do
          countt=countt+cnstp+1       !total count in the map is cnstp+1 per "protein"
       endif
       count=count+cnstp+1
    end do
    cnstp_threads(itask)=countt
    write(kprint,100) itask,countt
100 FORMAT(13x," OMP_Task=", i5,"   CNST_MAP_DIM=",i10)
    return
  end subroutine split_cnstr
 

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  subroutine lc_index_omp(indmax,ncx,ncy,ncz,nind,indxi,indxj,indxk,i_ind,j_ind, &
       k_ind,n_ind,ctoff,co,npoints)

    use unit, only:kprint
    implicit none
    !----------------------------
    ! 
    !----------------------------
    real*8 dist_ijk
    real*8 ctoff,co(3,3)

    integer ncx,ncy,ncz,indmax
    integer nind(2)
    integer  indxi(2,*),indxj(2,*),indxk(2,*)
    integer i_ind(*),j_ind(*),k_ind(*),n_ind(*)

    real*8 dx,dy,dz,sqcut,rmin

    integer i,j,k,n,il,npoints
    integer imax,jmax,kmax
    integer nxmax,nymax,nzmax
    integer warnx,warny,warnz
    REAL*8 map_rmin(npoints)

    !----------------------------
    ! initialisations
    !----------------------------
    sqcut=ctoff**2
    dx=2.d0/ncx
    dy=2.d0/ncy
    dz=2.d0/ncz
    imax=0
    jmax=0
    kmax=0

!$OMP PARALLEL DO PRIVATE(I,J,K) 
    do il=1,npoints
       i=i_ind(il)
       j=j_ind(il)
       k=k_ind(il)
       n=n_ind(il)
       map_rmin(il)=dist_ijk(i,j,k,dx,dy,dz,co)
    end do
!$OMP END PARALLEL DO 


    !----------------------------
    ! calcul de la demi-sphere
    ! des indices
    !----------------------------
    nind(1)=1
    indxi(1,1)=0
    indxj(1,1)=0
    indxk(1,1)=0

    nind(2)=1
    indxi(2,1)=0
    indxj(2,1)=0
    indxk(2,1)=0

    do il=1,npoints
       i=i_ind(il)
       j=j_ind(il)
       k=k_ind(il)
       n=n_ind(il)
       if(map_rmin(il).lt.sqcut) then
          nind(n)=nind(n)+1
          if(nind(n).gt.INDMAX) then
             write(kprint,*) "=== INDMAX is too small",indmax 
             write(kprint,*) "=== Decrease cut-off or increase ", &
                             "_INDMAX_ and _CELLMAX_ in config.H"
             stop
          endif
          indxi(n,nind(n))=i
          indxj(n,nind(n))=j
          indxk(n,nind(n))=k
          if(imax.lt.abs(i)) imax=abs(i)
          if(jmax.lt.abs(j)) jmax=abs(j)
          if(kmax.lt.abs(k)) kmax=abs(k)
          if(i.eq.0 .and. j.eq.0 .and. k.eq.0) nind(n)=nind(n)-1
       endif
    end do
    !----------------------------------------------------------------------------------------
    ! pour eviter que, lors du comptage des paires d'atomes une cellule n'apparaisse
    ! deux fois, il faut ajouter le test suivant : si ncx est pair, il faut imax <= ncx/2
    ! si ncx est impair, il faut imax <= (ncx+1)/2 (idem pour jmax et kmax)
    !----------------------------------------------------------------------------------------

    nxmax=(ncx+1)/2
    nymax=(ncy+1)/2
    nzmax=(ncz+1)/2
    warnx=0
    warny=0
    warnz=0
    if(imax.ge.nxmax) warnx=1
    if(jmax.ge.nymax) warny=1
    if(kmax.ge.nzmax) warnz=1
    if(warnx.eq.1 .or. warny.eq.1 .or. warnz.eq.1) then
       write(kprint,*) "des cellules risquent d'etre comptees deux fois"
       write(kprint,*)  "diminuez le cutoff ou :"
       if(warnx.eq.1) write(kprint,*) "augmentez ncx"
       if(warny.eq.1) write(kprint,*) "augmentez ncy"
       if(warnz.eq.1) write(kprint,*) "augmentez ncz"
       stop
    endif
    return
  end subroutine lc_index_omp

  !--------------------------------------------------------------------------
  subroutine mapp_lc(ncx,ncy,ncz,i_ind,j_ind,k_ind,n_ind,npoints)
    !--------------------------------------------------------------------------
    !    called once by mtsmd. Set up map for dist_ijk parallel comp. 
    !--------------------------------------------------------------------------

    !    args
    integer ncx,ncy,ncz,i_ind(*),j_ind(*),k_ind(*),n_ind(*),npoints
    ! local 
    integer i,j,k,il,n,istart,jstart,kstart

    il=0
    do  n=1,2
       istart=1-ncx
       if(n.eq.1) istart=0
       do i=istart,ncx-1
          jstart=1-ncy
          if(n.eq.1 .and. i.eq.0) jstart=0
          do j=jstart,ncy-1
             kstart=1-ncz
             if(n.eq.1 .and. i.eq.0 .and. j.eq.0) kstart=0
             do  k=kstart,ncz-1
                il= il+ 1
                i_ind(il)=i 
                j_ind(il)=j 
                k_ind(il)=k 
                n_ind(il)=n
             end do
          end do
       end do
    end do
    npoints=il
    return
  end subroutine mapp_lc
  SUBROUTINE comp_stress_add_tpg(nprot,virial,co,oc,xcm,ycm,zcm,fcax,fcay,fcaz,vco,masspp,tn0)
    
    !======================== DECLARATIONS ================================*
    
    IMPLICIT none
    
    !----------------------------- ARGUMENTS ------------------------------*
    
    REAL*8  xcm(*),ycm(*),zcm(*),fcax(*),fcay(*),fcaz(*)
    REAL*8  co(3,3),oc(3,3),vco(3,3),masspp(*),virial(3,3),st(3,3),tn0
    integer nprot
    
    !------------------------- LOCAL VARIABLES ----------------------------*
    
    INTEGER i
    REAL*8  xc,yc,zc
    
    !----------------------- EXECUTABLE STATEMENTS ------------------------*


    virial = 0.d0 
    st = 0.d0 
!$OMP  PARALLEL DO SCHEDULE(STATIC) REDUCTION(+:virial) 
    DO i=1,nprot
       xc=co(1,1)*fcax(i)+co(1,2)*fcay(i)+co(1,3)*fcaz(i)
       yc=co(2,1)*fcax(i)+co(2,2)*fcay(i)+co(2,3)*fcaz(i)
       zc=co(3,1)*fcax(i)+co(3,2)*fcay(i)+co(3,3)*fcaz(i)
       virial(1,1)=virial(1,1)+xc*xcm(i)
       virial(1,2)=virial(1,2)+yc*xcm(i)
       virial(1,3)=virial(1,3)+zc*xcm(i)
       virial(2,1)=virial(2,1)+xc*ycm(i)
       virial(2,2)=virial(2,2)+yc*ycm(i)
       virial(2,3)=virial(2,3)+zc*ycm(i)
       virial(3,1)=virial(3,1)+xc*zcm(i)
       virial(3,2)=virial(3,2)+yc*zcm(i)
       virial(3,3)=virial(3,3)+zc*zcm(i)
    END DO
!$OMP END PARALLEL DO 
    CALL DGEMM('N','T',3,3,3,1.0D0,virial,3,oc,3,0.0D0,st,3)
    CALL correc3x3(vco,st,masspp,tn0)
    return 
  end SUBROUTINE comp_stress_add_tpg
end module omp_integr

      



