c========================================================================
      subroutine do_dynamic_input (xp0,yp0,zp0,nato_slt,nmol_slv
     &     ,unitc,unite,avogad,kdynamic)
c========================================================================
c  Write Force field for solute (input to EE and dynamic programs ) 
c  Written by P. Procacci UNIFI 1998   
c========================================================================

*======================= DECLARATIONS ==================================

      use parst
      use cpropar

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      integer nato_slt,kdynamic,nmol_slv
      real*8  unitc,unite,avogad,xp0(*),yp0(*),zp0(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      integer la,lb,lc,naux,m,l11,l22
      real*8 xr1,yr1,zr1,xr2,yr2,zr2,xr3,yr3,zr3,x21,y21,z21,x32,y32
     &     ,z32,rs12,rs21,rs32,dcc2,cb,bb,gcpu,pi,aux,sign,phase

      gcpu = unite*avogad/1000.d0/4.184d0
      pi = dacos(-1.d0) 
      
c============write bonds table ===================================

      l11 = lbond-lbond_slv*nmol_slv
      write(kdynamic,*) l11, ' ===== BOND TABLE  ' 
      do m=1,l11
c--     compute actual distance 
        la=lbnd(1,m)
        lb=lbnd(2,m)
        xr1=xp0(la)
        yr1=yp0(la)
        zr1=zp0(la)
        xr2=xp0(lb)
        yr2=yp0(lb)
        zr2=zp0(lb)
        x21=xr2-xr1
        y21=yr2-yr1
        z21=zr2-zr1
        rs21=DSQRT(x21**2+y21**2+z21**2)
        WRITE(kdynamic
     &       ,'(1h ,a4,2x,a4,5x,i6,2x,i6,2x,f7.2,f10.3,10x,i4,f10.4)'
     &       )betb(lbnd(1,m)),betb(lbnd(2,m)),lbnd(1,m),lbnd(2,m)
     &       ,potbo(m,1)*gcpu,potbo(m,2),m,rs21
      end do 
      
c============write bends table ===================================

      l11 = lbend-lbend_slv*nmol_slv
      write(kdynamic,*) l11, ' ===== BEND TABLE  ' 
      do m=1,l11
        la=lbndg(1,m)
        lb=lbndg(2,m)
        lc=lbndg(3,m)
        xr1=xp0(la)
        yr1=yp0(la)
        zr1=zp0(la)
        xr2=xp0(lb)
        yr2=yp0(lb)
        zr2=zp0(lb)
        xr3=xp0(lc)
        yr3=yp0(lc)
        zr3=zp0(lc)
        x21=xr1-xr2
        y21=yr1-yr2
        z21=zr1-zr2
        x32=xr3-xr2
        y32=yr3-yr2
        z32=zr3-zr2
        rs12=x21**2+y21**2+z21**2
        rs32=x32**2+y32**2+z32**2
        dcc2=DSQRT(rs12*rs32)
        cb=(x21*x32+y21*y32+z21*z32)/dcc2
        bb=DACOS(cb)
        WRITE(kdynamic,1077) betb(lbndg(1,m)),betb(lbndg(2,m))
     &       ,betb(lbndg(3,m)),lbndg(1,m),lbndg(2,m),lbndg(3,m)
     &       ,potbe(m,1)*gcpu,potbe(m,2)*180.d0/pi,m,bb*180.d0/pi
1077    FORMAT(1h ,a4,2x,a4,2x,a4,5x,i6,2x,i6,2x,i6,2x,f7.2,f10.4,10x
     &       ,i4,f10.3)  
      end do   

      l11 = ltors-ltors_slv*nmol_slv
      l22=0  ! do this to account for moved i-torsions
      do m=1,litor
        if(ss_index(litr(1,m)).EQ.1) l22=l22+1
      end do
      write(kdynamic,*) l11+l22, l22, ' ===== TORSION TABLE  ' 
      
c============write proper torsion table ==========================

      DO m=1,l11
!       V_tors = sign(en)*abs(V0) ( 1+cos( abs(en)*phi-phase))
!       where phase=180 if V0 >0 and phase=0.0 if V0< 0         
        sign=potto(m,2)/abs(potto(m,2)) ! handle negative V0
        if(potto(m,1).gt.0.d0) THEN ! sign of v0 is phase 180/0
          phase=180.0 
        ELSE
          phase=0.0 
        END IF
        WRITE(kdynamic,'(1h ,4(a4,2x),2x,4(i4,1x),f9.4,i5,5x,f9.3,i4)'
     &       )betb(ltor(1,m)),betb(ltor(2,m)),betb(ltor(3,m))
     &       ,betb(ltor(4,m)),ltor(1,m),ltor(2,m),ltor(3,m),ltor(4,m)
     &       ,sign*abs(potto(m,1))*gcpu,int(abs(potto(m,2))),phase,m
      END DO
      
c============write improper torsion table ========================

      DO  m=1,litor
        if(ss_index(litr(1,m)).EQ.1) THEN
          if(potit(m,3).gt.0.d0) THEN 
10097       FORMAT( "#",a4,2x,3(a4,2x),2x,4(i4,1x),f9.4,f8.2,5x,i4) 
            WRITE(kdynamic,10097) betb(litr(1,m)),betb(litr(2,m))
     &           ,betb(litr(3,m)),betb(litr(4,m)),litr(1,m),litr(2,m)
     &           ,litr(3,m),litr(4,m),potit(m,1)*gcpu,potit(m,2)*180.d0
     &           /pi,m
          ELSE
            naux=DINT(potit(m,2)+0.5d0)
            WRITE(kdynamic,'(1h ,4(a4,2x),2x,4(i4,1x),f9.4,i5,5x,i4)'
     &           )betb(litr(1,m)),betb(litr(2,m)),betb(litr(3,m))
     &           ,betb(litr(4,m)),litr(1,m),litr(2,m),litr(3,m),litr(4,m
     &           ),potit(m,1)*gcpu,naux,m
          END IF
        end if
      END DO



c============write atom table  ===================================

      write(kdynamic,*) nato_slt, ' ===== ATOM TABLE  ' 
      do m=1,nato_slt
        write(kdynamic,3009) mass(m),xp0(m),yp0(m),zp0(m),chrge(m)
     &       *dsqrt(unitc),pnbd2(nbtype(m)),pnbd1(nbtype(m)),m
3009    FORMAT(f9.4,3f10.4,f10.5,2f10.4,i5)  
      end do
      l11 = int14p-int14p_slv*nmol_slv
      write(kdynamic,*) l11, ' ===== 14 CONTACT TABLE  ' 
      do m=1,l11
        write(kdynamic,3010) int14(1,m),int14(2,m),m
3010    format(2i6,i10)
      end do

*===============END OF EXECUTABLE STATEMENTS ============================

      RETURN
      END 
