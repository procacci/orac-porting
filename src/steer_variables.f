*======================================================================*
      subroutine steer_variables(nbonds_added,nbends_added
     &     ,nitors_added,atom_b1,atom_b2,atom_be1,atom_be2,atom_be3
     &     ,atom_it1,atom_it2,atom_it3,atom_it4,xp0,yp0,zp0
     &     ,bond_added,bend_added,tors_added,eqdied,n,m)
*======================================================================*
*  This subroutine computes the bond bend and tors coordinates
*  for steered dynamics. These coordinates are used by mts_test routine 
*  for computing the work along those coordinates to be used 
*  in the Jarzinky identity. 
*======================================================================*

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER nbonds_added,nbends_added ,nitors_added,atom_b1(*)
     &     ,atom_b2(*),atom_be1(*),atom_be2(*),atom_be3 (*),atom_it1(*)
     &     ,atom_it2(*),atom_it3(*),atom_it4(*)

      REAL*8 xp0(*),yp0(*),zp0(*),bond_added(*),bend_added(*)
     &     ,tors_added(*),eqdied(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,n,m
      REAL*8  x12,y12,z12,x23,y23,z23,x21,y21,z21,x32,y32,z32,x43,y43
     &     ,z43,rs12,rs32,cb,rsq21,rsq32,rsq43,rsp21,rsp32,rsp43,cb1,cb2
     &     ,cb3,sb1,sb2,sb3,aux,coa,xr1,xr2,xr3,xr4,yr1,yr2,yr3,yr4,zr1
     &     ,zr2,zr3,zr4,XN,YN,ZN
      

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      

c--   compute bonds
      DO i=1,nbonds_added
        x12 = xp0(atom_b1(i))- xp0(atom_b2(i))
        y12 = yp0(atom_b1(i))- yp0(atom_b2(i))
        z12 = zp0(atom_b1(i))- zp0(atom_b2(i))
        bond_added(i) = dsqrt(x12**2+y12**2+z12**2)
      END DO
c--   compute bends
      DO i=1,nbends_added
        x12=xp0(atom_be1(i))- xp0(atom_be2(i))
        y12=yp0(atom_be1(i))- yp0(atom_be2(i))
        z12=zp0(atom_be1(i))- zp0(atom_be2(i))
        x23=xp0(atom_be3(i))- xp0(atom_be2(i))
        y23=yp0(atom_be3(i))- yp0(atom_be2(i))
        z23=zp0(atom_be3(i))- zp0(atom_be2(i))
        rs12=x12**2+y12**2+z12**2
        rs32=x23**2+y23**2+z23**2
        cb=(x12*x23+y12*y23+z12*z23)/DSQRT(rs12*rs32)
        bend_added(i)=DACOS(cb)
      END DO
c--   compute tors
      DO i=1,nitors_added
        xr1=xp0(atom_it1(i)) 
        xr2=xp0(atom_it2(i)) 
        xr3=xp0(atom_it3(i)) 
        xr4=xp0(atom_it4(i)) 
        yr1=yp0(atom_it1(i)) 
        yr2=yp0(atom_it2(i)) 
        yr3=yp0(atom_it3(i)) 
        yr4=yp0(atom_it4(i)) 
        zr1=zp0(atom_it1(i)) 
        zr2=zp0(atom_it2(i)) 
        zr3=zp0(atom_it3(i)) 
        zr4=zp0(atom_it4(i)) 
        x21=xr2-xr1
        y21=yr2-yr1
        z21=zr2-zr1
        x32=xr3-xr2
        y32=yr3-yr2
        z32=zr3-zr2
        x43=xr4-xr3
        y43=yr4-yr3
        z43=zr4-zr3
        XN = (y21*z32-y32*z21)
        YN = (z21*x32-z32*x21)
        ZN = (x21*y32-x32*y21) 
        rsq21=x21**2+y21**2+z21**2
        rsq32=x32**2+y32**2+z32**2
        rsq43=x43**2+y43**2+z43**2
        rsp21=DSQRT(rsq21)
        rsp32=DSQRT(rsq32)
        rsp43=DSQRT(rsq43)
        cb1=(x21*x32+y21*y32+z21*z32)/(rsp21*rsp32)
        cb2=(x43*x32+y43*y32+z43*z32)/(rsp43*rsp32)
        cb3=(x21*x43+y21*y43+z21*z43)/(rsp21*rsp43)
        sb1=DSQRT(1.0d0-cb1**2)
        sb2=DSQRT(1.0d0-cb2**2)
        sb3=DSQRT(1.0d0-cb3**2)
        aux=sb1*sb2
        coa=(cb1*cb2-cb3)/aux
        tors_added(i)=DACOS(coa)
        aux=sign(tors_added(i),XN*x43+YN*Y43+ZN*z43)
        if(aux.lt.0.d0.and.dabs(aux).gt.1.57079)  THEN
          tors_added(i)=2*dacos(-1.d0)+aux
          if(abs(eqdied(i+n-m)-tors_added(i)).gt.3.14159D0)
     &         tors_added(i)=-2*3.14159D0+tors_added(i)
c          write(6,601) 180*tors_added(i)/dacos(-1.d0),tors_added(i),XN
c     &         *x43+YN*Y43+ZN*z43,aux,eqdied(i+n-m)
        ELSE IF(aux.lt.0.d0.and.dabs(aux).le.1.57079)  THEN
          tors_added(i)= aux
          if(abs(eqdied(i+n-m)-tors_added(i)).gt.3.14159D0)
     &         tors_added(i)=2*3.14159D0+tors_added(i)
c          write(6,602) 180*tors_added(i)/dacos(-1.d0),tors_added(i),XN
c     &         *x43+YN*Y43+ZN*z43,aux,eqdied(i+n-m)

        ELSE
          tors_added(i)=aux
          if(abs(eqdied(i+n-m)-tors_added(i)).gt.3.14159D0) 
     &         tors_added(i)=2*3.14159D0+tors_added(i)
          if(abs(eqdied(i+n-m)-tors_added(i)).gt.3.14159D0) 
     &         tors_added(i)=-4*3.14159D0+tors_added(i)
c          write(6,603) 180*tors_added(i)/dacos(-1.D0),tors_added(i),XN
c     &         *x43+YN*Y43+ZN*z43,aux,eqdied(i+n-m)
c601       FORMAT(" 1 steerv -----",5G12.4," 1.57079") 
c602       FORMAT(" 2 steerv -----",5G12.4," 1.57079") 
c603       FORMAT(" 3 steerv -----",5G12.4," 1.57079") 
        END IF
      END DO
      RETURN 
      END 
