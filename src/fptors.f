      SUBROUTINE fptors(ss_index,ltor,ltors,xp0,yp0,zp0,pota
     &     ,potb,utors_slt,utors_slv,fpx,fpy,fpz,mb,itask)

************************************************************************
*                                                                      *
*     Compute the contribution from proper torsions to                 *
*     macromolecular forces and energy.                                *
*                                                                      *
*     CO      :  Transform the positions to orthogo alized        (I)  *
*                crystallographic frame.                               *
*                >> real*8 CO(3,3) <<                                  *
*     LTOR    :  List of torsions for the macromolecule.          (I)  *
*                >> integer  LTOR(4,N1) <<                             *
*     LTORS   :  Number of torsions.                              (I)  *
*     N1      :  Physical column dimension of LTOR.               (I)  *
*     XP0     :  Coordinates of the macromolecule.                (I)  *
*     YP0        >> real*8 XP0(NATO), YP0(NATO), ZP0(NATO) <<          *
*     ZP0                                                              *
*                                                                      *
*     NATO    :  Number of atoms which the macromolecule is       (I)  *
*                composed of.                                          *
*     POTA    :  List of potential parameters : multiplicity      (I)  *
*                of the torsion.                                       *
*                >> real*8 POTA(N1) <<                                 *
*     POTB    :  List of potential parameters : torsion           (I)  *
*                force constant.                                       *
*                >> real*8 POTB(N1) <<                                 *
*     UTORS   :  Potential energy contribution from torsions.     (I)  *
*     FPX     :  List of forces acting on the macromolecule       (I)  *
*     FPY        atoms.                                                *
*     FPZ        >> real*8 FPX(NATO), FPY(NATO), FPZ(NATO) <<          *
*                                                                      *
*---- Last update 04/24/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*     ReWritten by Piero Procacci Unifi 2013                           *
*                                                                      *
*     EXTERNALS NONE                                                   *
*                                                                      *
************************************************************************
* N.B. The phase is not passed. ORAC uses the following conventions 
* to handle the phase: 
*    if Gamma=0 the sign of V0 is negative  
*    if Gamma=180 the sign of V0 is positive 
*    The actual sign of V0 specified in input is assimilated
*    in the multiplicity sign. So, if V0 is given a negative sign 
*    then n is negative. 
*    Phase of 60 and 90 are handled through the multiplicity
*    see rdaprm for details.
*======================= DECLARATIONS ==================================

#ifdef _OMP_
      use omp_integr, only:nthr2,array2_omp
#endif
      use cpropar, only:docking,minimize
      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER ltors,ltor(4,*),ss_index(*),mb,itask
      REAL*8  utors_slt,utors_slv,pota(*),potb(*),xp0(*),yp0(*),
     x        zp0(*),fpx(mb),fpy(mb),fpz(mb)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,l1,l2,l3,l4,ntph,tmp,type,nscrew
      REAL*8  xr1,xr2,xr3,xr4,yr1,yr2,yr3,yr4,zr1,zr2,zr3,zr4,x21,x32
     &     ,x43,y21,y32,y43,z21,z32,z43,rsq21,rsq32,rsq43,rsp21,rsp32
     &     ,rsp43,dummy,phase
      REAL*8  cb1,cb2,cb3,sb1,sb2,sb3,aux,aux1,aux2,auxa,auxb,auxc,
     x        dbx11,dbx31,dbx12,dbx22,dbx32,dbx13,dbx23,dbx33,dbx24,
     x        dbx34,dby11,dby31,dby12,dby22,dby32,dby13,dby23,dby33,
     x        dby24,dby34,dbz11,dbz31,dbz12,dbz22,dbz32,dbz13,dbz23,
     x        dbz33,dbz24,dbz34,uux1,uux2,uux3,uux4,uuy1,uuy2,uuy3,
     x     uuy4,uuz1,uuz2,uuz3,uuz4,aux3,auxn,dauxn,xn,yn,zn,aux4,bb
      REAL*8  coa,soa,dsoa,enrg,qforce,sig,utors(2),ranf,sin60
      LOGICAL screw_torsion,lwrite
      EXTERNAL ranf
      DATA sin60  /0.866025/ 


*==================== EXECUTABLE STATEMENTS ============================
      
      DO i=1,2
         utors(i)=0.0D0
      END DO

      screw_torsion=.false.
      lwrite=(.not.docking).AND.(.not.minimize) 
      nscrew=0
      dummy=3.0 
      DO i=1,ltors
#ifdef _OMP_
        if(array2_omp(i,itask)) THEN 
#endif
          if(abs(potb(i)).lt.1.d-12) go to 10
          l1=ltor(1,i)
          l2=ltor(2,i)
          l3=ltor(3,i)
          l4=ltor(4,i)

*-----------------------------------------------------------------------
*------------ No bonded interaction exists between solvent and solute --
*-----------------------------------------------------------------------
              
          type=ss_index(l1)
          xr1=xp0(l1)
          yr1=yp0(l1)
          zr1=zp0(l1)
          xr2=xp0(l2)
          yr2=yp0(l2)
          zr2=zp0(l2)
          xr3=xp0(l3)
          yr3=yp0(l3)
          zr3=zp0(l3)
          xr4=xp0(l4)
          yr4=yp0(l4)
          zr4=zp0(l4)
1009      x21=xr2-xr1
          y21=yr2-yr1
          z21=zr2-zr1
          x32=xr3-xr2
          y32=yr3-yr2
          z32=zr3-zr2
          x43=xr4-xr3
          y43=yr4-yr3
          z43=zr4-zr3
          rsq21=x21**2+y21**2+z21**2
          rsq32=x32**2+y32**2+z32**2
          rsq43=x43**2+y43**2+z43**2
          rsp21=DSQRT(rsq21)
          rsp32=DSQRT(rsq32)
          rsp43=DSQRT(rsq43)
          cb1=(x21*x32+y21*y32+z21*z32)/(rsp21*rsp32)
          cb2=(x43*x32+y43*y32+z43*z32)/(rsp43*rsp32)
          cb3=(x21*x43+y21*y43+z21*z43)/(rsp21*rsp43)
          sb1=DSQRT(abs(1.0d0-cb1**2))
          sb2=DSQRT(abs(1.0d0-cb2**2))
          sb3=DSQRT(abs(1.0d0-cb3**2))
          aux=sb1*sb2
          if(abs(aux*sb3).lt.1.d-25) THEN
c--           Verbose mode only in dynamic not in minimization
            if(lwrite) THEN 
              if(.not.screw_torsion) WRITE(6,776) i,l1,l2,l3,l4
776           FORMAT
     &             (' ****** WARNING(fptors)-  Angle singularity:',i10,/
     &             ' ****** is the first torsion to screw ',4i5)
            END IF
            screw_torsion=.true.
            nscrew=nscrew+1
            if(minimize.AND.nscrew.gt.(ltors+500)) THEN
              write(6,*) "Looping in fptors.",
     &             "Program STOPS to avoid NaN "
              STOP
            endif
c           slightly rattle coordinates...
            xr1=xr1-0.005+ 0.01*ranf(dummy)
            yr1=yr1-0.005+ 0.01*ranf(dummy) 
            zr1=zr1-0.005+ 0.01*ranf(dummy)
            xr2=xr2-0.005+ 0.01*ranf(dummy)
            yr2=yr2-0.005+ 0.01*ranf(dummy)
            zr2=zr2-0.005+ 0.01*ranf(dummy)
            xr3=xr3-0.005+ 0.01*ranf(dummy)
            yr3=yr3-0.005+ 0.01*ranf(dummy)
            zr3=zr3-0.005+ 0.01*ranf(dummy)
            xr4=xr4-0.005+ 0.01*ranf(dummy)
            yr4=yr4-0.005+ 0.01*ranf(dummy)
            zr4=zr4-0.005+ 0.01*ranf(dummy)
            GOTO 1009
          ENDIF 
          coa=(cb1*cb2-cb3)/aux
          tmp=IDNINT(pota(i))
          ntph=IABS(tmp)
          sig=DFLOAT(tmp/ntph)
          IF(coa.GT.1.0d0.AND.coa.LT.1.000001) coa=1.0d0
          IF(coa.LT.-1.0d0.AND.coa.GT.-1.000001) coa=-1.0d0
          aux1=cb1*sb2
          aux2=sb1*cb2
          auxa=-(aux1*coa+aux2)/aux
          auxb=-(aux2*coa+aux1)/aux
          auxc=sb3/aux
          aux1=rsp21*rsp32*sb1
          aux2=rsp32*rsp43*sb2
          aux3=rsp21*rsp43*sb3
          dbx11= x32/aux1-cb1*x21/(rsq21*sb1)
          dbx13=-x21/aux1+cb1*x32/(rsq32*sb1)
          dbx12=-dbx13-dbx11
          dbx22= x43/aux2-cb2*x32/(rsq32*sb2)
          dbx24=-x32/aux2+cb2*x43/(rsq43*sb2)
          dbx23=-dbx22-dbx24
          dbx31=x43/aux3-cb3*x21/(rsq21*sb3)
          dbx33=x21/aux3-cb3*x43/(rsq43*sb3)
          dbx32=-dbx31
          dbx34=-dbx33
          dby11= y32/aux1-cb1*y21/(rsq21*sb1)
          dby13=-y21/aux1+cb1*y32/(rsq32*sb1)
          dby12=-dby13-dby11
          dby22= y43/aux2-cb2*y32/(rsq32*sb2)
          dby24=-y32/aux2+cb2*y43/(rsq43*sb2)
          dby23=-dby22-dby24
          dby31=y43/aux3-cb3*y21/(rsq21*sb3)
          dby33=y21/aux3-cb3*y43/(rsq43*sb3)
          dby32=-dby31
          dby34=-dby33
          dbz11= z32/aux1-cb1*z21/(rsq21*sb1)
          dbz13=-z21/aux1+cb1*z32/(rsq32*sb1)
          dbz12=-dbz13-dbz11
          dbz22= z43/aux2-cb2*z32/(rsq32*sb2)
          dbz24=-z32/aux2+cb2*z43/(rsq43*sb2)
          dbz23=-dbz22-dbz24
          dbz31=z43/aux3-cb3*z21/(rsq21*sb3)
          dbz33=z21/aux3-cb3*z43/(rsq43*sb3)
          dbz32=-dbz31
          dbz34=-dbz33
          IF(ntph.EQ.1) THEN
              enrg=sig*(DABS(potb(i))-coa*potb(i))
              qforce=sig*potb(i)
          ELSE IF(ntph.EQ.2) THEN
              enrg=DABS(potb(i))-(2.0d0*coa**2-1.0d0)*potb(i)
              enrg=sig*enrg
              qforce=4.0d0*potb(i)*coa
              qforce=sig*qforce
          ELSE IF(ntph.EQ.3) THEN
              enrg=DABS(potb(i))-(4.0d0*coa**3-3.0d0*coa)*potb(i)
              enrg=sig*enrg
              qforce=(12.0d0*coa**2-3.0d0)*potb(i)
              qforce=sig*qforce
          ELSE IF(ntph.EQ.4) THEN
              enrg=DABS(potb(i))-
     x                   (8.0D0*coa**4-8.0D0*coa**2+1.0D0)*potb(i)
              enrg=sig*enrg
!             write(6,1001) 180*dacos(coa)/3.1459,enrg/abs(potb(i))
              qforce=(32.0D0*coa**3-16.0D0*coa)*potb(i)
              qforce=sig*qforce
          ELSE IF(ntph.EQ.5) THEN
              enrg=DABS(potb(i))-
     x                   (16.0D0*coa**5-20.0D0*coa**3+5.0D0*coa)*potb(i)
              enrg=sig*enrg
              qforce=(80.0D0*coa**4-60.0D0*coa**2 + 5.d0)*potb(i)
              qforce=sig*qforce
          ELSE IF(ntph.EQ.6) THEN
              enrg=DABS(potb(i))-potb(i)*
     x        (32.0D0*coa**6-48.0D0*coa**4+18.0D0*coa**2-1)
              enrg=sig*enrg
                   qforce=potb(i)*
     x        (192.0D0*coa**5-192.0D0*coa**3+36.0D0*coa)
              qforce=sig*qforce
          ELSE IF(ntph.EQ.7) THEN
              enrg=DABS(potb(i))-potb(i)*
     x        (64.0D0*coa**7-112.0D0*coa**5+56.0D0*coa**3-7.d0*coa)
              enrg=sig*enrg
                   qforce=potb(i)*
     x        (448.0D0*coa**6-560*coa**4+168.D0*coa**2-7.d0)
              qforce=sig*qforce
          ELSE IF(ntph.EQ.8) THEN
              enrg=DABS(potb(i))-potb(i)*
     x       (128.0D0*coa**8-256*coa**6+160.0D0*coa**4-32.d0*coa**2 + 1)
              enrg=sig*enrg
                   qforce=potb(i)*
     x        (1024D0*coa**7-1536*coa**5+640.D0*coa**3-64.d0*coa)
              qforce=sig*qforce
          ELSE IF(abs(ntph).EQ.14) THEN
!             phase is 60 and n=4
              bb=DACOS(coa)
              XN = (y21*z32-y32*z21)
              YN = (z21*x32-z32*x21)
              ZN = (x21*y32-x32*y21) 
              aux4=sign(bb,XN*x43+YN*Y43+ZN*z43)
              soa=sign(dsqrt(1.d0-coa**2.d0),aux4)
              dsoa=-coa/soa
              enrg=potb(i)*(1.d0 +
     x           0.5*(8.0D0*coa**4-8.0D0*coa**2+1.0D0))
              auxn = potb(i)*sin60*4.d0*(soa*coa**3-coa*soa**3)
              dauxn= potb(i)*sin60*4.d0*(dsoa*(coa**3-3.d0*coa*soa**2) -
     &             soa**3+3.d0*soa*coa**2)
              enrg=enrg + auxn  
              enrg=sig*enrg
              qforce=0.5d0*(32.0D0*coa**3-16.0D0*coa)*potb(i)  
              qforce= -(qforce + dauxn)
              qforce=sig*qforce
          ELSE IF(abs(ntph).EQ.12) THEN
!             phase is 60 and n=2
              bb=DACOS(coa)
              XN = (y21*z32-y32*z21)
              YN = (z21*x32-z32*x21)
              ZN = (x21*y32-x32*y21) 
              aux4=sign(bb,XN*x43+YN*Y43+ZN*z43)
              soa=sign(dsqrt(1.d0-coa**2.d0),aux4)
              dsoa=-coa/soa
              enrg=potb(i)*(1.d0 +
     x           0.5*(2.d0*coa**2+1.0D0))
              auxn = potb(i)*sin60*2.d0*soa*coa
              dauxn= potb(i)*sin60*2.d0*(soa+dsoa*coa)
              enrg=enrg + auxn  
              enrg=sig*enrg
              qforce=2.0D0*coa*potb(i)  
              qforce= -(qforce + dauxn)
              qforce=sig*qforce
          ELSE
              WRITE(6,'(a)') '====== ABORT Multiplicity parameter is ===
     x==='
              WRITE(6,'(10x,18h NTPH is equal to ,i5)') ntph
              WRITE(6,'(2e12.5)') pota(i),potb(i)
              STOP
          END IF
          utors(type)=utors(type)+enrg
          uux1=auxa*dbx11+auxc*dbx31
          uux2=auxa*dbx12+auxb*dbx22+auxc*dbx32
          uux3=auxa*dbx13+auxb*dbx23+auxc*dbx33
          uux4=auxb*dbx24+auxc*dbx34
          uuy1=auxa*dby11+auxc*dby31
          uuy2=auxa*dby12+auxb*dby22+auxc*dby32
          uuy3=auxa*dby13+auxb*dby23+auxc*dby33
          uuy4=auxb*dby24+auxc*dby34
          uuz1=auxa*dbz11+auxc*dbz31
          uuz2=auxa*dbz12+auxb*dbz22+auxc*dbz32
          uuz3=auxa*dbz13+auxb*dbz23+auxc*dbz33
          uuz4=auxb*dbz24+auxc*dbz34
          fpx(l1)=fpx(l1)+qforce*uux1
          fpx(l2)=fpx(l2)+qforce*uux2
          fpx(l3)=fpx(l3)+qforce*uux3
          fpx(l4)=fpx(l4)+qforce*uux4
          fpy(l1)=fpy(l1)+qforce*uuy1
          fpy(l2)=fpy(l2)+qforce*uuy2
          fpy(l3)=fpy(l3)+qforce*uuy3
          fpy(l4)=fpy(l4)+qforce*uuy4
          fpz(l1)=fpz(l1)+qforce*uuz1
          fpz(l2)=fpz(l2)+qforce*uuz2
          fpz(l3)=fpz(l3)+qforce*uuz3
          fpz(l4)=fpz(l4)+qforce*uuz4
10      CONTINUE
#ifdef _OMP_
        endif
#endif
      end do
      if(lwrite.and.nscrew.gt.0.and.itask.eq.1) WRITE(6,103) nscrew
103   FORMAT( 
     & ' ****** Total number of singularities: ',i8)
      utors_slt=utors(1)
      utors_slv=utors(2)

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
