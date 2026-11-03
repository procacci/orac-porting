*=======================================================================
      subroutine assign_agbnp_types(ntap,nbtype,epslj,sigm,betb
     &     ,pnbd1,pnbd2,tipo,kprint)
*=======================================================================
*  This code assigns agbnp atomic types. Receive orac types sjorg and  * 
*  ejorg (Lennard Jones parameters) and returns agbnp types as read in * 
*  form the agbnp.param library. Written by P. Procacci 2011           * 
*=======================================================================
      
*======================= DECLARATIONS ==================================

*----------------------- ARGUMENTS -------------------------------------

      CHARACTER*7 betb(*) 
      CHARACTER*1 sep(2),comm(2)
      INTEGER ntap,nbtype(*),tipo(*)
      REAL*8  epslj(*),sigm(*),pnbd1(*),pnbd2(*)

*-----------------------LOCAL VARIABLES --------------------------------
      INTEGER kuni
      LOGICAL  exist,lok,check1,check2
      CHARACTER*80 errmsg,line,strngs(40),filepath,agbnp_param_file
      CHARACTER*7 agbnp_type(200)
      CHARACTER*8 fmt
      DATA sep/' ',','/comm/'#','!'/
      REAL*8  diel_in, diel_out

      nex=0
      agbnp_param_file="./agbnp.param"
      INQUIRE(FILE=agbnp_param_file,EXIST=exist)
      IF(.NOT. exist) THEN
        errmsg=
     &       'File agbnp.param does not exist in current dir. Abort.'
        CALL xerror(errmsg,80,1,30)
        STOP
      ELSE
        CALL openf(kunit,agbnp_param_file,'FORMATTED','OLD',0)
      END IF
      diel_in=-1.d0
      diel_out=-1.d0
100   READ(kunit,'(a80)',END=600) line(1:80)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      if(strngs(1).EQ."dielectric_in") THEN 
       CALL fndfmt(2,strngs(2),fmt)
       READ(strngs(2),fmt,err=20) diel_in 
       GOTO 100
      end if
      if(strngs(1).EQ."dielectric_out") THEN 
       CALL fndfmt(2,strngs(2),fmt)
       READ(strngs(2),fmt,err=20) diel_out 
       GOTO 100
      end if
      if(nword.LT.10) goto 100
      nex = nex + 1
      agbnp_type(nex)=strngs(2)
      go to 100
600   write(kprint,1000) nex, diel_in,diel_out
1000  format(10x,
     &     '===================== AGBNP MINIMIZATION ===============',
     &/10x,'=                 ',6x,'                               =',
     &/10x,'=                 ',i6,' agbnp types found             =',
     &/10x,'=                 ',f10.3,' Dielectric in             =',
     &/10x,'=                 ',f10.3,' Dielectric out            =',
     &/10x,'=                 ',6x,'                               =')

c---  match orac types to agbnp types 
c      write(kprint,1001)
c1001  FORMAT(8x,"atom",2x,"agbnp-type",2x,"amber-type",9x,"eps",7x
c     &     ,"sigma")
      nsevere=0
      do i=1,ntap 
        lok=.false.
        tipo(i)=0
        do na=1,nex
          check1=TRIM(agbnp_type(na)).eq.TRIM(betb(nbtype(i)))
          check2=TRIM("*"//agbnp_type(na)).eq.TRIM(betb(nbtype(i))) 
          if(check1.or.check2) THEN
            tipo(i)=na
            epslj(i)= pnbd2(nbtype(i))
            sigm(i) = pnbd1(nbtype(i))/(2.d0)**(1.d0/6.d0) 
            lok=.true.
c            write(kprint,1002) i,tipo(i),nbtype(i),epslj(i),sigm(i)
c     &           ,agbnp_type(na),betb(nbtype(i)) 
c1002        format(3i12,2f15.4,2x,2A7) 
          END IF
        end do
        if(.not.lok) THEN
          errmsg=betb(nbtype(i))//" type unassigned; Program stop."  
          CALL xerror(errmsg,80,1,30)
          nsevere=nsevere+1
        END IF
      end do
      if(nsevere.gt.0) THEN 
        errmsg="unassigned types. Add unassigned types to agbnp.param"
        CALL xerror(errmsg,80,1,30)
        STOP
      ELSE
        WRITE(kprint,1003) 
1003    FORMAT(10x,
     &     '=     agbnp types assignment OK!                       ='
     &/10x,'========================================================')
      END IF  
      rewind(kunit)
      IF(diel_in.lt.0.d0.or.diel_out.LT.0.d0) THEN 
        errmsg=
     & "Dielectric constants undefined or negative. Check agbnp.param"
        CALL xerror(errmsg,80,1,30)
        STOP
      END IF
      RETURN
 20   CONTINUE
      iret=1
      errmsg='internal reading error: wrong format?? TAB character??'
      CALL xerror(errmsg,80,1,2)
      STOP
      END
             
