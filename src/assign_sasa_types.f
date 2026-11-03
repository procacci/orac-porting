*=======================================================================
      subroutine assign_sasa_types(ntap,nbtype,betb,kprint,sasa_nbtype
     &     ,sasa_betb,rvdw,psasa,ssasa,r2cut,rsolv)
*=======================================================================
*  This code assigns sasa atomic types. 
*=======================================================================
      
*======================= DECLARATIONS ==================================

      use unit, only:efact 
      IMPLICIT NONE

*----------------------- ARGUMENTS -------------------------------------

      CHARACTER*7 betb(*) 
      CHARACTER*1 sep(2),comm(2)
      INTEGER ntap,kprint,nbtype(*),sasa_nbtype(*)
      real*8  rvdw(*),psasa(*),ssasa(*),r2cut,rsolv

*-----------------------LOCAL VARIABLES --------------------------------
      INTEGER kuni,nex,i,iret,kunit,na,nsevere,nword
      LOGICAL  exist,lok,check1,check2
      CHARACTER*80 errmsg,line,strngs(40),filepath,sasa_param_file
      CHARACTER*7 sasa_betb(200)
      CHARACTER*8 fmt
      DATA sep/' ',','/comm/'#','!'/

      nex=0
      sasa_param_file="./sasa.param"
      r2cut=0.d0
      rsolv=1.4           ! default value for solvent Radius (Scheraga) 
      INQUIRE(FILE=sasa_param_file,EXIST=exist)
      IF(.NOT. exist) THEN
        errmsg=
     &       'File sasa.param does not exist in current dir. Abort.'
        CALL xerror(errmsg,80,1,30)
        STOP
      ELSE
        CALL openf(kunit,sasa_param_file,'FORMATTED','OLD',0)
      END IF
100   READ(kunit,'(a80)',END=600) line(1:80)
      IF(line(1:1) .EQ. '#') GOTO 100 
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      if(strngs(1).EQ."solvent_radius") THEN 
       CALL fndfmt(2,strngs(2),fmt)
       READ(strngs(2),fmt,err=20) rsolv
       GOTO 100
      end if
      if(nword.NE.5) goto 100 
      nex = nex + 1
      if(nex.gt.400) THEN 
        errmsg=" too many SASA types:  max is 400 "
        CALL xerror(errmsg,80,1,30)
        STOP
      END IF 
      sasa_betb(nex)=strngs(2)
      CALL fndfmt(2,strngs(3),fmt)
      READ(strngs(3),fmt,err=20) rvdw(nex)
      if((2.d0*rvdw(nex)+rsolv)**2.gt.r2cut) r2cut=(2*rvdw(nex)+rsolv)
     &     **2
      CALL fndfmt(2,strngs(4),fmt)
      READ(strngs(4),fmt,err=20) psasa(nex)
      CALL fndfmt(2,strngs(5),fmt)
      READ(strngs(5),fmt,err=20) ssasa(nex)
      ssasa(nex)=ssasa(nex)*4.184*1000.d0/efact   ! convert kcal to program units 
      go to 100
600   write(kprint,1000) nex
1000  format(10x,
     &     '===================== SASA MINIMIZATION ===============',
     &/10x,'=                 ',6x,'                               =',
     &/10x,'=                 ',i6,' sasa types found             =',
     &/10x,'=                 ',6x,'                               =')

c---  match orac types to sasa types 
c      write(kprint,1001)
c1001  FORMAT(8x,"atom",2x,"sasa-type",2x,"amber-type",9x,"eps",7x
c     &     ,"sigma")
      nsevere=0
      do i=1,ntap
        lok=.false.
        do na=1,nex
          check1=TRIM(sasa_betb(na)).eq.TRIM(betb(i))
          check2=TRIM("*"//sasa_betb(na)).eq.TRIM(betb(i)) 
          if(check1.or.check2) THEN 
            lok=.true.
            sasa_nbtype(i)=na
            exit 
          END IF
        end do
        if(.not.lok) THEN
          errmsg=betb(i)//" type unassigned; Program stop."  
          CALL xerror(errmsg,80,1,30)
          nsevere=nsevere+1
        END IF
      end do
      if(nsevere.gt.0) THEN 
        errmsg="unassigned types. Add unassigned types to sasa.param"
        CALL xerror(errmsg,80,1,30)
        STOP
      ELSE
        WRITE(kprint,1003) 
1003    FORMAT(10x,
     &     '=     sasa types assignment OK!                       ='
     &/10x,'========================================================')
      END IF  
      rewind(kunit)
      RETURN
 20   CONTINUE
      iret=1
      errmsg='internal reading error: wrong format?? TAB character??'
      CALL xerror(errmsg,80,1,2)
      STOP
      END
             
