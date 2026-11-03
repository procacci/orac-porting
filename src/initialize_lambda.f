!--------------------written by Procacci 2013---------------------------------
      subroutine initialize_lambda(nat_added,nat_removed,ladded,laddedq
     &     ,lremoved,lremovedq,atom_added,atom_removed,lambda,lambdaq
     &     ,listqq,nlistqq,lstretch,lstrtch,lbend,lbndg,int14p,int14
     &     ,lconstr,lcnstr,clewld,ss_alchemy)
!-----------------------------------------------------------------------------

*======================== DECLARATIONS ================================*

      IMPLICIT none


*----------------------- ARGUMENTS  OR EXTERNAL ------------------------

      INTEGER nat_added,nat_removed,atom_added(*),atom_removed(*)
     &     ,lstretch,lstrtch(2,*),lconstr,lcnstr(2,*),lbend,lbndg(3,*)
     &     ,int14p,int14(2,*),nlistqq,listqq(3,*),ss_alchemy(*)

      REAL*8 ladded(*),lremoved(*),laddedq(*),lremovedq(*),lambda(*)
     &     ,lambdaq(*)
      
      EXTERNAL  near0
      LOGICAL near0,clewld

*-------------------LOCAL VARIABLES ------------------------------------

      INTEGER  i,j,l, iatm,jatm,ia,ib,iat(2000),type(2000)
      CHARACTER*80 errmsg
      logical lskip
      
      
      
      if(nat_added.NE.0) THEN 
        do i=1,nat_added
          ss_alchemy(atom_added(i))=2
          lambda(atom_added(i))=1.d0
          lambdaq(atom_added(i))=1.d0
        end do
        if(.not.near0(abs(ladded(1)-1.d0)). and. .not.
     &       near0(abs(laddedq(1)-1.d0))) THEN
          errmsg='In alchemical transformation: initial value'
     &         //
     &         ' of lambda for added atoms is not 1. check alchemy file'
          CALL xerror(errmsg,80,1,1)
        end if
      ENDIF
      if(nat_removed.NE.0) THEN 
        do i=1,nat_removed
          ss_alchemy(atom_removed(i))=-2
          lambda(atom_removed(i))=0.d0
          lambdaq(atom_removed(i))=0.d0
        end do
        if(.not.near0(abs(lremoved(1))). and .
     &       .not.near0(abs(lremovedq(1)))) THEN
          errmsg='In alchemical transformation: initial value'
     &    // 
     &    ' of lambda for removed atoms is not 0. check alchemy file'
          CALL xerror(errmsg,80,1,1)
        end if
      ENDIF

!     prepare a list of dq dq intralchemical interactions to be 
!     subtracted out in mts_furier_alchemy (ferrf_alchemy)  

      if(.not.clewld) RETURN 


      do i=1,nat_added
        iat(i)=atom_added(i) 
        type(i)=1
      end do

      do i=1,nat_removed
        iat(i+nat_added)=atom_removed(i) 
        type(i+nat_added)=0
      end do

      
      nlistqq=0
      do i=1,nat_added+nat_removed-1
        iatm =iat(i)
        do j=i+1,nat_added+nat_removed 
          jatm= iat(j) 
          lskip=.false.
c---      check whether iat-jat is a 1-2 interaction
          do l=1,lstretch
            ia=lstrtch(1,l)
            ib=lstrtch(2,l)
            if(iatm.eq.ia.and.jatm.eq.ib.or.iatm.eq.ib.and.jatm.eq.ia)
     &           THEN 
              lskip=.true.
              exit 
            end if 
          end do
c---      check whether iat-jat is a 1-3 interaction
          if(lskip) goto 101
          do l=1,lconstr
            ia=lcnstr(1,l)
            ib=lcnstr(2,l)
            if(iatm.eq.ia.and.jatm.eq.ib.or.iatm.eq.ib.and.jatm.eq.ia)
     &           THEN 
              lskip=.true.
              exit 
            end if 
          end do
c---      check whether iat-jat is a 1-3 interaction
          if(lskip) goto 101
          do l=1,lbend
            ia=lbndg(1,l)
            ib=lbndg(3,l)
            if(iatm.eq.ia.and.jatm.eq.ib.or.iatm.eq.ib.and.jatm.eq.ia)
     &           THEN 
              lskip=.true.
              exit 
            end if 
          end do
c---      check whether iat-jat is a 1-4 interaction
          if(lskip) goto 101
          do l=1,int14p
            ia=int14(1,l)
            ib=int14(2,l)
            if(iatm.eq.ia.and.jatm.eq.ib.or.iatm.eq.ib.and.jatm.eq.ia)
     &           THEN 
              lskip=.true.
              exit 
            end if 
          end do
!   and if so skips the interaction
          if(.not.lskip) THEN
            nlistqq=nlistqq+1 
            listqq(1,nlistqq)= iatm
            listqq(2,nlistqq)= jatm
            listqq(3,nlistqq)= 1
!           removed species do not feel added species
            if(type(iatm).ne.type(jatm)) listqq(3,nlistqq)=0
          end if
101       continue
        end do
      end do
      RETURN 
      end 
