      SUBROUTINE rdcmac(pdb_unit,kprint,res,beta,x0,y0,z0,nato
     &     ,multiple,mol,nmol,iret,errmsg)
************************************************************************
*      Read pdb file                                                   *
************************************************************************

*======================= DECLARATIONS ==================================

      use parst
      use cpropar, only:mend,prsymb,read_remark,co,oc,volume,pdbfile
      use rem, only:cdist,para_index  
      
      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER pdb_unit,kprint,nato,res(nato),iret,mol,nmol
      REAL*8 x0(*),y0(*),z0(*)
      CHARACTER*7 beta(*)
      CHARACTER*80 errmsg
      LOGICAL multiple

      INTEGER  mmm
      parameter     (mmm=nores+slvunit) 

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8   xtp,ytp,ztp,xa,ya,za
      INTEGER  ires,i,ires_old,res_sequence,pdb_res,n,nwx,nwy,nwz,iline
     &     ,nword,idot,dot(3),pdb_res_old,ibeg(mmm),iend(mmm),jatm,j,is
     &     ,nwx_first,nwy_first,nwz_first,ibgdt(3),ieddt(3),nerror,i1,i2
     &     ,nh,pdb_res_count,ialias
      LOGICAL*1     lok,found,ltmp(sitslu+tsites)
      LOGICAL ok
      CHARACTER*10  integers
      CHARACTER*5   field(3)
      CHARACTER*4   beta4,char,charold 
      CHARACTER*180 line
      CHARACTER*1  beta1

*-------------------- DECLARATION OF A SCRATCH COMMON ------------------

      COMMON /dynam/ ibeg,iend,ltmp

c=======================================================================
c---  finds the start-of-residue end-of-residue pointers in 
c---  topologic array. Also find total expected number of residue 
c=======================================================================

      ires = 0
      res_sequence=0
      nh = 0
      jatm=0
      do i=1,nato
c---     this coordinates are overwritten if hydrogen are in the pdb
         beta1=beta(i)
	 if(beta1.eq.'h') then 
	     nh = nh + 1
	     x0(i)=1.d10 
	     y0(i)=1.d10 
	     z0(i)=1.d10
	 end if 
         ltmp(i)=.false.
         ires_old = ires
         ires=res(i)
         if(ires.ne.ires_old) then 
            res_sequence = res_sequence + 1
            ibeg(res_sequence) = i 
            if (i.gt.1) THEN 
               iend(res_sequence-1) = i-1
            END IF
         END IF
      END DO 

      iend(res_sequence) = nato

      pdb_res= 1
      pdb_res_count=1
      ok=.FALSE.

c======================================================================
c     Start parsing the PDB file 
c======================================================================

      nwx=0
      nwy=0
      nwz=0
      ibgdt(1)=31
      ieddt(1)=38
      ibgdt(2)=39
      ieddt(2)=46
      ibgdt(3)=47
      ieddt(3)=54
      field(1) = '31-38'
      field(2) = '39-46'
      field(3) = '47-45'
      integers = '0123456789' 
      nerror=0
       
      n=0
      iline=0

      IF(mol .EQ. 1) write(kprint,'(5x,a,a)') 
     &     'Reading the PDB file               ---->',trim(pdbfile)

 1    iline = iline + 1 
c-----   read line and bring to lowercase
         READ(pdb_unit,'(a180)',END=2,ERR=1001) line
         CALL up_low(line,80)
         IF(line(1:3).eq."end") goto 2
         IF(cdist.eq.3.and.(line(1:7).eq."remark")) THEN
           write(kprint,'(/,5x,a)')
     &     'Reading CO Matrix and REM index from PDB file     ---->'
           read(line,'(42x,9f10.5,11x,i4)',err=1001) co(1,1),co(1,2)
     &          ,co(1,3),co(2,1),co(2,2),co(2,3),co(3,1),co(3,2),co(3,3)
     &          ,para_index
           do i=1,3
             do j=1,3
               co(i,j)=0.5*co(i,j)
             end do
           end do
           CALL matinv(3,3,co,oc,volume)
           write(kprint,'(/,15x,a)') 
     &          ' ---  CO Matrix ---      ---  CO Matrix ---'
           write(kprint,'(3f10.5,5x,3f10.5)') ((co(i,j),j=1,3),(oc(i,j)
     &          ,j=1,3),i=1,3)
           write(kprint,'(" Para_index =", i4 )') para_index 
           write(6,*) line
           goto 8001
         end if
         IF(read_remark.and.(line(1:7).eq."remark")) THEN
           write(kprint,'(/,5x,a)')
     &     'Reading CO Matrix from PDB file    ---->'
           read(line,'(42x,9f10.5)',err=1001) co(1,1),co(1,2),co(1,3)
     &          ,co(2,1),co(2,2),co(2,3),co(3,1),co(3,2),co(3,3)
           do i=1,3
             do j=1,3
               co(i,j)=0.5*co(i,j)
             end do
           end do
           CALL matinv(3,3,co,oc,volume)
           write(kprint,'(/,15x,a)') 
     &          ' ---  CO Matrix ---      ---  CO Matrix ---'
           write(kprint,'(3f10.5,5x,3f10.5)') ((co(i,j),j=1,3),(oc(i,j)
     &          ,j=1,3),i=1,3)
         end if
8001     IF(line(1:4) .eq. 'atom' .or. line(1:6) .eq. 'hetatm'
     &        ) THEN 

c======================================================================
c           Check line to see if .pdb 
c======================================================================

c===========check position of (.) character for coordinates============

            idot=0
            dot(1) = 0
            dot(2) = 0
            dot(3) = 0

c---        if position are wrong, abort

            do is=31,54
               if(line(is:is).eq.'.') then 
                  idot = idot+1 
                  if(idot.gt.3) THEN 
                     write(kprint,100) iline
                     nerror=nerror+1       
                  end if 
                  dot(idot) = is
               endif   
            end do

            if(idot.lt.3) THEN 
               write(kprint,210) iline
               nerror=nerror+1       
            end if 

c---        if position are misaligned give warning 

            do i=1,3
               lok = dot(i).ge.ibgdt(i).and.dot(i).le.ieddt(i)
               IF (.not.lok) THEN 
                  write(kprint,130) iline,field(i)
                  nerror=nerror+1       
               END IF
            END DO
            IF(line(35:35).ne.'.') THEN
               nwx = nwx + 1
               if (nwx.eq.1) nwx_first=iline
            ENDIF 
            IF (line(43:43).ne.'.') THEN 
               nwy = nwy + 1
               if (nwy.eq.1) nwy_first=iline
            ENDIF 
            IF (line(51:51).ne.'.') THEN 
               nwz = nwz + 1
               if (nwz.eq.1) nwz_first=iline
            ENDIF 

c===========check if field 26-26 is a number================================= 

            lok = .false.
            do i=1,10
               if(line(26:26).eq.integers(i:i)) THEN
                  lok = .true. 
               end if
            end do
            IF (.not.lok) THEN 
               write(kprint,120) iline
               STOP
            END IF

            if(nerror.gt.100) THEN 
               write(kprint,310) 
               STOP
            END IF
               
c----       tests passed: Good chance that line is in genuin pdb format  

c======================================================================
c           Now do internal read of LINE in PBD format
c======================================================================

            IF(ok) pdb_res_old = pdb_res
            READ(line,4,err=3) char,pdb_res,xtp,ytp,ztp 
            IF(.NOT. ok) pdb_res_old=pdb_res
            ok=.TRUE.

            IF (pdb_res .NE. pdb_res_old)THEN
               pdb_res_count=pdb_res_count+1
            END IF
            
c----       check if .pdb residue number is something meaningful

            if(pdb_res_count .GT. res_sequence) then
               IF(.NOT. multiple) THEN
                  write(kprint,10) iline,pdb_res_count,res_sequence
                  STOP
               END IF
               BACKSPACE pdb_unit
               GOTO 2
            END IF

c--         justify left character string 

            call jusleft(char)

c---        Now we check the Residue sequence of the .pdb file. 
c----       check if CHAR matches some BETA in residue ires

            ialias=0
            charold=char 
25          found = .false.
            do j=ibeg(pdb_res_count),iend(pdb_res_count) 
               beta4=beta(j)
               if (char.eq.beta4) then 
                  found = .true. 

c----             j-th atom of residue pdb_res matches!!!!

                  if(.not.ltmp(j)) then 
                     ltmp(j)=.true.
                   else

c----                 unique labels must be specified in the .pdb 

                      write(kprint,75) iline,char
                      STOP
                   END IF 
                  go to 101
               end if
            end do

!           now try aliases  
            if(ialias.EQ.0) THEN 
              call alias(char,prsymb(mend(pdb_res_count)))
              ialias=1
              go to 25 
            END IF
            
c---        AAAAArgh label not found: write error message
            
            write(kprint,20) iline,line,charold,pdb_res_count
     &           ,prsymb(mend(pdb_res_count))
            STOP
 101        continue
            if (found) then 
               n = n + 1 
               x0(j)=xtp
               y0(j)=ytp
               z0(j)=ztp
            end if   

c---        if the atom was an hydrogen decrement number of h atoms in
c---        tipological sequence and increment nh of pbd sequence            

            beta1=char
            if(beta1.eq.'h') then 
               nh = nh-1
            end if
            go to 1
         END IF
         if(pdb_res_count .le.res_sequence) go to 1
      CONTINUE

 2    continue

c--   check for bad alignment 

      IF(nerror.gt.0) THEN
         write(kprint,66)  nerror
         STOP
      END IF

c---  check number of atoms 
      
      n = n + nh 

      IF (n.ne.nato.and.n.gt.0) THEN  
         write(kprint,40) n,nato,n-nh,nh
         write(kprint,3100) 
3100     FORMAT(//4x," ----- List of unassigned heavy atoms ------",/
     &    ,12x,"atom",2x,"label",9x,"ires",4x,"res")     
         j=0
         do i=1,nato
           beta1=beta(i)
           if(beta1.NE."h") j=j+1
           if(.not.ltmp(i)) THEN
             if(beta1.NE."h") THEN
               write(kprint,3101) j,beta(i),res(i),prsymb(mend(res(i)))
3101           FORMAT(6x,i9,5x,a7,2x,i5,5x,a7)
             END IF
           END IF
         end do
         STOP
      END IF
      IF (pdb_res_count.ne.res_sequence) THEN 
         write(kprint,45) pdb_res_count, res_sequence
         STOP
      END IF


      IF (n.eq.0) THEN
         write(kprint,50) 
         STOP 
      END IF

c======================================================================
c     Atoms number is OK; trasform coordinates and exit; 
c======================================================================

      IF(mol .EQ. nmol) write(kprint,290) n*nmol,nh

c======================================================================
c     Signals misaligned floating point in .pdb coordinates
c======================================================================

c---  jump here because of internal reading error; set iret to 1
c---  prints lines where possible misalignment of reals could 
c---  have been found and exit 

      IF(nwx.gt.0) THEN
         WRITE (kprint,65) 
         write(kprint,70) nwx,nwx_first
         WRITE (kprint,65) 
      END IF
      IF(nwy.gt.0) THEN
         WRITE (kprint,65) 
         write(kprint,76) nwy,nwy_first
         
         WRITE (kprint,65) 
      END IF
      IF(nwz.gt.0) THEN
         WRITE (kprint,65) 
         write(kprint,80) nwz,nwz_first
         
         WRITE (kprint,65) 
      END IF
      write(kprint,'(5x,"=",15x,a,23x,"=",/5x,75("="),/)') 
     &     "Successful reading of the pdb file " 
      RETURN

 3    write(kprint,55) iline
      STOP

 1001 write(kprint,60) iline
*======================= END statements  ================================
      STOP
4     FORMAT(12x,a4,6x,i4,4x,3f8.3)
10    FORMAT(//'****ERROR IN PDB FILE at LINE: ',i6,
     &     '; RESIDUE NUMBER IS ',i5,' BEYOND EXPECTED NUMBER',i6/    
     &     '; MAYBE  MISSING OR MISPLACED <END> KEYWORD?  ')    
20    FORMAT(//'****ERROR IN PDB FILE at LINE: ',i6,/,A80,/
     &     ' LABEL --> ',A4,' OF RESIDUE',I5," "A6,' UNKNOWN'/)
40    FORMAT(//'****ERROR IN PDB FILE:', / 
     &  '* The total numer of atoms read in the PDB was ',i6,/
     &  '* while the expected number of atoms according to the tpg',
     &  '  file was ',i6,/ 
     &  '* Number of found heavy atoms are ', i6,/
     &  '* Number of (added) hydrogen atoms are ', i6) 
45         FORMAT(//'****ERROR IN PDB FILE:', / 
     &   ' The total numer of residues in the PDB was ',i4,/
     &   ' while the expected number of residues according', 
     &   ' to the input file was ',i4)
50    FORMAT(//'****ERROR IN PDB FILE: no ATOM or HETATM',
     &     ' string found')
55    FORMAT(/'****ERROR in PDB: internal reading error at line:',i5)
60    FORMAT('****ERROR in reading PDB at line:',i5) 
65    FORMAT('***************************************************')
66    FORMAT( / '-----  Stop while reading the PDB file:',    
     &     / '----- ',i3,' fatal errors were found')
70    FORMAT('***WARNING:IN PDB FILE X coordinates are misaligned:'
     &     , /,'   dot (.) char should be at column 35; I found ',i5,
     &     /,'   lines with misaligned X, starting at line: ',i5)     
      
75    FORMAT(//'****ERROR IN PDB FILE at LINE: ',i6,
     &     '; ATOM ',a4, ' already assigned')
76    FORMAT('***WARNING:IN PDB FILE Y coordinates are misaligned:'
     &     ,/,'   dot (.) char should be at column 35; I found ',i5,
     &     /,'   lines with misaligned Y, starting at line: ',i5)     
80    FORMAT('***WARNING:IN PDB FILE Z coordinates are misaligned:'
     &     ,/'   dot (.) char should be at column 35; I found ',i5,
     &     /,'   lines with misaligned Z, starting at line: ',i5)     
100   FORMAT('****ERROR IN PDB FILE at LINE: ',i6,  
     &     ' unexpected FORMAT: more than 3 reals')
120   FORMAT('****ERROR IN PDB FILE at LINE: ',i6,  
     &     ' - at column 26 a number is expected')
130   FORMAT('****ERROR IN PDB FILE at LINE: ',i6,  
     &     ' Character (.) not found in field ', A5)
210   FORMAT('****ERROR IN PDB FILE at LINE: ',i6,  
     &     ' unexpected FORMAT: less  than 3 reals')
280   FORMAT
     &     (/'                    Reading PDB file.......') 
310   FORMAT(// ' --- more than 100 errors. Check alignment'/)
      
290   FORMAT
     &     (//5x,75("="),
     &      /,5x,'=',4x,'PDB file parsed:',i6,' atoms; - ',i6
     &     ,' hydrogens to be',' assigned      =')

      END
************************************************************************
      SUBROUTINE jusleft(STRING)
      implicit none 
C     justify to the left 4 character string  
*======================= DECLARATIONS ==================================
      CHARACTER*4 STRING
      INTEGER  i,istart 
*======================= Exec statements  ==============================
      i=0
      istart=1
 1    i=i+1
      if (string(i:i).ne.' ') then
         istart = i
         go to 2
      end if
      go to 1
2     continue
      string = string(istart:4)		   	 	 

*======================= END statements  ================================

      RETURN
      END 

      subroutine alias(char,res) 
      implicit none 


      character*4 char,charold 
      character*7 res
      
!     protein atoms alias list follows
      charold=char
!===================General aliases
      if(char.eq."ht1")THEN
        char="h1"               ! n-termini AMBER
      else if(char.eq."h1")THEN
        char="ht1"              ! n-termini OPLS
      end if
      if(char.eq."ht2")THEN
        char="h2"
      else if(char.eq."h2") THEN
        char="ht2"
      end if
      if(char.eq."ht3")THEN
        char="h3" 
      else if(char.eq."h3")THEN
        char="ht3"
      end if
      if(char.eq."hn1") char="ht1"   ! n-termini
      if(char.eq."hn2") char="ht2"
      if(char.eq."hn3") char="ht3" 
      if(char.eq."ot1") THEN
        char="o"                ! c-termini 
      else IF (char.eq."o") THEN
        char="ot1"              ! c-termini
      end if
      if(char.eq."ot2") THEN
        char="oxt"                ! c-termini
      else IF (char.eq."oxt") THEN
        char="ot2"              ! c-termini
      end if
      if(char.eq."hn")  THEN
        char="h"                ! backbone AMBER
      else
        if(char.eq."h")   THEN
          char="hn"             ! backbone OPLS
        end if
      end if
!===================other  aliases
      if(res(1:3).eq."glu".or.res(1:3).eq."trp".or.res(1:3)
     &   .eq."arg".or.res(1:3)
     &   .eq."gln".or.res(1:3).eq."leu".or.res(1:3).eq."asp".or.res(1:3)
     &   .eq."asn".or.res(1:3).eq."tyr".or.res(1:3).eq."met".or.res(1:3)
     &   .eq."pro".or.res(1:3).eq."ser".or.res(1:3).eq."lys".or.res(1:3)
     &   .eq."hsd".or.res(1:3).eq."hse".or.res(1:3).eq."hsp".or.res(1:3)
     &   .eq."phe" ) THEN 
        if(char.eq."hb3") char="hb1" 
        if(char.eq."hd3") char="hd1" 
        if(char.eq."he3") char="he1" 
        if(char.eq."hg3") char="hg1"  
      end if
      if(res(1:3).eq."cys".OR.res(1:3).eq."ser") THEN 
        if(char.eq."hg1") THEN
          char="hg"
        else if(char.eq."hg") THEN
          char="hg1"
        end if
      endif
      if(res(1:3).eq."arg") THEN
        if(char.eq."1hh1") THEN
          char="hh11"
        endif
        if(char.eq."2hh1") THEN
          char="hh12"
        endif
        if(char.eq."2hh2") THEN
          char="hh22"
        endif
        if(char.eq."1hh2") THEN
          char="hh21"
        endif
      endif
      if(res(1:3).eq."ile") THEN   
        if(char.eq."cd")  THEN
          char="cd1"
        endif
        if(char.eq."hd1") THEN
          char="hd11"
        else if(char.eq."hd11") THEN
          char="hd1"
        end if
        if(char.eq."hd2") THEN
          char="hd12"
        else if(char.eq."hd12") THEN
          char="hd2"
        end if
        if(char.eq."hd3") THEN
          char="hd13"
        else if(char.eq."hd13") THEN
          char="hd3"
        end if
      END IF
      if(res(1:3).eq."val") THEN
        if(char.eq."1hg1") char="hg11"
        if(char.eq."2hg1") char="hg12"
        if(char.eq."3hg1") char="hg13"
        if(char.eq."1hg2") char="hg21"
        if(char.eq."2hg2") char="hg22"
        if(char.eq."3hg2") char="hg23"
      END IF
      if(res(1:3).eq."gln") THEN 
        if(char.eq."1he2") char="he21"
        if(char.eq."2he2") char="he22"
      END IF
      if(res(1:3).eq."thr") THEN
        if(char.eq."1hg2") char="hg21"
        if(char.eq."2hg2") char="hg22"
        if(char.eq."3hg2") char="hg23"
      END IF
      if(res(1:3).eq."asn") THEN
        if(char.eq."1hd2") char="hd21"
        if(char.eq."2hd2") char="hd22"
        if(char.eq."3hd2") char="hd23"
      END IF
      if(res(1:3).eq."leu") THEN
        if(char.eq."1hd1") char="hd11"
        if(char.eq."2hd1") char="hd12"
        if(char.eq."3hd1") char="hd13"
        if(char.eq."1hd2") char="hd21"
        if(char.eq."2hd2") char="hd22"
        if(char.eq."3hd2") char="hd23"
      END IF
      if(res(1:3).eq."ile") THEN
        if(char.eq."1hg1") char="hg11"
        if(char.eq."2hg1") char="hg12"
        if(char.eq."1hg2") char="hg21"
        if(char.eq."2hg2") char="hg22"
        if(char.eq."3hg2") char="hg23"
        if(char.eq."1hd1") char="hd11"
        if(char.eq."2hd1") char="hd12"
        if(char.eq."3hd1") char="hd13"
      END IF
      if(res(1:3).eq."gly") THEN 
        if(char.eq."ha3") char="ha1"
      end if
!      write(6,10) charold,char,res
!10    format("name ", a4, " aliased to ", a4, " in residue ",A7)  
      return
      end 
