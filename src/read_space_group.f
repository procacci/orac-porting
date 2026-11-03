      SUBROUTINE read_space_group(kprint,kgroup,cgroup,gnmol,cg,tg,iret
     &     ,errmsg)

************************************************************************
*   Time-stamp: <97/07/12 21:46:35 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Mon Feb 10 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*
      
      INTEGER kprint,kgroup,gnmol,iret
      REAL*8 cg(3,3,*),tg(3,*)
      CHARACTER*80 cgroup,errmsg

*----------------------- VARIABLES IN COMMON --------------------------*

*------------------------- LOCAL VARIABLES ----------------------------*

      CHARACTER*80 line,strngs(40),strngs1(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)

      INTEGER nword,nword1,iline,iok,nsevere,n,i,j
      DATA sep/' ',','/comm/'#','!'/

*----------------------- EXECUTABLE STATEMENTS ------------------------*


      WRITE(kprint,'(5x,a)')
     &     'Reading Space Group File           ---->'
      nsevere=0
      DO i=1,80
         line(i:i)=' '
      END DO
      
      iline = 0
100   CONTINUE
      READ(kgroup,'(a78)',END=600) line(1:78)
      iline = iline + 1
      IF(line(1:1) .EQ. '#') GOTO 100
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      
c========search for matching group======================================
c =
      
      IF(strngs(1).eq.'Space'.and.strngs(2).eq.'Group'.and.strngs(3)
     &     .eq.'Symmetry') THEN 
         
         CALL parse(cgroup,sep,2,comm,strngs1,40,nword1,iret,errmsg)
         iok=0
         do i=1,nword-3
            if(strngs(3+i).eq.strngs1(i)) iok = iok+1 
         end do
         if(iok.eq.(nword1)) THEN 
            READ(kgroup,'(a78)',END=600) line(1:78)
            iline = iline + 1 
            CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
            CALL fndfmt(1,strngs(1),fmt)
            READ(strngs(1),fmt,err=22) gnmol
c--------      group found; start loopin over molecules per cell
            do n=1,gnmol
c---              parsing rotation matrix for molecule n 
               DO i=1,3
                  READ(kgroup,'(a78)',END=600) line(1:78)
                  iline = iline + 1
                  CALL parse(line,sep,2,comm,strngs,40,nword
     &                 ,iret,errmsg)
                  DO J=1,3
                     if(nword.ne.3) THEN 
                        call int_str(iline,strngs(1),n) 
                        errmsg=
     &                       '*** ERROR while reading SPACE GROUP'//
     &                       'FILE atLINE: '//strngs(1)//'; Expecting'
     &                       //' 3 reals.'    
                     ELSE
                        CALL fndfmt(2,strngs(j),fmt)
                        READ(strngs(j),fmt,err=22) cg(i,j,n)
                     END IF
                  END DO 
               END DO
c---              parsing fractional translation
               READ(kgroup,'(a78)',END=600) line(1:78)
               CALL parse(line,sep,2,comm,strngs,40,nword
     &              ,iret,errmsg)
               iline = iline + 1
               if(nword.ne.3) THEN 
                  call int_str(iline,strngs(1),n) 
                  errmsg= '*** ERROR in SPACE GROUP FILE at'
     &                 //'LINE: '//strngs(1)/
     &                 /'; Expecting 3 reals'
               ELSE
                  DO j=1,3
                     CALL fndfmt(2,strngs(j),fmt)
                     READ(strngs(j),fmt,err=22) tg(j,n)
                  END DO
               END IF
            END DO
         ELSE 
            iok=0
            GO TO 100 
         END IF   
      ELSE
         go to 100
      END IF
600   IF(iok.eq.0) THEN 
         errmsg='*** ERROR in SPACE GROUP FILE: group not found' 
         call xerror(errmsg,80,1,2)
         STOP
      ELSE IF(iok.gt.0.and.nsevere.gt.0) THEN 
         errmsg='*** ERRORS while reading SPACE GROUP file'
         call xerror(errmsg,80,1,2)
      ELSE IF (strngs1(nword1).eq.'noprint') THEN 
         WRITE(kprint,'(5x,a31,a10)')
     &        'Space Group File read in !!--> ',cgroup
         go to 21
      ELSE
****     Write out data on symmetry
         
         WRITE(kprint,30000)
         WRITE(kprint,30010) cgroup
         WRITE(kprint,30020)
         WRITE(kprint,30030) gnmol
         DO n=1,gnmol
            WRITE(kprint,30050)
            WRITE(kprint,30040) ((cg(i,j,n),j=1,3),i=1,3)
            WRITE(kprint,30040) (tg(i,n),i=1,3)
         END DO
         WRITE(kprint,30060) 
         go to 21
      END IF
      
22    n=0
      CALL int_str (iline,strngs(1),n)
      errmsg='In SETUP_BOX: reading error in '
     &     //'SPACE GROUP file  at line:'//strngs(1)
      call xerror(errmsg,80,1,2)
      STOP
21    CONTINUE

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*
      
30000 FORMAT(//
     &'                    =========================================='/
     &'                    =             Space Group                =')
30010 FORMAT(
     &    20x, '=',10x,a20,10x,'=')
30020 FORMAT(
     &'                    =                                        ='/
     &'                    =          Symmetry Operations           ='/
     &'                    =                                        =')    
30030 FORMAT(20x,'=',18x,i4,18x,'=')
30040 FORMAT(20x,'=',5x,3f10.2,5x,'=')
30050 FORMAT(
     &'                    =                                        =')    
30060 FORMAT(
     &'                    =                                        ='/
     &'                    ==========================================')    
      RETURN
      END
