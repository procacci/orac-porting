      SUBROUTINE getchr(ffield,fexcl,nfatm,nexatm,nfield,n1,n2,n3,
     x           iret,errmsg)

************************************************************************
*                                                                      *
*     Read chromophore information: atoms forming the chromophore      *
*     excluded fields.                                                 *
*                                                                      *
*======================================================================*
*                                                                      *
*     UNIT    :  List of labels for all the residues.                  *
*                >> character*8 UNIT(N1) <<                            *
*     NUNIT   :  Number of labels in the list of residues.             *
*     N1      :  Physical dimension of UNIT.                           *
*================== RETURN CODES ======================================*
*                                                                      *
*     IRET    :  Return code.                                          *
*     ERRMSG  :  Error message.                                        *
*                                                                      *
*----- Last update  09/04/89 ------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS   Free format package, RDUNIT, CSORT, ATOMS,           *
*                 BOND, BACKBN, HBOND, NEAR0.                          *
*                                                                      *
*                                                                      *
************************************************************************

      use unit
      
      IMPLICIT none
      INTEGER n1,n2,n3,ffield(n1,*),fexcl(n3,*),nfatm(*),nexatm(*),
     x        nfield
      CHARACTER*80 errmsg,line
      INTEGER iret

      INTEGER i,n,naux,nword,ncount
      CHARACTER*80 strngs(40)
      CHARACTER*8 fmt
      CHARACTER*1 sep(2),comm(2)
      REAL*8  dummy,dummy1,fior
      COMMON /card / jrec,nrec,jerr,nerr,nline,istrt
      INTEGER jrec,nrec,jerr,nerr,nline,istrt(80)
      DATA sep/' ',','/comm/'#','!'/

      nfield=0
      line(79:80)='  '

!=======================================================================
!==== Scan the records =================================================
!=======================================================================

100   READ(knlist,'(a78)',END=600) line(1:78)
      IF(jerr.EQ.1) CALL wrenc(kprint,line)
      CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
      IF(strngs(1) .EQ. 'END') GOTO 500

!=======================================================================
!==== Case: 1) record contains keyword "chromo" ========================
!=======================================================================

      IF(strngs(1) .EQ. 'chromo') THEN
          nfield=nfield+1
          nexatm(nfield)=0
200       READ(knlist,'(a78)',END=600) line(1:78)
          IF(jerr.EQ.1) CALL wrenc(kprint,line)
          IF(line(1:1) .EQ. '#') GOTO 200 
          CALL parse(line,sep,2,comm,strngs,40,nword,iret,errmsg)
250       CONTINUE
          IF(strngs(1) .NE. 'end') THEN
              IF(strngs(1) .EQ. 'group') THEN
300               READ(knlist,'(a78)',END=600) line(1:78)
                  IF(jerr.EQ.1) CALL wrenc(kprint,line)
                  IF(line(1:1) .EQ. '#') GOTO 300 
                  CALL parse(line,sep,2,comm,strngs,40,nword,
     x                       iret,errmsg)

!=======================================================================
!==== Checks if labels correspond to a keyword =========================
!=======================================================================

                  IF(strngs(1) .EQ. 'end') GOTO 250
                  IF(strngs(1) .EQ. 'exclude') GOTO 250
                  IF(strngs(1) .EQ. 'group') GOTO 250

!=======================================================================
!==== Finally reads the atom list ======================================
!=======================================================================

                  DO i=1,nword
                      IF(strngs(i) .NE. '-') THEN
                          nfatm(nfield)=nfatm(nfield)+1
                          CALL fndfmt(1,strngs(i),fmt)
                          READ(strngs(i),fmt) 
     x                    ffield(nfatm(nfield),nfield)
                      ELSE IF(i .NE. nword .AND. i .NE. 1) THEN
                          naux=ffield(nfatm(nfield),nfield)
                          CALL fndfmt(1,strngs(i+1),fmt)
                          READ(strngs(i+1),fmt) ncount
                          DO n=naux+1,ncount-1
                              nfatm(nfield)=nfatm(nfield)+1
                              ffield(nfatm(nfield),nfield)=n
                          END DO
                      ELSE
                          iret=1
                          errmsg=
     x                'In getchr: Symbol ''-'' cannot be '
     x           //'at the end or beginning of a line. Abort.'
                          RETURN
                      END IF
                  END DO
                  GOTO 300

!=======================================================================
!==== Case: 1) record contains keyword "exclude" =======================
!=======================================================================


              ELSE IF(strngs(1) .EQ. 'exclude') THEN
400               READ(knlist,'(a78)',END=600) line(1:78)
                  IF(jerr.EQ.1) CALL wrenc(kprint,line)
                  IF(line(1:1) .EQ. '#') GOTO 400 
                  CALL parse(line,sep,2,comm,strngs,40,nword,
     x                       iret,errmsg)

!=======================================================================
!==== Checks if labels correspond to a keyword =========================
!=======================================================================

                  IF(strngs(1) .EQ. 'end') GOTO 250
                  IF(strngs(1) .EQ. 'exclude') GOTO 250
                  IF(strngs(1) .EQ. 'group') GOTO 250

!=======================================================================
!==== Finally reads the list of excluded atoms =========================
!=======================================================================

                  DO i=1,nword
                      IF(strngs(i) .NE. '-') THEN
                          nexatm(nfield)=nexatm(nfield)+1
                          CALL fndfmt(1,strngs(i),fmt)
                          READ(strngs(i),fmt)
     x                    fexcl(nexatm(nfield),nfield)
                      ELSE IF(i .NE. nword .AND. i .NE. 1) THEN
                          naux=fexcl(nexatm(nfield),nfield)
                          CALL fndfmt(1,strngs(i+1),fmt)
                          READ(strngs(i+1),fmt) ncount
                          DO n=naux+1,ncount-1
                              nexatm(nfield)=nexatm(nfield)+1
                              fexcl(nexatm(nfield),nfield)=n
                          END DO
                      ELSE
                          iret=1
                          errmsg=
     x                'In getchr: Symbol ''-'' cannot be '
     x           //'at the end or beginning of a line. Abort.'
                          RETURN
                      END IF
                  END DO
                  GOTO 400

              ELSE
                  iret=1
                  errmsg=
     x'In getchr: Subcommand '// strngs(1)(1:11)//' not available'//
     x' . Abort'
                  RETURN
              END IF
          ELSE
              GOTO 100
          END IF
      END IF

500   CONTINUE
      iret=0
      RETURN

600   CONTINUE
      iret=1
      errmsg=
     x'In getchr: End of file found. Abort!'
      RETURN
      END
