      SUBROUTINE pckfrg(nbun,mend,natm,beta,stpr,stsc,nstpr,nstsc,
     x           nstrg,ma,mb,bsitp,asitp,nbsitp,nasitp,nrigg,
     x           mc,md,iret,errmsg)

************************************************************************
*                                                                      *
*     PCKFRG makes a list of protein fragments from the topology       *
*     file.                                                            *
*                                                                      *
*======================================================================*
*                                                                      *
*     NBUN     : Number of residues                              (I)   *
*     MEND     : Sequence of protein residues                    (I)   *
*                >>INTEGER MEND(*) <<                                  *
*     NATM     : Number of atoms per residue                     (I)   *
*                >>INTEGER NATM(*) <<                                  *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

!==== DECLARATIONS: ===================================================*

      IMPLICIT none

!---- ARGUMENTS: ------------------------------------------------------*

      INTEGER nbun,nrigg,ma,mb,mc,md,iret
      INTEGER mend(*),natm(*),nstpr(ma,*),nstsc(ma,*),nstrg(*)
      INTEGER bsitp(mc,*),asitp(mc,*),nbsitp(*),nasitp(*)
      CHARACTER*7 beta(*),stpr(4,ma,*),stsc(mb,ma,*)
      CHARACTER*80 errmsg

!---- LOCAL VARIABLES: ------------------------------------------------*

      INTEGER i, j, k, n, l, m, unit, unitm, unitp, ia, nr, ncount
      INTEGER nato,ofset,ofsetm,ofsetp,natop,natom
      CHARACTER*7 a(4),b(40)
      LOGICAL ok

!==== EXECUTABLE STATEMENTS: ==========================================*

      ofset=0
      ofsetm=0
      ofsetp=natm(mend(1))
      ncount=1
      n=0
      iret=0
      DO i=1,nbun
          unit=mend(i)
          nato=natm(unit)
          IF(i .EQ. 1) THEN
              unitp=mend(i+1)
              unitm=0
              natop=natm(unitp)
              natom=0
          ELSE IF(i .EQ. nbun) THEN
              unitp=0
              unitm=mend(i-1)
              natop=0
              natom=natm(unitm)
          ELSE
              unitp=mend(i+1)
              unitm=mend(i-1)
              natop=natm(unitp)
              natom=natm(unitm)
          END IF

          DO j=1,nstrg(unit)
              DO m=1,nstpr(j,unit)
                  a(m)=stpr(m,j,unit)
              END DO
              DO m=1,nstsc(j,unit)
                  b(m)=stsc(m,j,unit)
              END DO
              nasitp(ncount)=nstpr(j,unit)
              nbsitp(ncount)=nstsc(j,unit)

!----------

              IF(i .NE. 1 .AND. i .NE. nbun) THEN
		  ok = .TRUE.
                  DO l=1,nstpr(j,unit)
                      IF( a(l)(1:1) .EQ. '-') THEN
                          DO k=1,natom
                          IF(a(l)(2:6) .EQ. beta(k+ofsetm)(1:6)) THEN
                                 asitp(l,ncount)=k+ofsetm
				 GOTO 200
                          END IF
                          END DO
                      ELSE IF(a(l)(1:1) .EQ. '+') THEN
                          DO k=1,natop
                          IF(a(l)(2:6) .EQ. beta(k+ofsetp)(1:6)) THEN
                                 asitp(l,ncount)=k+ofsetp
				 GOTO 200
                          END IF
                          END DO
                      ELSE
                          DO k=1,nato
                          IF(a(l)(1:6) .EQ. beta(k+ofset)(1:6)) THEN
                                 asitp(l,ncount)=k+ofset
				 GOTO 200
                          END IF
                          END DO
                      END IF
                      ok = .FALSE.
200                   CONTINUE
                  END DO 

                  DO l=1,nstsc(j,unit)
                      IF( b(l)(1:1) .EQ. '-') THEN
                          DO k=1,natom
                          IF(b(l)(2:6) .EQ. beta(k+ofsetm)(1:6)) THEN
                                 bsitp(l,ncount)=k+ofsetm
				 GOTO 210
                          END IF
                          END DO
                      ELSE IF(b(l)(1:1) .EQ. '+') THEN
                          DO k=1,natop
                          IF(b(l)(2:6) .EQ. beta(k+ofsetp)(1:6)) THEN
                                 bsitp(l,ncount)=k+ofsetp
				 GOTO 210
                          END IF
                          END DO
                      ELSE
                          DO k=1,nato
                          IF(b(l)(1:6) .EQ. beta(k+ofset)(1:6)) THEN
                                 bsitp(l,ncount)=k+ofset
				 GOTO 210
                          END IF
                          END DO
                      END IF
                      ok = .FALSE.
210                   CONTINUE
                  END DO 
                  IF(.NOT.ok) THEN
                      WRITE(6,*) i
                      errmsg='In PCKFRG: Fragment not found. Abort!'
                      iret=1
                      RETURN
                  END IF

              ELSE IF(i .EQ. 1 ) THEN
		  ok=.TRUE.
                  DO l=1,nstpr(j,unit)
                      IF(a(l)(1:1) .EQ. '-') THEN
                          GOTO 100 

                      ELSE IF(a(l)(1:1) .EQ. '+') THEN
                          DO k=1,natop
                          IF(a(l)(2:6) .EQ. beta(k+ofsetp)(1:6)) THEN
                                 asitp(l,ncount)=k+ofsetp
                                 GOTO 300
                          END IF
                          END DO
                      ELSE
                          DO k=1,nato
                          IF(a(l)(1:6) .EQ. beta(k+ofset)(1:6)) THEN
                                 asitp(l,ncount)=k+ofset
                                 GOTO 300
                          END IF
                          END DO
                      END IF
                      ok = .FALSE.
300                   CONTINUE
                  END DO
                  DO l=1,nstsc(j,unit)
                      IF(b(l)(1:1) .EQ. '-') THEN
                          GOTO 100

                      ELSE IF(b(l)(1:1) .EQ. '+') THEN
                          DO k=1,natop
                          IF(b(l)(2:6) .EQ. beta(k+ofsetp)(1:6)) THEN
                                 bsitp(l,ncount)=k+ofsetp
                                 GOTO 310
                          END IF
                          END DO
                      ELSE
                          DO k=1,nato
                          IF(b(l)(1:6) .EQ. beta(k+ofset)(1:6)) THEN
                                 bsitp(l,ncount)=k+ofset
                                 GOTO 310
                          END IF
                          END DO
                      END IF
                      ok = .FALSE.
310                   CONTINUE
                  END DO


		  IF(.NOT.ok) THEN
                      WRITE(6,*) i
                      errmsg='In PCKFRG: torsion not found. Abort!'
                      iret=1
                      RETURN
                  END IF

              ELSE IF(i .EQ. nbun) THEN
		  ok=.TRUE.
                  DO l=1,nstpr(j,unit)
                      IF( a(l)(1:1) .EQ. '-') THEN
                          DO k=1,natom
                          IF(a(l)(2:6) .EQ. beta(k+ofsetm)(1:6)) THEN
                                 asitp(l,ncount)=k+ofsetm
                                 GOTO 400
                          END IF
                          END DO

                      ELSE IF(a(l)(1:1) .EQ. '+') THEN
                          GOTO 100
                      ELSE
                          DO k=1,nato
                          IF(a(l)(1:6) .EQ. beta(k+ofset)(1:6)) THEN
                                 asitp(l,ncount)=k+ofset
                                 GOTO 400
                          END IF
                          END DO

                      END IF
                      ok = .FALSE.
400                   CONTINUE
                  END DO 
                  DO l=1,nstsc(j,unit)
                      IF( b(l)(1:1) .EQ. '-') THEN
                          DO k=1,natom
                          IF(b(l)(2:6) .EQ. beta(k+ofsetm)(1:6)) THEN
                                 bsitp(l,ncount)=k+ofsetm
                                 GOTO 410
                          END IF
                          END DO

                      ELSE IF(b(l)(1:1) .EQ. '+') THEN
                          GOTO 100
                      ELSE
                          DO k=1,nato
                          IF(b(l)(1:6) .EQ. beta(k+ofset)(1:6)) THEN
                                 bsitp(l,ncount)=k+ofset
                                 GOTO 410
                          END IF
                          END DO

                      END IF
                      ok = .FALSE.
410                   CONTINUE
                  END DO 

		  IF(.NOT.ok) THEN
                      WRITE(6,*) i
                      errmsg='In PCKFRG: torsion not found. Abort!'
                      iret=1
                      RETURN
                  END IF
              END IF
	      ncount=ncount+1
              IF(ncount .GT. md) THEN
                  errmsg=' In PCKFRG, physical dimensions of the'//
     x                     ' work arrays are insufficient. Abort!'
                  iret=1
                  RETURN
              END IF

100           CONTINUE
          END DO

          ofset=ofset+natm(unit)
          IF(unitp.EQ.0) THEN 
            ofsetp=0
          ELSE
            ofsetp=ofsetp+natm(unitp)
          END IF
          IF(unitm.EQ.0) THEN 
              ofsetm=0
          ELSE 
              ofsetm=ofsetm+natm(unitm)
          END IF

      END DO
      
      nrigg=ncount-1

!---- JUMP BACK TO CALLING ROUTINE: -----------------------------------*

      RETURN
      END
