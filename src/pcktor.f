      SUBROUTINE pcktor(nbun,mend,natm,beta,ctor,ntor,o1,stor,
     x           stors,iret,errmsg)

************************************************************************
*                                                                      *
*     PCKTOR makes a list of torsions from the topology file.          *
*     The starting inputs are a list of atom labels, the list          *
*     of residues for the protein and for each residue its             *
*     torsions.                                                        *
*                                                                      *
*======================================================================*
*                                                                      *
*     NBUN     : Number of residues                              (I)   *
*     MEND     : Sequence of protein residues                    (I)   *
*                >>INTEGER MEND(*) <<                                  *
*     NATM     : Number of atoms per residue                     (I)   *
*                >>INTEGER NATM(*) <<                                  *
*     RES      :                                                       *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

!==== DECLARATIONS: ===================================================*

      use parst

      IMPLICIT none

!---- ARGUMENTS: ------------------------------------------------------*

      INTEGER o1,perm,nbun,stors,iret
      INTEGER mend(*),natm(*),ntor(*),stor(4,*)
      CHARACTER*7 beta(*),ctor(4,o1,*)
      CHARACTER*80 errmsg

!---- LOCAL VARIABLES: ------------------------------------------------*

      INTEGER ma
      PARAMETER (ma=m3)
      INTEGER i, j, k, n, l, m, unit, unitm, unitp, ia, nr, ncount
      INTEGER nato,ofset,ofsetm,ofsetp,natop,natom
      INTEGER iperm(6,4),ipt(2,4),ipi(6,4)
      INTEGER torsx(4,ma)
      INTEGER na,nb,nc,nd,kc(2)
      COMMON /rag1/ torsx
      CHARACTER*7  a(4),blank
      LOGICAL ok
      DATA blank /'       '/
      DATA kc /2,1/

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

          DO j=1,ntor(unit)
              a(1)=ctor(1,j,unit)
              a(2)=ctor(2,j,unit)
              a(3)=ctor(3,j,unit)
              a(4)=ctor(4,j,unit)
              IF(i .NE. 1 .AND. i .NE. nbun) THEN
		  ok = .TRUE.
                  DO l=1,4
                      IF( a(l)(1:1) .EQ. '-') THEN
                          DO k=1,natom
                          IF(a(l)(2:6) .EQ. beta(k+ofsetm)(1:6)) THEN
                                 torsx(l,ncount)=k+ofsetm
				 GOTO 200
                          END IF
                          END DO
                      ELSE IF(a(l)(1:1) .EQ. '+') THEN
                          DO k=1,natop
                          IF(a(l)(2:6) .EQ. beta(k+ofsetp)(1:6)) THEN
                                 torsx(l,ncount)=k+ofsetp
				 GOTO 200
                          END IF
                          END DO
                      ELSE
                          DO k=1,nato
                          IF(a(l)(1:6) .EQ. beta(k+ofset)(1:6)) THEN
                                 torsx(l,ncount)=k+ofset
				 GOTO 200
                          END IF
                          END DO
                      END IF
                      ok = .FALSE.
200                   CONTINUE
                  END DO 
		  IF(.NOT.ok) THEN
                      WRITE(6,*) i,a(1),a(2),a(3),a(4)
                      errmsg='In PCKTOR: torsion not found. Abort!'
                      iret=1
                      RETURN
                  END IF

              ELSE IF(i .EQ. 1 ) THEN
		  ok=.TRUE.
                  DO l=1,4
                      IF(a(l)(1:1) .EQ. '-') THEN
                          GOTO 100 

                      ELSE IF(a(l)(1:1) .EQ. '+') THEN
                          DO k=1,natop
                          IF(a(l)(2:6) .EQ. beta(k+ofsetp)(1:6)) THEN
                                  torsx(l,ncount)=k+ofsetp
				  GOTO 300
                          END IF
                          END DO
                      ELSE
                          DO k=1,nato
                          IF(a(l)(1:6) .EQ. beta(k+ofset)(1:6)) THEN
                                  torsx(l,ncount)=k+ofset
				  GOTO 300
                          END IF
                          END DO
                      END IF
                      ok = .FALSE.
300                   CONTINUE
                  END DO 
		  IF(.NOT.ok) THEN
                    WRITE(6,*) i,a(1),a(2),a(3),a(4)
                      errmsg='In PCKTOR: torsion not found. Abort!'
                      iret=1
                      RETURN
                  END IF

              ELSE IF(i .EQ. nbun) THEN
		  ok=.TRUE.
                  DO l=1,4
                      IF(a(l)(1:1) .EQ. '-') THEN
                          DO k=1,natom
                          IF(a(l)(2:6) .EQ. beta(k+ofsetm)(1:6)) THEN
                                  torsx(l,ncount)=k+ofsetm
				  GOTO 400
                          END IF

                          END DO
                      ELSE IF(a(l)(1:1) .EQ. '+') THEN
                          GOTO 100
                      ELSE
                          DO k=1,nato
                          IF(a(l)(1:6) .EQ. beta(k+ofset)(1:6)) THEN
                                  torsx(l,ncount)=k+ofset
				  GOTO 400
                          END IF
                          END DO
                      END IF
                      ok = .FALSE.
400                   CONTINUE
                  END DO 

		  IF(.NOT.ok) THEN
                      errmsg='In PCKTOR: torsion not found. Abort!'
                      iret=1
                      RETURN
                  END IF
              END IF
	      ncount=ncount+1
              IF(ncount .GT. ma) THEN
                  errmsg=' In PCKTOR, physical dimensions of the'//
     x                     ' work arrays are insufficient. Abort!'
                  iret=1
                  RETURN
              END IF

100           CONTINUE
              WRITE(6,*) i,a(1),a(2),a(3),a(4),OK
          END DO

          ofset=ofset+natm(unit)
          IF(unitp .EQ. 0) THEN
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
      stors=ncount-1
      DO i=1,stors
          DO j=1,4
              stor(j,i) = torsx(j,i)
          END DO 
      END DO 
	  
                    
!---- JUMP BACK TO CALLING ROUTINE: -----------------------------------*

      RETURN
      END
