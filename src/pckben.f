      SUBROUTINE pckben(nbun,mend,natm,beta,cben,nben,o1,sben,sbens,
     x                  iret,errmsg)

************************************************************************
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

      use parst

      IMPLICIT none

!==== DECLARATIONS: ===================================================*

!---- ARGUMENTS: ------------------------------------------------------*

      INTEGER o1,perm,nbun,sbens,iret
      INTEGER mend(*),natm(*),nben(*),sben(3,*)
      CHARACTER*7 beta(*),cben(3,o1,*)
      CHARACTER*80 errmsg

!---- LOCAL VARIABLES: ------------------------------------------------*

      INTEGER ma
      PARAMETER (ma=m1)
      INTEGER i, j, k, n, l, m, unit, unitm, unitp, ia, nr, ncount
      INTEGER nato,ofset,ofsetm,ofsetp,natop,natom
      INTEGER bendx(3,ma)
      INTEGER na,nb,nc,nd,kc(2)
      COMMON /rag1/ bendx
      CHARACTER*7  a(3),blank
      LOGICAL ok
      DATA blank /'       '/
      DATA kc /2,1/

!==== EXECUTABLE STATEMENTS: ==========================================*

      ofset=0
      ofsetm=0
      ofsetp=natm(mend(1))
      ncount=1
      n=0
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

          DO j=1,nben(unit)
              a(1)=cben(1,j,unit)
              a(2)=cben(2,j,unit)
              a(3)=cben(3,j,unit)
              IF(i .NE. 1 .AND. i .NE. nbun) THEN
		  ok = .TRUE.
                  DO l=1,3
                      IF( a(l)(1:1) .EQ. '-') THEN
                          DO k=1,natom
                          IF(a(l)(2:6) .EQ. beta(k+ofsetm)(1:6)) THEN
                                 bendx(l,ncount)=k+ofsetm
				 GOTO 200
                          END IF
                          END DO
                      ELSE IF(a(l)(1:1) .EQ. '+') THEN
                          DO k=1,natop
                          IF(a(l)(2:6) .EQ. beta(k+ofsetp)(1:6)) THEN
                                 bendx(l,ncount)=k+ofsetp
				 GOTO 200
                          END IF
                          END DO
                      ELSE
                          DO k=1,nato
                          IF(a(l)(1:6) .EQ. beta(k+ofset)(1:6)) THEN
                                 bendx(l,ncount)=k+ofset
				 GOTO 200
                          END IF
                          END DO
                      END IF
                      ok = .FALSE.
200                   CONTINUE
                  END DO 
		  IF(.NOT.ok) THEN
                      WRITE(6,*) i
                      errmsg='In PCKBEN: bending not found. Abort!'
                      iret=1
                      RETURN
                  END IF

              ELSE IF(i .EQ. 1 ) THEN
		  ok=.TRUE.
                  DO l=1,3
                      IF(a(l)(1:1) .EQ. '-') THEN
                          GOTO 100 

                      ELSE IF(a(l)(1:1) .EQ. '+') THEN
                          DO k=1,natop
                          IF(a(l)(2:6) .EQ. beta(k+ofsetp)(1:6)) THEN
                                  bendx(l,ncount)=k+ofsetp
				  GOTO 300
                          END IF
                          END DO
                      ELSE
                          DO k=1,nato
                          IF(a(l)(1:6) .EQ. beta(k+ofset)(1:6)) THEN
                                  bendx(l,ncount)=k+ofset
				  GOTO 300
                          END IF
                          END DO
                      END IF
                      ok = .FALSE.
300                   CONTINUE
                  END DO 
		  IF(.NOT.ok) THEN
                      WRITE(6,*) i
                      errmsg='In PCKBEN: bending not found. Abort!'
                      iret=1
                      RETURN
                  END IF

              ELSE IF(i .EQ. nbun) THEN
		  ok=.TRUE.
                  DO l=1,3
                      IF(a(l)(1:1) .EQ. '-') THEN
                          DO k=1,natom
                          IF(a(l)(2:6) .EQ. beta(k+ofsetm)(1:6)) THEN
                                  bendx(l,ncount)=k+ofsetm
				  GOTO 400
                          END IF

                          END DO
                      ELSE IF(a(l)(1:1) .EQ. '+') THEN
                          GOTO 100
                      ELSE
                          DO k=1,nato
                          IF(a(l)(1:6) .EQ. beta(k+ofset)(1:6)) THEN
                                  bendx(l,ncount)=k+ofset
				  GOTO 400
                          END IF
                          END DO
                      END IF
                      ok = .FALSE.
400                   CONTINUE
                  END DO 

		  IF(.NOT.ok) THEN
                      WRITE(6,*) i
                      errmsg='In PCKBEN: bending not found. Abort!'
                      iret=1
                      RETURN
                  END IF
              END IF
	      ncount=ncount+1
              IF(ncount .GT. ma) THEN
                  errmsg=' In PCKBEN, physical dimensions of the'//
     x                     ' work arrays are insufficient. Abort!'
                  iret=1
                  RETURN
              END IF

100           CONTINUE
          END DO

          ofset=ofset+natm(unit)
          IF(unitp .EQ. 0) THEN
             ofsetp=0
          ELSE
             ofsetp=ofsetp+natm(unitp)
          END IF
          IF(unitm .EQ. 0) THEN 
              ofsetm=0
          ELSE 
              ofsetm=ofsetm+natm(unitm)
          END IF

      END DO
      sbens=ncount-1
      DO i=1,sbens
          DO j=1,3
              sben(j,i) = bendx(j,i)
          END DO 
      END DO 
	  
                    
!---- JUMP BACK TO CALLING ROUTINE: -----------------------------------*

      RETURN
      END
