      SUBROUTINE addtor(beta,atres,rtor,xtor,xntor,stor,stors,n1,
     x           impr,connct,n2,iret,errmsg)

************************************************************************
*                                                                      *
*                                                                      *
*     This subroutine adds extra torsions to the list of protein       *
*     torsions.                                                        *
*                                                                      *
*======================================================================*
*                                                                      *
*                                                                      *
*     BETA    : List of atomic labels for the protein/macromolecule.   *
*               >> CHARACTER*7 BETA(*) <<                              *
*     ATRES   : List of pointers to the initial and final atoms        *
*               of the residues of the protein/macromolecule.          *
*               >> INTEGER ATRES(2,*) <<                               *
*               ATRES(1,NRES) provides the address of the first        *
*               atom of the residue NRES. ATRES(2,NRES) gives          *
*               the address of the final atom.                         *
*     RTOR    : List of torsions to be added by atom nuber.            *
*               >> INTEGER RTOR(4,*) <<                                *
*     XTOR    : List of torsions to be added by atom label.            *
*               >> CHARACTER*7 XTOR(4,*) <<                            *
*     XNTOR   : Number of new extra torsions.                          *
*     STOR    : In input list of existing torsions. In output          *
*               include the old and extra torsions. By atom number.    *
*               >> INTEGER STOR(4,*) <<                                *
*     STORS   : In input, number of existing torsions. In output,      *
*               number of the old and extra torsions.                  *
*     N1      : Physical dimension of STOR.                            *
*     IMPR    : Logical flag. It is .TRUE. if the torsion is           *
*               improper.                                              *
*     CONNCT  : List of links between atoms.                           *
*               >> INTEGER CONNCT(N2,*) <<                             *
*     N2      : First dimension of CONNCT                              *
*                                                                      *
*================== RETURN CODES ======================================*
*                                                                      *
*     IRET    :  Return code.                                          *
*     ERRMSG  :  Error message.                                        *
*                                                                      *
*----- Last update  01/05/92 ------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi UC Berkeley 1992                       *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER iret,n1,n2,xntor,stors
      INTEGER atres(2,*),rtor(4,*),stor(4,*),connct(n2,*)
      CHARACTER*7 beta(*),xtor(4,*)
      CHARACTER*80 errmsg
      LOGICAL impr

*-------------------- LOCAL VARIABLES ----------------------------------

      CHARACTER*7 char1,char2,char3,char4
      INTEGER n,m,m1,m2,m3,m4,i,j,k,l,tors
      LOGICAL ok,cc

*==================== EXECUTABLE STATEMENTS ============================

      tors=stors
      iret=0
      DO n=1,xntor
          tors=tors+1
          IF(tors .GT.n1) THEN
             iret=1
             errmsg=' Error in ADDTOR: Physical dimensions of torsion'//
     x       ' list is insufficient. Abort.'
             RETURN
          END IF
          m1=rtor(1,n)
          m2=rtor(2,n)
          m3=rtor(3,n)
          m4=rtor(4,n)
          char1=xtor(1,n)
          char2=xtor(2,n)
          char3=xtor(3,n)
          char4=xtor(4,n)
          ok=.FALSE.

!=======================================================================
!----- Match atoms in the list with those of the protein ---------------
!=======================================================================

          DO m=atres(1,m1),atres(2,m1)
              IF(char1 .EQ. beta(m)) THEN
                  i=m
                  ok=.TRUE.
              END IF
          END DO
          IF(.NOT. ok) THEN
              iret=1
              errmsg=' Error in ADDTOR: First atom of the link NOT'//
     x               ' found in the protein list. Abort.'
              RETURN
          END IF
       
          ok=.FALSE.
          DO m=atres(1,m2),atres(2,m2)
              IF(char2 .EQ. beta(m)) THEN
                  j=m
                  ok=.TRUE.
              END IF
          END DO
          IF(.NOT. ok) THEN
              iret=1
              errmsg=' Error in ADDTOR: Second atom of the link NOT'//
     x               ' found in the protein list. Abort.'
              RETURN
          END IF
          ok=.FALSE.
          DO m=atres(1,m3),atres(2,m3)
              IF(char3 .EQ. beta(m)) THEN
                  k=m
                  ok=.TRUE.
              END IF
          END DO
          IF(.NOT. ok) THEN
              iret=1
              errmsg=' Error in ADDTOR: Third atom of the link NOT'//
     x               ' found in the protein list. Abort.'
              RETURN
          END IF
          ok=.FALSE.
          DO m=atres(1,m4),atres(2,m4)
              IF(char4 .EQ. beta(m)) THEN
                  l=m
                  ok=.TRUE.
              END IF
          END DO
          IF(.NOT. ok) THEN
              iret=1
              errmsg=' Error in ADDTOR: Fourth atom of the link NOT'//
     x               ' found in the protein list. Abort.'
              RETURN
          END IF

!=======================================================================
!----- Add the new torsions to the list by atoms number ----------------
!=======================================================================

          IF(impr) THEN
              IF(cc(i,j,connct,n2) .AND. cc(i,k,connct,n2) 
     x        .AND. cc(i,l,connct,n2) ) THEN
                  stor(1,tors)=i
                  stor(2,tors)=j
                  stor(3,tors)=k
                  stor(4,tors)=l
              ELSE
                  iret=1
              errmsg=' Error in ADDTOR: Not an improper torsion.'//   
     x               ' Connections should be 1-2, 1-3, 1-4.Abort.'
                  RETURN
              END IF
          ELSE
              IF(cc(i,j,connct,n2) .AND. cc(j,k,connct,n2)
     x        .AND. cc(k,l,connct,n2) ) THEN
                  stor(1,tors)=i
                  stor(2,tors)=j
                  stor(3,tors)=k
                  stor(4,tors)=l
              ELSE
                  iret=1
              errmsg=' Error in ADDTOR: Not a proper torsion.'//   
     x               ' Connections should be 1-2, 2-3, 3-4.Abort.'
                  RETURN
              END IF
          END IF

      END DO
      stors=tors

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
