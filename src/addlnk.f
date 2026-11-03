      SUBROUTINE addlnk(beta,atres,connct,n1,n2,linkn,linka,linkp,iret,
     x                  errmsg)

************************************************************************
*                                                                      *
*                                                                      *
*     This subroutine adds a link to the list of interesidue           *
*     links.                                                           *
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
*     CONNCT  : List of links between atoms.                           *
*               >> INTEGER CONNCT(N1,N2) <<                            *
*     N1      : First dimension of CONNCT                              *
*     N2      : Second dimension of CONNCT                             *
*     LINKN   : Collection of couples of residues that are link        *
*               together.                                              *
*               >> INTEGER LINKN(2,*) <<                               *
*     LINKA   : List of atom labels which are linked together          *
*               for the residue in LINKN.                              *
*               >> CHARACTER*7 LINKA(2,*) <<                           *
*     LINKP   : Number of links that are to be added to CONNCT.        *
*                                                                      *
*                                                                      *
*================== RETURN CODES ======================================*
*                                                                      *
*     IRET    :  Return code.                                          *
*     ERRMSG  :  Error message.                                        *
*                                                                      *
*----- Last update  25/10/92 ------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CECAM, Universite de Paris-Sud         *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER iret,n1,n2,linkp
      INTEGER atres(2,*),connct(n1,n2),linkn(2,*)
      CHARACTER*7 beta(*),linka(2,*)
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      CHARACTER*7 char1,char2,chari
      INTEGER n,m,m1,m2,i,j,coordi,coordj
      LOGICAL ok

*==================== EXECUTABLE STATEMENTS ============================

      DO n=1,linkp
          m1=linkn(1,n)
          m2=linkn(2,n)
          char1=linka(1,n)
          char2=linka(2,n)
          ok=.FALSE.
          DO m=atres(1,m1),atres(2,m1)
              IF(char1 .EQ. beta(m)) THEN
                  i=m
                  ok=.TRUE.
              END IF
          END DO
          IF(.NOT. ok) THEN
              iret=1
              write(chari,'(i5)') m1 
              errmsg=' Error: 1st atom ' // char1 // ' of res ' // chari
     &             //'not'//' found in the protein list. Abort.'
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
              write(chari,'(i5)') m2
              errmsg=' Error: 2nd atom '// char2 //' of res '//chari
     &             //'not'//' found in the protein list. Abort.'
              RETURN
          END IF
          coordi            =connct(i,1) + 1
         IF(coordi .GT. n2-1) THEN
             iret=1
             errmsg=' Error in ADDLNK: connection tables exceed '//
     x              'physical dimensions. Abort.'
          END IF
          connct(i,1)       =coordi
          connct(i,1+coordi)=j
          coordj            =connct(j,1) + 1
         IF(coordj .GT. n2-1) THEN
             iret=1
             errmsg=' Error in ADDLNK: connection tables exceed '//
     x              'physical dimensions. Abort.'
          END IF
          connct(j,1)       =coordj
          connct(j,1+coordj)=i
      END DO

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
