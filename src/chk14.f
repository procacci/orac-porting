      SUBROUTINE chk14(label,nlabel,alfa,alfb,i1,i2)

************************************************************************
*                                                                      *
*     CHK14 will return the two atom types of a 1-4 interaction.       *
*     The routine checks if amongst the nonbonded labels there are     *
*     two whose the first n-1 non-blank characters are identical       *
*     respectively to ALFA and ALFB and the n-th non-blank character   *
*     is equal to '*'. If for both ALFA and ALFB such label does       *
*     not exist on return i1 and i2 are zero.                          *
*                                                                      *
*     LABEL   :  List of non-bonded labels.                       (I)  *
*                >> character*7 LABEL(*) <<                            *
*     NLABEL  :  Number of non-bonded labels.                     (I)  *
*     ALFA    :  Label type of atom 1 of the actual 1-4           (I)  *
*                interaction.                                          *
*                >> character*7 ALFA <<                                *
*     ALFB    :  Label type of atom 2 of the actual 1-4           (I)  *
*                interaction.                                          *
*                >> character*7 ALFB <<                                *
*     I1      :  Atomic type of the atom 1. It is returned      (I/O)  *
*                as zero if the serch is unsuccesfull.                 *
*     I2      :  Atomic type of the atom 4. It is returned      (I/O)  *
*                as zero if the serch is unsuccesfull.                 *
*                                                                      *
*---- Last update 05/19/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNALS NONE.                                                  *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT CHARACTER*80(a-z)

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nlabel,i1,i2
      CHARACTER*7 alfa,alfb,label(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j
      CHARACTER*7 a,b,blank7
      CHARACTER*1 blank,star
      DATA blank/' '/star/'*'/
      DATA blank7/'       '/

*==================== EXECUTABLE STATEMENTS ============================

      a=blank7
      b=blank7
      j=0
      DO 30 i=1,7
          IF(alfa(i:i).EQ.blank(1:1)) THEN
              j=i
              GOTO 10
          END IF
30    CONTINUE
10    CONTINUE
      DO 40 i=1,j-1
          a(i:i)=alfa(i:i)
40    CONTINUE
      a(j:j)=star
      j=0
      DO 50 i=1,7
          IF(alfb(i:i).EQ.blank(1:1)) THEN
              j=i
              GOTO 20
          END IF
50    CONTINUE
20    CONTINUE
      DO 60 i=1,j-1
          b(i:i)=alfb(i:i)
60    CONTINUE
      b(j:j)=star
      DO 70 i=1,nlabel
          IF(label(i).EQ.a) THEN
              i1=i
          END IF
          IF(label(i).EQ.b) THEN
              i2=i
          END IF
70    CONTINUE

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
