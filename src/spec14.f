      SUBROUTINE spec14(int,intp,label,type,labnb,labnbp,nato,lst14)

************************************************************************
*                                                                      *
*                                                                      *
*     SPEC14 will return an integer list of interaction types          *
*     associated to each of the third neighbour interactions.          *
*                                                                      *
*     INT     :  List of 1-4 (or third neighbour) interactions.   (I)  *
*                >> integer INT(2,*) <<                                *
*     INTP    :  Number of 1-4 interactions.                      (I)  *
*     LABEL   :  List of label types for all the solute atoms.    (I)  *
*                >> charcater*7 LABEL(NATO) <<                         *
*     TYPE    :  List of types for the solute atoms.              (I)  *
*                >> integer TYPE(NATO) <<                              *
*     LABNB   :  List of labels types for the model potential.    (I)  *
*                >> character*7 LABNB(*) <<                            *
*     LABNBP  :  Number of types in the potetial model.           (I)  *
*     NATO    :  Number of solute atoms.                          (I)  *
*     LST14   :  List of interaction types for 1-4 interactions.  (O)  *
*                >> integer*4 LST14(*) <<                              *
*                                                                      *
*---- Last update 05/19/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNAL CHK14.                                                  *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER intp,nato,labnbp,int(2,*),type(nato),lst14(*)
      CHARACTER*7 label(nato),labnb(*)

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,i1,i2,j1,j2,ij

*==================== EXECUTABLE STATEMENTS ============================

      DO 10 i=1,intp
          i1=int(1,i)
          i2=int(2,i)
          j1=type(i1)
          j2=type(i2)
          CALL chk14(labnb,labnbp,label(i1),label(i2),j1,j2)
          i1=MAX0(j1,j2)
          i2=MIN0(j1,j2)
          ij=i1*(i1-1)/2+i2
          lst14(i)=ij
10    CONTINUE

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
