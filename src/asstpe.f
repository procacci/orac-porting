      SUBROUTINE asstpe(label,type,mass,natom,ntype,typep,massp,
     x           bug,nbug)

************************************************************************
*                                                                      *
*     The subroutine ASSTPE assign to each atom of the protein         *
*     its type number and its mass.                                    *
*                                                                      *
*     LABEL     :  List of label type for each atom           (INPUT)  *
*                  >> character*7 LABEL(NATOM) <<                      *
*     TYPE      :  List of all the label types in the         (INPUT)  *
*                  model potential                                     *
*                  >> character*7 TYPE(NTYPE) <<                       *
*     MASS      :  List of all the mass types in the          (INPUT)  *
*                  model potential                                     *
*                  >> real*8 MASS(NTYPE) <<                            *
*     NATOM     :  Number of atoms forming the protein        (INPUT)  *
*                  >> integer NATOM <<                                 *
*     NTYPE     :  Number of atom types in the                (INPUT)  *
*                  model potential                                     *
*                  >> integer NTYPE <<                                 *
*     TYPEP     :  List of types for each atom               (OUTPUT)  *
*                  >> integer TYPEP(NATOM)                             *
*     MASSP     :  List of masses for each atom              (OUTPUT)  *
*                  >> real*8 MASSP(NATOM) <<                           *
*     IRET      :  Return code                               (OUTPUT)  *
*                  >> integer iret <<                                  *
*     ERRMSG    :  Error message                             (OUTPUT)  *
*                  >> character*80 ERRMSG <<                           *
*                                                                      *
*---- Last update 04/01/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*                                                                      *
*     ETERNALS :   NONE                                                *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT CHARACTER*80(a-z)

*----------------------- ARGUMENTS -------------------------------------

      INTEGER natom,ntype,typep(natom),iret,bug(*),nbug
      CHARACTER*7 label(natom),type(ntype)
      REAL*8 mass(ntype),massp(natom)
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER m,n
      LOGICAL err,ok

*==================== EXECUTABLE STATEMENTS ============================

      m=1
      err=.false.
      
      nbug=0
      DO m=1,natom
         DO n = 1,ntype
            IF(label(m).eq.type(n)) THEN
               typep(m)=n
c              assign mass only if unassigned from tpg
               if(massp(m).lt.1.d-30) THEN
                 massp(m)=mass(n)
               ENDIF 
               GO TO 100 
            END IF
         END DO
         nbug=nbug+1
         bug(nbug)=m 
c---     one thousend Bugs: we had enough 
         if(nbug.gt.1000) RETURN 
 100     CONTINUE 
      END DO

*==================== END EXECUTABLE STATEMENTS ========================

      RETURN
      END
