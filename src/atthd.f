      SUBROUTINE atthd(list,hlist,rlist,nlist,xh,yh,zh,beta,nato,
     x           wkx,wky,wkz,iret,errmsg)

************************************************************************
*                                                                      *
*     ATTHD will provide on return a set of hydrogen coordinates       *
*     appended to the list of the solute coordinates.                  *
*                                                                      *
*                                                                      *
*     CO      :  Rotation matrix from simulation box frame to     (I)  *
*                orthogonalized frame.                                 *
*                >> real*8 CO(3,3) <<                                  *
*     OC      :  Inverse of CO.                                   (I)  *
*                >> real*8 OC(3,3) <<                                  *
*     LIST    :  Lists of atoms bonded to hydrogens.              (I)  *
*                >> integer LIST(3,*) <<                               *
*                LIST(1,I) = I-th atom bonded to hydrogens.            *
*                LIST(2,I) = Number of hydrogens bonded to the         *
*                            I-th atom.                                *
*                LIST(3,I) = Number of non-hydrogens bonded to         *
*                            the I-th atom.                            *
*     HLIST   :  List of the hydrogens.                           (I)  *
*                >> integer HLIST(3,*) <<                              *
*     RLIST   :  List of hydrogens second neighbours.             (I)  *
*                >> integer RLIST(3,*) <<                              *
*     NLIST   :  Number of atoms bonded to hydrogens.             (I)  *
*     XH      :  Solute coordinates.                            (I/O)  *
*     YH         >> real*8 XH(*), ... <<                               *
*     ZH                                                               *
*     BETA    :  List of labels for the macromolecule atoms.      (I)  *
*                >> character*7 BETA(*) <<                             *
*     NATO    :  Number of solute atoms.                          (I)  *
*     WKX     :  Work arrays.                                     (I)  *
*     WKY        >> real*8 WKX(*), ... <<                              *
*     WKZ                                                              *
*     IRET    :  Return code.                                     (O)  *
*     ERRMSG  :  Error message.                                   (O)  *
*                                                                      *
*---- Last update 06/21/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNAL ADDH                                                    *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER nlist,iret,nato,list(3,*),hlist(3,*),rlist(3,*)
      CHARACTER*7 beta(*)
      REAL*8  xh(*),yh(*),zh(*),wkx(*),wky(*),wkz(*)
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,m,n,blist(3),nbl,nhb
      REAL*8  xb(3),yb(3),zb(3),xc(3),yc(3),zc(3)

*==================== EXECUTABLE STATEMENTS ============================

      DO 10 n=1,nato
         wkx(n)=xh(n)
         wky(n)=yh(n)
         wkz(n)=zh(n)
10    CONTINUE
      DO 20 i=1,nlist
          j=list(1,i)
          nhb=list(2,i)
          nbl=list(3,i)
          DO 30 n=1,nbl
              blist(n)=rlist(n,i)
30        CONTINUE
          CALL addh(j,nbl,nhb,blist,wkx,wky,wkz,beta(j),xb,yb,zb,iret,
     x                errmsg)
          IF(iret.EQ.1) RETURN
          DO 40 n=1,nhb
              xc(n)=xb(n)
              yc(n)=yb(n)
              zc(n)=zb(n)
40        CONTINUE
          DO 50 n=1,nhb
              m=hlist(n,i)
              xh(m)=xc(n)
              yh(m)=yc(n)
              zh(m)=zc(n)
50        CONTINUE
20    CONTINUE

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
