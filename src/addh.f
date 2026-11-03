      SUBROUTINE addh(ia,nbl,nhb,blist,xh,yh,zh,type,x0,y0,z0,iret,
     x                errmsg)

************************************************************************
*                                                                      *
*     Compute the position of the hydrogen atoms bonded to a given     *
*     atom of the solute.                                              *
*                                                                      *
*     IA      :  Index of the central atom in the solute list.    (I)  *
*     NBL     :  Number of atoms bonded to IA other than hydrogens(I)  *
*     NHB     :  Number of hydrogen bonded to IA.                 (I)  *
*     BLIST   :  List of atoms bonded to IA.                      (I)  *
*                >> integer BLIST(*) <<                                *
*     XH      :  Cartesian coordinates of the solute atoms.       (I)  *
*     YH         >> real*8 XH(*), .... <<                              *
*     ZH                                                               *
*     TYPE    :  Type label of the atom IA.                       (I)  *
*     X0      :  Cartesian coordinates of the hydrogen atoms      (O)  *
*     Y0         >> real*8 X0(*), ... <<                               *
*     Z0                                                               *
*                                                                      *
*---- Last update 06/20/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNAL ADD1H, ADD2H, ADD3H                                     *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER ia,nbl,nhb,iret,blist(*)
      REAL*8  xh(*),yh(*),zh(*),x0(*),y0(*),z0(*)
      CHARACTER*7 type
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      REAL*8 xa(3),ya(3),za(3),ax(3),pax(3)
      REAL*8 length,oh,n2h,n3h,sh,c3h,c2h
      CHARACTER*1 oo,nn,ss,cc,char
      INTEGER i,j,ntot
      DATA oh,n2h,n3h,sh,c2h,c3h/0.960d0,1.01d0,1.01d0,1.336d0,
     x     1.08d0,1.09d0/
      DATA oo,nn,ss,cc/'o','n','s','c'/

*==================== EXECUTABLE STATEMENTS ============================

      ntot=nbl+nhb
      IF(ntot.GT.4) THEN
          iret=1
          errmsg=' IN ADDH: Cannot handle atoms with more than 4 bonds*
     x ABORT. '
          RETURN
      END IF

*=======================================================================
*---- Compute group axis -----------------------------------------------
*=======================================================================

      DO 10 i=1,nbl
          j=blist(i)
          xa(i)=xh(ia)-xh(j)
          ya(i)=yh(ia)-yh(j)
          za(i)=zh(ia)-zh(j)
10    CONTINUE
      DO 20 i=1,3
          ax(i)=0.0D0
20    CONTINUE
      DO 30 j=1,nbl
          ax(1)=ax(1)+xa(j)
          ax(2)=ax(2)+ya(j)
          ax(3)=ax(3)+za(j)
30    CONTINUE

*=======================================================================
*---- Calculate bond length from atom type -----------------------------
*=======================================================================

      char=type(1:1)
      IF(char.EQ.oo) THEN
          length=oh
      ELSE IF(char.EQ.nn) THEN
          IF(ntot.EQ.3) THEN
              length=n2h
          ELSE IF(ntot.EQ.4) THEN
              length=n3h
          END IF
      ELSE IF(char.EQ.ss) THEN
          length=sh
      ELSE IF(char.EQ.cc) THEN
          IF(ntot.EQ.3) THEN
              length=c2h
          ELSE IF(ntot.EQ.4) THEN
              length=c3h
          END IF
      END IF

*=======================================================================
*---- Compute position of the hydrogens w.r.t. the atom IA -------------
*=======================================================================

      IF(nhb.EQ.1) THEN
          CALL add1h(ax,length,x0,y0,z0)
      ELSE IF(nhb.EQ.3) THEN
          CALL add3h(ax,length,x0,y0,z0)
      ELSE IF(nhb.EQ.2) THEN
          IF(nbl.EQ.2) THEN
              pax(1)=ya(1)*za(2)-za(1)*ya(2)
              pax(2)=za(1)*xa(2)-xa(1)*za(2)
              pax(3)=xa(1)*ya(2)-ya(1)*xa(2)
              CALL add2h(ax,pax,1,length,x0,y0,z0)
          ELSE IF(nbl .EQ. 0) THEN
              CALL addh2o(x0,y0,z0)
          ELSE
              CALL add2h(ax,pax,0,length,x0,y0,z0)
          END IF
      END IF

*=======================================================================
*---- Compute cartesian coordinates of the hydrogens -------------------
*=======================================================================

      DO 40 i=1,nhb
          x0(i)=x0(i)+xh(ia)
          y0(i)=y0(i)+yh(ia)
          z0(i)=z0(i)+zh(ia)
40    CONTINUE

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
