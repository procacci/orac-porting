      SUBROUTINE serchd(x0,y0,z0,nato,connct,n1,list,hlist,rlist,
     x                  nlist,iret,errmsg)

************************************************************************
*                                                                      *
*     SERCHD will search for those hydrogens whose coordinates have    *
*     not been read in and create 3 lists containing the atoms         *
*     bonded to those hydrogens, the hydrogens themselves and          *
*     their second neghbours.                                          *
*                                                                      *
*     X0      :  Solute coordinates in simulation box frame.      (I)  *
*     Y0         >> real*8 X0(*), ... <<                               *
*     Z0                                                               *
*     NATO    :  Number of solute atoms.                          (I)  *
*     CONNCT  :  Connection table.                                (I)  *
*                >> integer CONNCT(N1,10) <<                           *
*     N1      :  Physical row dimension of CONNCT.                (I)  *
*     LIST    :  Lists of atoms bonded to hydrogens.              (O)  *
*                >> integer LIST(3,*) <<                               *
*                LIST(1,I) = I-th atom bonded to hydrogens.            *
*                LIST(2,I) = Number of hydrogens bonded to the         *
*                            I-th atom.                                *
*                LIST(3,I) = Number of non-hydrogens bonded to         *
*                            the I-th atom.                            *
*     HLIST   :  List of the hydrogens.                           (O)  *
*                >> integer HLIST(3,*) <<                              *
*     RLIST   :  List of hydrogens second neighbours.             (O)  *
*                >> integer RLIST(3,*) <<                              *
*     NLIST   :  Number of atoms bonded to hydrogens.             (O)  *
*     IRET    :  Return code.                                     (O)  *
*     ERRMSG  :  Error message.                                   (O)  *
*                                                                      *
*---- Last update 06/21/89 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                      *
*     EXTERNAL none                                                    *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      INTEGER n1,nlist,nato,iret,connct(n1,*),list(3,*),hlist(3,*),
     x        rlist(3,*)
      REAL*8  x0(*),y0(*),z0(*)
      CHARACTER*80 errmsg

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j,n,m,map,coord
      REAL*8 xg,yg,zg,tg

*==================== EXECUTABLE STATEMENTS ============================

      DO 10 i=1,nato
          DO 20 j=1,3
              list(j,i)=0
              hlist(j,i)=0
              rlist(j,i)=0
20        CONTINUE
10    CONTINUE
      map=0
      DO 30 i=1,nato
          xg=x0(i)
          yg=y0(i)
          zg=z0(i)
          tg=xg+yg+zg
          IF(tg.GT.1.0D+5) THEN
              coord=connct(i,2)
              DO 40 j=1,map
                  IF(coord.EQ.list(1,j)) THEN
                      list(2,j)=list(2,j)+1
                      hlist(list(2,j),j)=i
                      GOTO 200
                  END IF
40            CONTINUE
              map=map+1
              list(1,map)=coord
              list(2,map)=list(2,map)+1
              hlist(list(2,map),map)=i
200           CONTINUE
          END IF
30    CONTINUE
      nlist=map
      DO 50 i=1,nlist
          j=list(1,i)
          coord=connct(j,1)
          list(3,i)=coord-list(2,i)
          map=0
          DO 60 n=1,coord
              DO 70 m=1,list(2,i)
                  IF(connct(j,n+1).EQ.hlist(m,i)) THEN
                      GOTO 100
                  END IF
70            CONTINUE
              map=map+1
              rlist(map,i)=connct(j,n+1)
100           CONTINUE
60        CONTINUE
          IF(map.NE.list(3,i)) THEN
              iret=1
              errmsg=' IN SERCHD : Number of 2nd neighbours does not mat
     xch. ABORT * '
	      WRITE(6,1001) j,i,list(2,i),list(3,i),map
1001          FORMAT(" ATOM", I5,2x,"(",I5,"-th atom with H)  has ",I3
     &             " h neighbors and ",I3, 
     &             " non-h neighbors while expecting",I3)     
c	      WRITE(6,1002) list(1,i-1),list(2,i-1),
c     x                                list(3,i-1)
              RETURN
          END IF
50    CONTINUE

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
