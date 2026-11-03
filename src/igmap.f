      SUBROUTINE igmap(ngrp,grppt,ingrpp,ingrp,ma,mapnl,errmsg,iret)

************************************************************************
*   Time-stamp: <97/09/05 16:30:44 marchi>                             *
*   This routine computes the non-bonded interactions occurring        *
*   among atoms of the same group.                                     *
*                                                                      *
*   NGRP:     Number of group in the solute                            *
*   GRPPT:    First and last atom of the group                         *
*              >GRPPT(2,*)<                                            *
*   INGRPP:   Number of nb interaction among atoms of the same         *
*             group.                                                   *
*   INGRP:    List of the interactions.                                *
*             >INGRP(2,*)<                                             *
*   MA:       Physical dimension of INGRP.                             *
*   MAPNL:    Map of excluded interactions.                            *
*             >MAPNL(*)<                                               *
*   ERRMSG:   Error message.                                           *
*   IRET:     Error flag.                                              *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Tue Sep 19 1995 -                                     *
*                                                                      *
************************************************************************
                        
*---- This subroutine is part of the program ORAC ----*
                        
                        
      IMPLICIT none     
      INTEGER ngrp,ingrpp,grppt(2,*),ingrp(2,*)
      INTEGER mapnl(*),ma,iret
                        
      INTEGER i,j,i1,j1,noff,n,m,m1,map
      CHARACTER*80 errmsg
                        
      n=0               
      map=0             
      iret=0            
      DO i=1,ngrp
          DO i1=grppt(1,i),grppt(2,i)
              m=mapnl(1+n)
              DO j1=i1+1,grppt(2,i)
                 DO m1=1,m
                    IF(j1 .EQ. mapnl(1+n+m1)) THEN
                       GOTO 100
                    END IF
                 END DO
                 map=map+1
                 IF(map .GT. ma) THEN
                    iret=1
                    errmsg=' In IGMAP: Physical dimension of ingrp'
     x                   //' exceeded. Abort'
                 END IF
                 ingrp(1,map)=i1
                 ingrp(2,map)=j1
100              CONTINUE
              END DO
              noff=mapnl(1+n)+1
              n=n+noff
          END DO
      END DO
      ingrpp=map
      RETURN
      END
