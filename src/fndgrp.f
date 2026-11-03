      SUBROUTINE fndgrp(nprot,protl,atomp)

************************************************************************
*                                                                      *
*          => nprot,protl
*         <=  atomp                                       *
*                                                                      *
************************************************************************

           IMPLICIT none

!-------------Arguments -------------------------------------------------
           INTEGER protl(*),atomp(*),nprot

!---------------Local variables ----------------------------------------

           INTEGER count,i,m,i1,j

           count=0
           DO i=1,nprot
              m=protl(count+1)
              DO i1=1,m
                 j=protl(count+1+i1)
                 atomp(j)=i
              END DO
              count=count+m+1
           END DO

           RETURN
           END
