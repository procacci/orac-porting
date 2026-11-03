      SUBROUTINE inicmp(ss_index,x,y,z,xcm,ycm,zcm,mass,nprot,protl)

************************************************************************
*                                                                      *
*     Initialize the center of mass coordinates of the solute          *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      REAL*8  mass(*),x(*),y(*),z(*)
      REAL*8  xcm(*),ycm(*),zcm(*)
      REAL*8  mtot
      INTEGER i,j,count,nprot,m,i1,ss_index(*)
      INTEGER protl(*)

*==================== EXECUTABLE STATEMENTS ============================


      count=0
      DO j=1,nprot

C    Computes centres of mass in reduced coordinates in the unshrinked frame
         mtot = 0.0D0
         xcm(j) = 0.0D0
         ycm(j) = 0.0D0
         zcm(j) = 0.0D0
         m=protl(1+count)
         DO i=1,m
            i1=protl(1+count+i)
            mtot = mtot + mass(i1)
            xcm(j) = xcm(j) + x(i1)*mass(i1)
            ycm(j) = ycm(j) + y(i1)*mass(i1)
            zcm(j) = zcm(j) + z(i1)*mass(i1)
         END DO
         xcm(j) = xcm(j)/mtot
         ycm(j) = ycm(j)/mtot
         zcm(j) = zcm(j)/mtot
         count=count+m+1
      END DO 

*================= END OF EXECUTABLE STATEMENTS ========================

      RETURN
      END
