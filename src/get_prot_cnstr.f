      SUBROUTINE get_prot_cnstr(nato,protl,nprot,cnstp,cnstpp,cnst_protl
     &     ,cnst_protp,tot_prot1,tot_prot2,mask,work,ma)

************************************************************************
*   Time-stamp: <97/12/08 09:14:36 marchi>                             *
*                                                                      *
*                                                                      *
*                                                                      *
*======================================================================*
*                                                                      *
*              Author:  Massimo Marchi                                 *
*              CEA/Centre d'Etudes Saclay, FRANCE                      *
*                                                                      *
*              - Thu Feb 27 1997 -                                     *
*                                                                      *
************************************************************************

*---- This subroutine is part of the program ORAC ----*


*======================== DECLARATIONS ================================*

      IMPLICIT none

*----------------------------- ARGUMENTS ------------------------------*

      INTEGER ma,cnst_protp,tot_prot1,tot_prot2
      INTEGER protl(*),nprot,nato,cnstp(2,*),cnstpp,cnst_protl(*)
     &     ,work(ma,*)
      LOGICAL mask(*)

*------------------------- LOCAL VARIABLES ----------------------------*

      INTEGER i,k,j,ka,lc,atomp,n,m,la,lb

*----------------------- EXECUTABLE STATEMENTS ------------------------*

      DO i=1,nato
         work(i,1)=0
      END DO

      DO k=1,cnstpp
         la=cnstp(1,k)
         lb=cnstp(2,k)
         work(la,1)=work(la,1)+1
         work(la,1+work(la,1))=k

         work(lb,1)=work(lb,1)+1
         work(lb,1+work(lb,1))=k
      END DO

      DO k=1,cnstpp
         mask(k)=.TRUE.
      END DO

      tot_prot1=0
      tot_prot2=0
      
      cnst_protp=0
      DO i=1,nprot
         atomp=protl(1+tot_prot1)
         n=0
         IF(atomp .GT. 1) THEN
            cnst_protp=cnst_protp+1
            DO j=1,atomp
               lc=protl(1+tot_prot1+j)
               m=work(lc,1)
               DO ka=1,m
                  k=work(lc,1+ka)
                  IF(mask(k)) THEN
                     n=n+1
                     mask(k)=.FALSE.
                     cnst_protl(1+tot_prot2+n)=k
                  END IF
               END DO
            END DO
            cnst_protl(1+tot_prot2)=n
            tot_prot2=tot_prot2+n+1
         END IF
         tot_prot1=tot_prot1+atomp+1
      END DO

*----------------- END OF EXECUTABLE STATEMENTS -----------------------*

      RETURN
      END
