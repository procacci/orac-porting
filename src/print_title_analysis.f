      SUBROUTINE print_title_analysis(fstep)

************************************************************************
*                                                                      *
*     Print titles for this run. There is no argument in output.       *
*                                                                      *
*                                                                      *
*---- Last update 02/05/93 --------------------------------------------*
*                                                                      *
*     Written by Massimo Marchi CECAM Orsay 1993                       *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit
      use parst
      use cpropar

      IMPLICIT none
      REAL*8   fstep

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j
      REAL*8 ulen,volm,cut,aux1,aux2

*==================== EXECUTABLE STATEMENTS ============================

      WRITE(kprint,1000) 
      WRITE(kprint,2000) 
      WRITE(kprint,4000) 
      WRITE(kprint,3000)
      WRITE(kprint,4000)
      IF(slt_exist) THEN
         IF(.NOT. slv_exist) THEN
            DO i=1,3
               WRITE(kprint,5000) ntap
            END DO
         END IF
      END IF
      IF(slv_exist) THEN
         IF(.NOT. slt_exist) THEN
            WRITE(kprint,6000) nmol*nato_slv
         ELSE
            WRITE(kprint,7000) ntap-nmol*nato_slv,nmol*nato_slv
         END IF
      END IF
      WRITE(kprint,4000)
      WRITE(kprint,2000) 
      WRITE(kprint,1000) 

      WRITE(kprint,10000)

      WRITE(kprint,11000) fstep
      WRITE(kprint,12000)
      IF(gofr) THEN
         WRITE(kprint,18000)

         WRITE(kprint,19000) DFLOAT(gofr_ncomp)*fstep,DFLOAT(gofr_navg)
     &        *fstep,DFLOAT(gofr_nprint)*fstep
      END IF
      
      IF(anxrms) THEN
         WRITE(kprint,20000) 
      END IF 
      IF(avg_str) THEN
         WRITE(kprint,21000) DFLOAT(navg_str_xrms)*fstep
      END IF

      WRITE(kprint,22000)

      WRITE(kprint,8000) nmol,nato_slv,nprot-nmol,ntap,lbond,lstretch
     &     ,lconstr,lbend,ltors,litor,int14p,ngrp,nbun,nprot

      WRITE(kprint,9000)

*================= END OF EXECUTABLE STATEMENTS ========================

1000  FORMAT(///)
2000  FORMAT('********************************************************',
     &'***********************')
3000  FORMAT('*           Analysis of an M.D. Trajectory for a System ',
     &'Composed of           *')
4000  FORMAT('*                                                       ',
     &'                      *')
5000  FORMAT('*                           ',i6,2x,'   Solute Atoms     ',
     &'                      *')
6000  FORMAT('*                           ',i6,2x,'   Solvent Atoms    ',
     &'                      *')
7000  FORMAT('*            ',i6,2x,'   Solute and   ',i6,2x,
     &'  Solvent Atoms                  *')
8000  FORMAT(//
     &     '     *****************************************************',
     &     '*******************'/
     &     '     *                            TOPOLOGY  List          ',
     &     '                  *'/
     &     '     *                                                    ',
     &     '                  *'/
     &     '     *          ',i6,' Solvent Mol  ',i6,' Atoms Each',i6,
     &     ' Solute Mol.     *'/
     &     '     *                                                    ',
     &     '                  *'/
     &     '     *                            For the Sysyem          ',
     &     '                  *'/
     &     '     *                                                    ',
     &     '                  *'/
     &     '     *          ',i6,' Atoms        ',i6,' Bonds     ',i6,
     &     ' FLexible Bonds  *'/
     &     '     *          ',i6,' Rigid Bonds  ',i6,' Angles    ',i6,
     &     ' P-Torsions      *'/
     &     '     *          ',i6,' I-Torsions   ',i6,' 1-4 Inter.',i6,
     &     ' Atomic Groups   *'/
     &     '     *          ',i6,' Units        ',i6,' Molecules      ',
     &     '                  *'/
     &     '     *                                                    ',
     &     '                  *'/
     &     '     *****************************************************',
     &     '*******************'/)
9000  FORMAT(//80('=')////)
10000 FORMAT(/'==========================    S Y N O P S I S   =======',
     &     '========================'///)

11000 FORMAT('                Coordinates on the Trajectory File were',
     &' collected'/
     &     '                with a frequency of  ',f10.3,'  fs'/)
12000 FORMAT('                The following properties will be',
     &' computed:'/)
13000 FORMAT('                          VORONOI Volumes     ')
14000 FORMAT('                          VORONOI Areas       ')
15000 FORMAT('                          VORONOI Contact Surfaces')
16000 FORMAT('                          VORONOI Neighbors   ')
17000 FORMAT('                               Freq.  ',f10.3,'  fs'//)
18000 FORMAT('                          PAIR CORRELATIONS Solvent'/
     &       '                          PAIR CORRELATIONS Solute')
19000 FORMAT('                            Freq. calc. ',f10.3,'  fs',
     &      /'                            Avg. over   ',f10.3,'  fs',
     &      /'                            Freq. prt.  ',f10.3,'  fs'/)
20000 FORMAT('                    X-RMS instantaneous deviation  ',    
     & ' from X-Ray')
21000 FORMAT('                    X-RMS deviation of the average ',
     & ' from X-Ray'/
     &     '                            Freq.       ',f10.3,'  fs')
22000 FORMAT(//'=====================================================',
     &     '=========================='///)

      RETURN
      END
