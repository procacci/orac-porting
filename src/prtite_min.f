      SUBROUTINE prtite_min

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

*-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,j
      LOGICAL hoover1
      REAL*8 ulen,cut,aux1

*==================== EXECUTABLE STATEMENTS ============================


      hoover1=.FALSE.
      IF(hoover .OR. thermos) hoover1=.TRUE.

      IF(slt_exist) THEN
          IF(.NOT. slv_exist) THEN
              WRITE(kprint,22000)
              WRITE(kprint,24000)
              WRITE(kprint,26000)
              DO 10 i=1,3
                  WRITE(kprint,28100) ntap
10            CONTINUE
              WRITE(kprint,26000)
              WRITE(kprint,24000)
          END IF
      END IF
      IF(slv_exist) THEN
          IF(.NOT. slt_exist) THEN
              WRITE(kprint,22000)
              WRITE(kprint,24000)
              WRITE(kprint,26000)
              DO 20 i=1,3
                  WRITE(kprint,28000) nmol
20            CONTINUE
              WRITE(kprint,26000)
              WRITE(kprint,24000)
          ELSE
              WRITE(kprint,22000)
              WRITE(kprint,24000)
              WRITE(kprint,26000)
              DO 30 i=1,3
                  WRITE(kprint,28300) ntap-nmol*nato_slv
                  WRITE(kprint,28200) nmol
30            CONTINUE
              WRITE(kprint,26000)
              WRITE(kprint,24000)
          END IF
      END IF

      WRITE(kprint,6000)
      cut=rspoff
      WRITE(kprint,490)
      IF(nflag(1) .EQ. 0) THEN
          IF(slvpdb) THEN
              WRITE(kprint,550)
          ELSE
              WRITE(kprint,560)
          END IF
      ELSE IF(nflag(1) .EQ. 1) THEN
          WRITE(kprint,570)
      ELSE
          WRITE(kprint,580)
      END IF
      WRITE(kprint,1000)  nrject,maxstp,maxrun
      WRITE(kprint,2000)  nprop,nsave,nplot,nascii,nconf

      WRITE(kprint,56000) cut
      WRITE(kprint,58000) time

      ulen=volume
      WRITE(kprint,54002) ntap
      WRITE(kprint,54005) ((co(i,j),j=1,3),i=1,3)
      WRITE(kprint,54003) ulen
      WRITE(kprint,3500)
      WRITE(kprint,530) nupdte,rspcut+cut
      IF(hydbnd) THEN
          aux1=DACOS(hacut)*180.0D0/pi
          WRITE(kprint,540) hrcut,aux1
      END IF
      IF(debug) THEN
         WRITE(kprint,610)
         DO i=1,nbun,6
            IF(nbun-i+1 .GE. 6) THEN
               WRITE(kprint,620) (prsymb(mend(i+j-1)),i+j-1,j=1,6)
            ELSE
               WRITE(kprint,620) (prsymb(mend(i+j-1)),i+j-1,j=1,nbun
     &              -i+1)
            END IF 
         END DO
         WRITE(kprint,630)
         WRITE(kprint,640) (beta(i),i=1,ntap)
      END IF
      WRITE(kprint,55000) ntap,lbond,lstretch,lconstr,lbend
     &              ,ltors,litor,int14p,ngrp,nbun,nprot

      WRITE(kprint,500)

*================= END OF EXECUTABLE STATEMENTS ========================

19001 FORMAT(5x,3i5,2f10.3)
14900 FORMAT(8x,3e15.5)
10000 FORMAT(//20x,19hnumber of molecules,3x,i5//
     x     20x,19hnumber of atoms    ,3x,i3)
12000 FORMAT(2(5x,i3),4(3x,f10.4))
14000 FORMAT(//10x,11hbond matrix)
15000 FORMAT(5x,4e15.6)
16000 FORMAT(//10x,8hmolecule,i6//20x,6(5x,8a1))
18000 FORMAT(12x,8a1,6(2x,f11.4))
20000 FORMAT(//10x,8hmolecule,i6,3x,13his monoatomic)
22000 FORMAT(///)
24000 FORMAT(79(1h*))
26000 FORMAT(1h*,77x,1h*)
28000 FORMAT(37h*                    minimization for,
     x  i5,37h polyatomic molecules               *)
28100 FORMAT(54h*            Minimization for an aggregate composed of,
     x i5,20h atoms             *)
28200 FORMAT(33h*         and solvated in        ,
     x  i5,41h polyatomic molecules                   *)
28300 FORMAT(54h*            Minimization for an aggregate composed of,
     x i5,20h atoms             *)
32000 FORMAT(/t10,a24,2x,a80/)
40000 FORMAT(/3x,33hpotential parameters  interaction,
     x      3x,13hlennard-jones,
     x     38h  eps(kj/mol)  sigma(ang)   charge(au),
     x     14h   cutoff(ang)  )
44000 FORMAT(10x,8a1,1h-,8a1,i2,1h-,i2,i4,5x,4f11.4)
50000 FORMAT(10x,8a1,1h-,8a1,i2,1h-,i2,i4,5x,2f11.4,
     x    11x,f11.4)
51000 FORMAT(/10x,16hmolecular mass =,f10.3,' a.m.u. '//)
52000 FORMAT(/5x,
     x' atom    type       mass          x (ang)      y (ang',
     x')      z (ang)'//)
56000 FORMAT(/
     x'          potential truncation distance =',f7.4,' ang')
58000 FORMAT(/10x,11htime-step =,1pe11.4,4h fs.)
60000 FORMAT(/10x,23hg(r) table resolution =,f7.4,4h ang)
62000 FORMAT(/10x,23hs(k) table resolution =,f7.4,8h ang**-1)
64000 FORMAT(/10x,17htail corrections /
     x   t15,8henergy =,-3pf14.4,7h kj/mol/
     x   t15,8hvirial =,f14.4,7h kj/mol/
     x   t15,8hewald  =,f14.4,7h kj/mol/
     x   t15,8hewald  =,f14.4,7h kj/mol)
500       FORMAT(//80('=')//)
490   FORMAT(
     x'==============  P A R A M E T E R S    F O R   T H E',
     x'    R U N  ================='//)
1000  FORMAT(10x,' Reject =',i6,1x,'Steps',5x,'Run =',i6,1x,
     x'Steps'/10x,' Maxrun =',i6,1x,'Steps'/)
2000  FORMAT(10x,' Save Subaverages every ',i6,1x,'Steps'/
     x10x,' Save Restart file every',i6,1x,'Steps'/
     x10x,' Save PLOT file every   ',i6,1x,'Steps'/
     x10x,' Save PDB  file every   ',i6,1x,'Steps'/
     x10x,' Save DUMP file every   ',i6,1x,'Steps')
1005  FORMAT(10x,' Reject =',i6,1x,'H-Steps',5x,'Run =',i6,1x,
     x'H-Steps'/10x,' Maxrun =',i6,1x,'H-Steps'/)
2005  FORMAT(10x,' Save Subaverages every ',i6,1x,'M-Steps'/
     x10x,' Save Restart file every',i6,1x,'H-Steps'/
     x10x,' Save PLOT file every   ',i6,1x,'M-Steps'/
     x10x,' Save PDB  file every   ',i6,1x,'M-Steps'/
     x10x,' Save DUMP file every   ',i6,1x,'M-Steps')
550   FORMAT(10x,' The Minimization will be run starting from a',
     x' PDB file '/)
560   FORMAT(10x,' The Minimization will be run from scratch.',
     x/20x,' Initial atomic coordinates will be generated. ')
570   FORMAT(10x,' The Minimization will be continued from a',
     x' RESTART file.'/30x,'Averages will be zeroed. ')
580   FORMAT(10x,' The Minimization will be run starting from a',
     x' RESTART file.'/25x,'Averages will NOT be zeroed. ')
3000  FORMAT(/
     x'-------------------------  Non Bonded Potential --------------',
     x'------------------'/)
510   FORMAT(11x,'EWALD       alpha  = ',f10.2,3x,' rkcut  = ',f10.2/
     &       11x,'            rspon  = ',f10.2,2x,' rspoff = ',f10.2/)
520   FORMAT(11x,
     x'NOEWALD     rspon  = ',f10.2,2x,' rspoff = ',f10.2/)
525   FORMAT(11x,
     &     'r-RESPA: N0-time = ',f8.3,' fs ',/,11x,
     &     '         N1-time = ',f8.3,' fs ',/,11x, 
     &     '          M-time = ',f8.3,' fs ',/,11x, 
     &     '          L-time = ',f8.3,' fs ',/,11x, 
     &     '          H-time = ',f8.3,' fs '/)
526   FORMAT(11x,
     &     'r-RESPA:  M-cut  = ',f8.3,' AA ', '--' ,2f6.2,/,11x,
     &     '          L-cut  = ',f8.3,' AA ', '--' ,2f6.2/,11x, 
     &     '          H-cut  = ',f8.3,' AA ', '--' ,2f6.2/)
515   FORMAT(11x,
     &     'PME    :  alpha  = ',f8.3,' AA^-1     grid = ',3i5/
     &11x,'          order  = ',i6/)
527   FORMAT(11x,
     &     'r-RESPA: shl-pme = ',6x,a1)
528   FORMAT(11x,
     &     'r-RESPA: shl-Ewa = ',8x,a1)
590   FORMAT(11x,
     x'H-BOND      hrson  = ',f10.2,2x,'hrsoff = ',f10.2/
     x'                       hanon  = ',f10.2,2x,'hanoff = ',
     xf10.2,1x,'nhskip = ',i2/)
3500  FORMAT(/
     x'----------------------------- Neighbor List -----------------',
     x'-------------------'/)
530   FORMAT(11x,'UPDATE every ',i5,1x,'Steps'//
     x11x,'Dispersion+Electrostatics  rspcut  = ',f8.2/)
540   FORMAT(11x,'Hydrogen bond              hrcut   = ',f8.2,1x,
     x' hacut  = ',f8.2/)
4000  FORMAT(//
     x'>>>>>>>>>>>>>>>>>>>>>>>>>>>>     SOLUTE    <<<<<<<<<<<<<<',
     x'<<<<<<<<<<<<<<<<<<<<<<<'//)
5000  FORMAT(//
     x'>>>>>>>>>>>>>>>>>>>>>>>>>>>>     SOLVENT   <<<<<<<<<<<<<<',
     x'<<<<<<<<<<<<<<<<<<<<<<<'//)
6000  FORMAT(/)
600   FORMAT(
     x'           Number of units   = ',i5,2x,' Number of groups  = ',
     xi5/)
610   FORMAT(/
     x'                    >>  Sequence of solute segments <<'//)
620   FORMAT(4x,6(a5,'(',i4,') '))
630   FORMAT(//
     x'                       >> Sequence of atom labels <<'//)
640   FORMAT(8(a5))

54000 FORMAT(/
     x'-------------------------- Statistical Ensemble -',
     x'-------------------------------')
54002 FORMAT(/
     x'           Number of atoms   =',i6)
54003 FORMAT(/
     x'           Volume            =',f15.3,'  Ang**3')
54004 FORMAT(/
     x'           Set Temperature   =',f15.3,'  Kelvin')
54005 FORMAT(/
     &     '                             CO-matrix'
     &     / /23x,3f10.5/23x,3f10.5/23x,3f10.5)
54011 FORMAT(/
     x'                   ** N V T **',
     x'                                     ')
54012 FORMAT(/
     x'           Number of atoms   =',i6)
54013 FORMAT(/
     x'           Volume            =',f15.3,'  Ang**3')
54014 FORMAT(/
     x'           Temperature       =',f15.3,'  Kelvin')
54015 FORMAT(/
     x     '                    CO-matrix'//23x,3f10.5/23x
     x     ,3f10.5/23x,3f10.5)
54016 FORMAT(/
     &     '           Thermostat mass   =',f15.3
     &     ,' a.m.u. A**2')
54021 FORMAT(/
     x'                   ** N P H **',
     x'                                     ')
54028 FORMAT(/
     x'              ** isotropic N P H **')
54022 FORMAT(/
     x'           Number of atoms   =',i6)
54023 FORMAT(/
     x'           Volume            =',f15.3,'  Ang**3')
54024 FORMAT(/
     x'           Temperature       =',f15.3,'  Kelvin')
54025 FORMAT(/
     x     '                    CO-matrix'//23x,3f10.5/23x
     x     ,3f10.5/23x,3f10.5)
54026 FORMAT(/
     x'           Barostat mass     =',f15.3,'  u.m.a.')
54027 FORMAT(/
     x'           External pressure =',f15.3,'  MPa   ')
54031 FORMAT(/
     x'                   ** N P T **',
     x'                                     ')
54038 FORMAT(/
     x'              ** isotropic N P T **')
54032 FORMAT(/
     x'           Number of atoms   =',i6)
54033 FORMAT(/
     x'           Volume            =',f15.3,'  Ang**3')
54034 FORMAT(/
     x'           Temperature       =',f15.3,'  Kelvin')
54035 FORMAT(/
     x     '                    CO-matrix'//23x,3f10.5/23x
     x     ,3f10.5/23x,3f10.5)
54036 FORMAT(/
     x'           Barostat mass     =',f15.3,'  u.m.a.')
54037 FORMAT(/
     x'           External pressure =',f15.3,'  MPa   ')
54039 FORMAT(/
     x'           Thermostat mass   =',f15.3
     &     ,' a.m.u. A**2')
54041 FORMAT(/
     x'                 ** Berendsen **')
54042 FORMAT(/
     x'           Number of atoms   =',i6)
54043 FORMAT(/
     x'           Volume            =',f15.3,'  Ang**3')
54044 FORMAT(/
     x'           Temperature       =',f15.3,'  Kelvin')
54045 FORMAT(/
     x     '                    CO-matrix'//23x,3f10.5/23x
     x     ,3f10.5/23x,3f10.5)
54046 FORMAT(/
     x'           P. coupling       =',f15.3,'  fs    ')
54047 FORMAT(/
     x'           External pressure =',f15.3,'  MPa   ')
54049 FORMAT(/
     x'           T. coupling       =',f15.3,'  fs    ')
54080 FORMAT(/
     &' 1 Thermostat Mass (cofm)    =',f15.3
     &     ,' a.m.u. A**2'/
     &' 2 Thermostat Mass (solute)  =',f15.3
     &     ,' a.m.u. A**2'/
     &' 3 Thermostat Mass (solvent) =',f15.3
     &     ,' a.m.u. A**2'/)
54117 FORMAT(/
     &'       No. of Hoover chains  =',i6)
55000 FORMAT(//
     &     '     *****************************************************',
     &     '*******************'/
     &     '     *                        System  TOPOLOGY  List      ',
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

      RETURN
      END
