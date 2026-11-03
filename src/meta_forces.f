      SUBROUTINE meta_forces(xp0,yp0,zp0,fpx,fpy,fpz,co,oc,ntap)
************************************************************************
*
*     Compute the forces from the history dependent potential.
*     The outer loop is on hills deposed during the run.
*     The inner loop is on the collective variables. 
*     If the distance of any of the RCs from the
*     center of a hill is greater than the chosen width along that
*     direction (the hill is "far"), then
*     the contribution of that hill to the forces along the RCs is zero.
*     Therefore, the algorithm can skip to the next hill in time without
*     calculating forces. Moreover, if the last hill was "far" in a
*     given RC, the algorithm check first the distance along that RC. If
*     the latter is still greater than the width, again it skips to the
*     next hill.
*
*     The functional form of the hills (in the 1D case) is given by:
*     f(x)  =  h * (1 + n * |x|/w) * (1 - |x|/w)**n       (|x| <= w)
*     f(x)  =  0                                          (|x| >  w)
*     centered at x=0, where w is the range and the "height"=f(0)=h.
*     This function (Lucy's function, LC) is the simplest possible
*     function satisfying n+2 desiderable conditions:
*     i) it is normalizable
*     ii) it has a finite range w
*     iii) it has a maximum at the origin
*     iv) to n+2) it has n-1 continuous dervative everywhere.
*     Since we need one continuous derivative only, n=2 and 
*     f(x)  =  h * (1 + 2 * |x|/w) * (1 - |x|/w)**2       (|x| <= w)
*     f(x)  =  0                                          (|x| >  w)  
*
*     The derivative is then simply
*     f'(x) = (6 * h / w**3) * x * (|x| - w)
*     
*     For the bidimensional case, one LF is given by
*     f(x,y) = h * (1 + 2 * |x|/w) * (1 - |x|/w)**2 * 
*                  (1 + 2 * |y|/w) * (1 - |y|/w)**2      
*                  (|x| <=w_{x}, |y| <= w_{y})
*     and therefore 
*     f(x,y) = f(x) * f(y) / h
*
*     The derivative with respect to x is then given by
*     f'(x,y) = f'(x) * f(y) / h
*     
*     In N dimensions:
*     f({\bf x})  = (1/h)**(N-1) * \prod_i^N (f_{i}(x))
*     f'({\bf x}) = (1/h)**(N-1) * f'(x) * \prod_{i=/j} f_i(x) 
************************************************************************
      
      
C======================= DECLARATIONS ==================================

      use unit
      use meta
      
      IMPLICIT none     
      
C----------------------- ARGUMENTS -------------------------------------
      
      INTEGER ntap
      REAL*8 xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*)
      REAL*8 ddx(ntap),ddy(ntap),ddz(ntap)
      REAL*8 co(3,3),oc(3,3)

C-------------------- LOCAL VARIABLES ----------------------------------     
      INTEGER i,irc,irc2,t,zer0
      REAL*8 dr0, dra, ff
      REAL*8 dr(metadim_nrc), wi(metadim_nrc), wi3(metadim_nrc),
     &     f(metadim_nrc),lucy(metadim_nrc), drdx(ntap,metadim_nrc),
     &     drdy(ntap,metadim_nrc),drdz(ntap,metadim_nrc)
      
C==================== EXECUTABLE STATEMENTS ============================
     
      
      IF (meta_run) THEN
        
! for each reaction coordinate, compute its value and the value of its
! derivatives with respect to the cartesian coordinates of the atoms
! involved
        DO irc = 1,meta_nrc
          IF (meta_rctype(irc) .EQ. 'bond') THEN
            CALL meta_bond(meta_atoms(1,irc),meta_atoms(2,irc),ntap,xp0 
     &           ,yp0,zp0,co,oc,meta_rc(irc),ddx,ddy,ddz) 
          ENDIF
          IF (meta_rctype(irc) .EQ. 'bend') THEN
            CALL meta_bend(meta_atoms(1,irc),meta_atoms(2,irc)
     &           ,meta_atoms(3,irc),ntap,xp0,yp0,zp0,co,oc,meta_rc(irc)
     &           ,ddx,ddy,ddz) 
          ENDIF
          IF (meta_rctype(irc) .EQ. 'tors') THEN
            CALL meta_tors(meta_atoms(1,irc),meta_atoms(2,irc)
     &           ,meta_atoms(3,irc),meta_atoms(4,irc),ntap,xp0,yp0,zp0
     &           ,co,oc,meta_rc(irc),ddx,ddy,ddz) 
          ENDIF
! store the derivatives in the matrices drdx,drdy,drdz
          DO i = 1,ntap
            drdx(i,irc) = ddx(i)
            drdy(i,irc) = ddy(i)
            drdz(i,irc) = ddz(i)
          ENDDO
          wi(irc) = 1.d0/meta_width(irc)
          wi3(irc) = wi(irc)**3
          f(irc) = 0.d0
        ENDDO
        
        IF (.NOT. meta_analysis) THEN ! during an analysis run, do not compute forces 
          
          IF (meta_n.GT.0) THEN
            zer0 = 0
          
! loop on the hills
            DO t = 1,meta_n           
              
! when the dimensionality of the RCs space is greater than 1, check
! first in the direction of the last far hill
              IF (meta_nrc.GT.1) THEN
                IF (zer0.GT.0) THEN ! if the last hill was far in RC "zer0", check if it is still far in that RC
                  dr0 = meta_rc(zer0) - meta_list(t,zer0)
                  IF ((meta_rctype(zer0) .EQ. 'tors') .OR.
     &                 (meta_rctype(zer0) .EQ. 'bend')) THEN
                    IF (dr0 .GT. 180.0) dr0 = dr0 - 360.0
                    IF (dr0 .LE. -180.0) dr0 = dr0 + 360.0
                  ENDIF
                  dra = DABS(dr0)
                  IF (dra .GT. meta_width(zer0)) GOTO 1000 ! do not calculate forces, goto the next hill
                ENDIF
              ENDIF
              
! compute the value of a monodimensional LF for each reaction coordinate
! goto the next hill if zero
              DO irc = 1,meta_nrc
                dr(irc) = meta_rc(irc) - meta_list(t,irc)
                dr0 = dr(irc)
                IF ((meta_rctype(irc) .EQ.
     &               'tors') .OR. (meta_rctype(irc) .EQ. 'bend')) THEN 
                  IF (dr0 .GT. 180.0) dr0 = dr0 - 360.0
                  IF (dr0 .LE. -180.0) dr0 = dr0 + 360.0
                ENDIF
                dra = DABS(dr0)
                IF (dra .LE. meta_width(irc)) THEN
                  zer0 = 0
                  lucy(irc) = (1.d0 + 2.d0 * wi(irc) * dra) * (1.d0 - ! this is a unnormalized monodimensional LF with height=1
     &                 wi(irc) *dra)**2  
                ELSE            ! do not calculate forces, goto the next hill
                  zer0 = irc
                  GOTO 1000        
                ENDIF
              ENDDO            
              
! compute the metadynamics forces for each reaction coordinate and add
! them to the total forces
              DO irc = 1,meta_nrc     
                dr0 = dr(irc)
                IF ((meta_rctype(irc) .EQ. 'tors') .OR.
     &               (meta_rctype(irc) .EQ. 'bend')) THEN
                  IF (dr0 .GT. 180.0) dr0 = dr0 - 360.0
                  IF (dr0 .LE. -180.0) dr0 = dr0 + 360.0
                ENDIF
                dra = DABS(dr0)
                ff = - dr0 * (dra - meta_width(irc))
                DO irc2 = 1,meta_nrc
                  IF (irc2 .NE. irc) ff = ff * lucy(irc2)
                ENDDO
! ff is the contribution of the t-th hill to the metadynamics force
! along i-th reaction coordinate
                f(irc) = f(irc) + ff
              ENDDO
              
1000          CONTINUE
              
            ENDDO           
            
            DO irc = 1,meta_nrc
              f(irc) = f(irc) * 6.d0 * meta_height * wi3(irc)
            ENDDO
            
            DO i = 1,ntap
              DO irc = 1,meta_nrc
                fpx(i) = fpx(i) + f(irc) * drdx(i,irc)
                fpy(i) = fpy(i) + f(irc) * drdy(i,irc) ! sum the components due to the i-th collective variable to the corrections of the cartesian forces
                fpz(i) = fpz(i) + f(irc) * drdz(i,irc)
              ENDDO
            ENDDO
            
          ENDIF
          
        ENDIF
        
      ENDIF
      
      
      RETURN
      END
      
      
      
      
