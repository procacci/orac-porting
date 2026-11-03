      REAL*8 FUNCTION meta_pot(rc)
! compute the total biasing potential in the point rc

C======================= DECLARATIONS ==================================
      
      use unit
      use meta
      
      IMPLICIT none   
      
      real*8 rc(metadim_nrc)
      
C-------------------- LOCAL VARIABLES ----------------------------------

      INTEGER i,irc,irc2,t,zer0
      REAL*8 rc0, dr0, dra, lucy
      REAL*8 dr(metadim_nrc), wi(metadim_nrc)

C==================== EXECUTABLE STATEMENTS ============================
      
      meta_pot = 0.d0
      
      DO irc = 1,meta_nrc
        wi(irc) = 1.d0/meta_width(irc)
      ENDDO
      
      IF(meta_n.GT.0) THEN
        zer0 = 0
        
        DO t = 1,meta_n      

          IF(meta_nrc.GT.1) THEN
            IF(zer0.GT.0) THEN  ! if the last hill was far in CV "zer0", check if it is still far in that CV
              rc0 = meta_list(t,zer0)
              dr0 = rc(zer0) - rc0
              if((meta_rctype(zer0).eq.'tors').or.
     &             (meta_rctype(zer0).eq.'bend'))then 
                if(dr0.gt.180.0) dr0 = dr0 - 360.0
                if(dr0.le.-180.0) dr0 = dr0 + 360.0
              endif
              dra = DABS(dr0)
              IF(dra.GT.meta_width(zer0)) GOTO 1000 ! goto the next hill
            ENDIF
          ENDIF
          
          DO irc = 1,meta_nrc     
            rc0 = meta_list(t,irc)
            dr(irc) = rc(irc) - rc0
            dr0 = dr(irc)
            if((meta_rctype(irc).eq.'tors').or.(meta_rctype(irc).eq
     &           .'bend'))then 
              if(dr0.gt.180.0) dr0 = dr0 - 360.0
              if(dr0.le.-180.0) dr0 = dr0 + 360.0
            endif
            dr(irc) = dr0
            dra = DABS(dr0)
            IF(dra.LE.meta_width(irc)) THEN
              zer0 = 0
            ELSE                ! do not calculate forces, goto the next hill
              zer0 = irc
              GOTO 1000        
            ENDIF
          ENDDO            
          
          lucy = 1.d0
          DO irc = 1,meta_nrc     
            dra = dabs(dr(irc))
            lucy = lucy * (1.d0 + 2.d0 * wi(irc) * dra) * (1.d0 ! this is a unnormalized monodimensional LF with height=1
     &           -wi(irc) *dra)**2  
          ENDDO
          
          meta_pot = meta_pot + lucy
          
1000      CONTINUE

        ENDDO           

        meta_pot = meta_height * meta_pot
        
      ENDIF
      
      RETURN
      END
      
