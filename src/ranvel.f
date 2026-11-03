      SUBROUTINE ranvel(temp,xmp,ntap,vpx,vpy,vpz,xp,yp,zp,linit
     &     ,massinfty)

************************************************************************
* RANVEL is the initialization/scaling routine for velocities          *
* when using r-RESPA integration algorithm.                            *
* This subroutine generates the (t=0) initial velocities and scale     * 
* velocities at the selected temperature. Inizilaization is done       *
* according to a gaussian distribution: for solvent molecules,         *
* due care is taken that the starting configuration for velocities     *
* satisfies equipartition theorem.                                     *
* This routine also scale velocities during thermalization. For        *
* solvent molecules scaling is done separately on transl.,and rot+vib  *
* degrees of freedom to accelerate equipartition and equilibration     *
*                                                                      *
*     LINIT   :  Does randomization if .TRUE.                     (I)  *
*     CO,OC   :  Transform the coordinates to orthogonal          (I)  *
*                crystallographic frame, and viceversa                 *
*                >> real*8  CO(3m3) OC(3,3)  <<                        *
*     XP      :  Solute coordinates, packed by molecule.          (I)  *
*     YP         >> real*8  XP0(*), ... <<                             *
*     ZP                                                               *
*     VPX     :  Solute velocities , packed by molecule.          (O)  *
*     VPY        >> real*8  XP0(*), ... <<                             *
*     VPZ                                                              *
*     XMP     :  List of solute atomic masses, packed by molecule.(I)  *
*                >> real*8  MASSB(*) <<                                *
*     NTAP    :  Number of atomic sites on the solute.            (I)  *
*                                                                      *
*     TEMP    :  Total temperature of the system                  (O)  *
************************************************************************

*======================= DECLARATIONS ==================================

      use unit
      use parst

      IMPLICIT none

*----------------------- ARGUMENTS -------------------------------------

      logical   linit
      real*8       temp,xmp(*),vpx(*),vpy(*),vpz(*),xp(*),yp(*),zp(*)
     &     ,co(3,3),oc(3,3),massinfty
      integer      ntap

*----- LOCAL VARIABLES  AND WORK ARRAYS  -------------------------------

      INTEGER   i,l,imax,nf
      real*8    vt(3),vbx,vby,vbz,scale
      real*8    ekt,bgg,vmax,v2,ek,xmtot,xtemp
*==================== EXECUTABLE STATEMENTS ============================

c---  ekt is the selected kinetic energy (temperature temp)
c---  per degrees of freedom per mole

      ekt=0.5*gascon*temp
      vmax=0.0D0
      i=0
      if(.not.linit) then
         do i=1,ntap
            v2=xmp(i)*(vpx(i)**2+vpy(i)**2+vpz(i)**2)
            if(v2.gt.vmax) then
               vmax=v2
               imax=i
            end if   
         end do   
         if (vmax*efact/gascon.gt.20.d0*temp) then 
            LINIT=.true.
            write(6,10002) imax,vmax*efact/gascon
10002       format('####WARNING:atom ',i5,' is "hot"; tmax =',f10.1)  
         end if   
      end if
      if(linit) CALL gauss(temp,xmp,ntap,1,vpx,vpy,vpz)

c---     remove net linear momentum  

      vbx = 0.d0
      vby = 0.d0
      vbz = 0.d0
      xmtot=0.d0
      do l = 1, ntap
        if(xmp(l).gt.massinfty) THEN 
          vpx(l)=0.d0
          vpy(l)=0.d0
          vpz(l)=0.d0
        END IF
        vbx = vbx + vpx(l)*xmp(l)
        vby = vby + vpy(l)*xmp(l)
        vbz = vbz + vpz(l)*xmp(l)
        xmtot=xmp(l)+xmtot
      end do
      bgg=1.d0/(xmtot)
      vbx = vbx*bgg
      vby = vby*bgg
      vbz = vbz*bgg
      
      do l=1,ntap
         if(xmp(l).lt.massinfty) THEN
           vpx(l) =  vpx(l)-vbx
           vpy(l) =  vpy(l)-vby
           vpz(l) =  vpz(l)-vbz
         end if
      end do   
      
      ek=0.d0
      do i=1,ntap
         ek=ek+xmp(i)*(vpx(i)**2+vpy(i)**2+vpz(i)**2)
      end do   

c======= scale velocities =============================================

      ek=0.5*ek*efact
      if(ek.gt.0.d0) THEN  
        scale=dsqrt(ekt*float(3*ntap)/ek)
        do l=1,ntap
          vpx(l)=vpx(l)*scale
          vpy(l)=vpy(l)*scale
          vpz(l)=vpz(l)*scale
        end do   
      endif
      
      ek=0.d0
      do i=1,ntap
         ek=ek+xmp(i)*(vpx(i)**2+vpy(i)**2+vpz(i)**2)
      end do
      nf=3*ntap
      xtemp=ek*efact/(gascon*nf)
      
      return
      end
