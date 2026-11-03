*=========================================================================  
      SUBROUTINE comp_work(nbondsadded,nbendsadded,nitorsadded,strbonds
     &     ,strbends,strtors,bondadded,bendadded,torsadded,eqdist
     &     ,eqdist1,eqang,eqang1,eqdied,eqdied1,force_cost,force_ang
     &     ,force_died,strvel,ben_vel,tor_vel,strtime0
     &     ,strtime1,strtime0r,strtime1r,radfact,timesteer,nplotstn0
     &     ,ninn0,kplotst0,workbo,workbe,workto,rtime,iaux)
*=========================================================================  
*    Computes and print works along reaction coordinates 
*   
*======================= DECLARATIONS ====================================

      IMPLICIT NONE

*----------------------- ARGUMENTS -------------------------------------
      
      INTEGER nbondsadded,nbendsadded,nitorsadded,nplotstn0,ninn0
     &     ,kplotst0,iaux
      REAL*8  bondadded(*),bendadded(*),torsadded(*),eqdist(*),eqdist1(
     &     *),eqang(*),eqang1(*),eqdied(*),eqdied1(*),force_cost(*)
     &     ,force_ang(*),force_died(*),strvel(*),ben_vel(*),tor_vel(*)
     &     ,workbo(*),workbe(*),workto(*)
     &     ,strtime0,strtime1,strtime0r,strtime1r,rtime,timesteer

      INTEGER strbonds(*),strbends(*),strtors(*)
      
*----------------------- LOCAL VARIABLES -------------------------------

      INTEGER i,j,iaux2
      REAL*8  workbot,workbet,worktot,diffbo(500),diffbe(500),diffto(500
     &     ),cost,radfact,totwork1,totwork2,lambo(500),lambe(500)
     &     ,lamto(500),lambda
      LOGICAL steeringf,steeringr 
*     ---- INITIALIZE AND SAVE BETWEEN CALLS ------------------------
      data diffbo, diffbe, diffto  /500*0.D0, 500*0.D0, 500*0.D0 /

C -- diffbo, etc are passed to steer_print and used as the value of last
C     call, so they must be saved when leaving and returning to this
C     routine:
      save diffbo, diffbe, diffto 


*==================== EXECUTABLE STATEMENTS ============================

      cost=1.d0 
      radfact=0.0174533d0 

      call steer_print(nbondsadded,force_cost,strvel
     &     ,strbonds,bondadded,eqdist,eqdist1,strtime0,strtime1
     &     ,rtime,workbo,lambo,diffbo,timesteer,cost,strtime0r
     &     ,strtime1r)

      call steer_print(nbendsadded,force_ang,ben_vel
     &     ,strbends,bendadded,eqang,eqang1,strtime0,strtime1,rtime
     &     ,workbe,lambe,diffbe,timesteer,radfact,strtime0r,strtime1r)
      
      call steer_print(nitorsadded,force_died,tor_vel
     &     ,strtors,torsadded,eqdied,eqdied1,strtime0,strtime1
     &     ,rtime,workto,lamto,diffto,timesteer,radfact,strtime0r
     &     ,strtime1r)

      IF(nbondsadded.gt.0) lambda = (lambo(1)-eqdist(1))/(eqdist1(1)
     &     -eqdist(1))
      IF(nbondsadded.eq.0.and.nbendsadded.gt.0) lambda = (lambe(1)
     &     /radfact-eqang(1))/(eqang1(1)-eqang(1))
      IF(nbondsadded.eq.0.and.nbendsadded.eq.0.and.nitorsadded.gt.0)
     &     lambda = (lamto(1)/radfact-eqdied(1))/(eqdied1(1)-eqdied(1))

c---  PRINTOUTS 
      steeringf= rtime.gt.strtime0.and.rtime.lt.strtime1
      steeringr= rtime.gt.strtime1r.and.rtime.lt.strtime0r

      iaux2=0
      IF((mod(ninn0,nplotstn0).EQ.0.AND.(steeringf.OR.steeringr)).
     &     or.
     &     (rtime.eq.strtime0.or.rtime.eq.strtime1)) THEN 
        workbot=0.d0
c---    do work on bonds
        do  i=1,nbondsadded
          if(strbonds(i).GT.0) THEN 
            iaux2=iaux2+1
            workbot=workbot+ workbo(i)
            if(iaux.lt.0) THEN 
              write(kplotst0,1015) i,rtime,lambo(i),bondadded(i)
     &             ,workbo(i)
            else
              write(kplotst0,1016) i,rtime,lambo(i),bondadded(i)
     &             ,workbo(i),iaux
            endif
1015        FORMAT(" bond n. ", i3,f13.1,4f15.5) 
1016        FORMAT(" bond n. ", i3,f13.1,3f15.5, "  Replica -> ",i4)   
          END IF
        end do

c---    do work on bends
        workbet=0.d0
        do i=1,nbendsadded
          if(strbends(i).GT.0) THEN 
            iaux2=iaux2+1
            workbet=workbet+ workbe(i)
            if(iaux.lt.0) THEN 
              write(kplotst0,1025)i,rtime,57.2958*lambe(i),57.2958
     &           *bendadded(i),workbe(i)
            else 
              write(kplotst0,1026)i,rtime,57.2958*lambe(i),57.2958
     &           *bendadded(i),workbe(i),iaux
            end if
1025        FORMAT(" bend n. ", i3,f13.1,3f15.5) 
1026        FORMAT(" bend n. ", i3,f13.1,3f15.5, "  Replica -> ",i4)   
          END IF
        end do

c---  do work on tors
        worktot=0.d0
        do i=1,nitorsadded
          if(strtors(i).GT.0) THEN 
            iaux2=iaux2+1
            worktot=worktot+ workto(i)
            if(iaux.lt.0) THEN 
              write(kplotst0,1035)i,rtime,57.2958*lamto(i),57.2958
     &             *torsadded(i),workto(i)
            else
              write(kplotst0,1036)i,rtime,57.2958*lamto(i),57.2958
     &             *torsadded(i),workto(i),iaux
            end if
1035        FORMAT(" tors n. ", i3,g13.2,3f15.5) 
1036        FORMAT(" tors n. ", i3,g13.2,3f15.5, "  Replica -> ",i4)   
          end if
        end do
        Totwork1 = workbot+workbet+worktot
        if(iaux2.gt.1) write(kplotst0,1045) rtime,lambda,workbot,workbet
     &       ,worktot,totwork1
      END IF
1045  FORMAT (" WORKfunction",f12.1,5f15.6)   
      
      RETURN
      END
