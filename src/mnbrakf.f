      SUBROUTINE mnbrakf(ax,bx,cx,fa,fb,fc,func,xp0,yp0,zp0,fpx,fpy,fpz)

!======================= DECLARATIONS ==================================

      IMPLICIT none

!----------------------- ARGUMENTS -------------------------------------

      REAL*8  xp0(*),yp0(*),zp0(*),fpx(*),fpy(*),fpz(*)

C     Routine taken from NUMERICAL RECIPIES

      EXTERNAL func
      REAL*8 func
      REAL*8 fa,fb,fc,ax,bx,cx,gold,glimit,tiny
      REAL*8 dum,fu,u,q,r,ulim
      parameter(gold=1.618034d0,glimit=100.d0,tiny=1.d-20)


      fa=func(ax,xp0,yp0,zp0,fpx,fpy,fpz)
      fb=func(bx,xp0,yp0,zp0,fpx,fpy,fpz)
      if(fb.gt.fa) then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
          fb=fa
          fa=dum
        endif
        cx=bx+gold*(bx-ax)
        fc=func(cx,xp0,yp0,zp0,fpx,fpy,fpz)
1       if(fb.ge.fc) then
          r=(bx-ax)*(fb-fc)
          q=(bx-cx)*(fb-fa)
          u=bx-((bx-cx)*q-(bx-ax)*r)/
     .         (2.d0*dsign(dmax1(dabs(q-r),tiny),q-r))
          ulim=bx+glimit*(cx-bx)
          if((bx-u)*(u-cx).gt.0.d0) then
            fu=func(u,xp0,yp0,zp0,fpx,fpy,fpz)
            if(fu.lt.fc) then
              ax=bx
              fa=fb
              bx=u
              fb=fu
              return
            else if(fu.gt.fb) then
              cx=u
              fc=fu
              return
            endif
            u=cx+gold*(cx-bx)
            fu=func(u,xp0,yp0,zp0,fpx,fpy,fpz)
          else if((cx-u)*(u-ulim).gt.0.d0) then
            fu=func(u,xp0,yp0,zp0,fpx,fpy,fpz)
            if(fu.lt.fc) then
              bx=cx
              cx=u
              u=cx+gold*(cx-bx)
              fb=fc
              fc=fu
              fu=func(u,xp0,yp0,zp0,fpx,fpy,fpz)
            endif
          else if((u-ulim)*(ulim-cx).ge.0.d0) then
            u=ulim
            fu=func(u,xp0,yp0,zp0,fpx,fpy,fpz)
          else
            u=cx+gold*(cx-bx)
            fu=func(u,xp0,yp0,zp0,fpx,fpy,fpz)
          endif
          ax=bx
          bx=cx
          cx=u
          fa=fb
          fb=fc
          fc=fu
          goto 1
        endif
        return
        end
