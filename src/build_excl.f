      SUBROUTINE build_excl(nat, xp0, yp0, zp0, jret)
*=======================================================================
*     Build the flat exclusion list (1-2, 1-3, 1-4) for the pair-wise
*     nonbonded force, and check that every excluded pair is within the
*     M-shell cutoff.  Called once from mtsmd, at the top, after the
*     topology is assembled and the coordinates are loaded.
*
*     The lists lbnd, int13, int14 are in the parst module and are the
*     full-system, correctly-offset lists (solute + solvent), built by
*     join/add_solvent_tpg.  No de-duplication is needed: clt14
*     guarantees int14 is disjoint from int13 and unique, and 1-2 / 1-3
*     are disjoint by construction.  The one residual case (a 1-2 pair
*     that is also a 1-4, possible only in a 3-membered ring) is
*     filtered below.
*
*     On output, nexcl, excl_i1, excl_i2 are in the cpropar module.
*     iret = 1 if a dimension or distance check fails.
*=======================================================================
      use parst
      use cpropar
      use unit, only: kprint
      IMPLICIT none
      INTEGER nat,jret
      REAL*8 xp0(nat),yp0(nat),zp0(nat) 
      INTEGER i1, i2, k, m, lj, li, lij
      INTEGER ig, ia, ib
      LOGICAL ldum
      REAL*8  xa, ya, za, rsq, r0_excl
*=======================================================================

      jret = 0

*----- allocate
      IF (ALLOCATED(excl_i1)) DEALLOCATE(excl_i1)
      IF (ALLOCATED(excl_i2)) DEALLOCATE(excl_i2)
      nexcl = lbond + int13p + int14p
      ALLOCATE(excl_i1(nexcl), excl_i2(nexcl), stat=jret)
      IF (jret .NE. 0) THEN
         CALL xerror('build_excl: allocation failed',80,1,2)
         RETURN
      END IF
      nexcl = 0

*----- 1-2 (bonds): full strength
      DO m = 1, lbond
         nexcl = nexcl + 1
         excl_i1(nexcl) = lbnd(1,m)
         excl_i2(nexcl) = lbnd(2,m)
      END DO

*----- 1-3 (bends): full strength
      DO m = 1, int13p
         nexcl = nexcl + 1
         excl_i1(nexcl) = int13(1,m)
         excl_i2(nexcl) = int13(2,m)
      END DO

*----- 1-4 (torsion endpoints): full strength (fnb14 supplies the fudge)
*     skip a 1-4 pair that is also a 1-2 (3-membered ring pathological case)
      DO m = 1, int14p
         i1 = int14(1,m)
         i2 = int14(2,m)
         ldum = .FALSE.
         DO k = 1, lbond
            IF ((lbnd(1,k).EQ.i1.AND.lbnd(2,k).EQ.i2).OR.
     &          (lbnd(1,k).EQ.i2.AND.lbnd(2,k).EQ.i1)) THEN
               ldum = .TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT.ldum) THEN
            nexcl = nexcl + 1
            excl_i1(nexcl) = i1
            excl_i2(nexcl) = i2
         END IF
      END DO

      WRITE(kprint,'(5x,a,i8,a,i8,a)')
     &   'Exclusion list built: nexcl = ',nexcl,
     &   '  (nominal ',lbond+int13p+int14p,')'

*----- optional: assert every intragroup pair is excluded (so ingrpp=0
*     is exact on the GPU).  For TIP3P water and AMBER proteins this
*     should not fire.
      DO ig = 1, ngrp
         DO ia = grppt(1,ig), grppt(2,ig)
            DO ib = ia+1, grppt(2,ig)
               ldum = .FALSE.
               DO k = 1, nexcl
                  IF ((excl_i1(k).EQ.ia.AND.excl_i2(k).EQ.ib).OR.
     &                (excl_i1(k).EQ.ib.AND.excl_i2(k).EQ.ia)) THEN
                     ldum = .TRUE.
                     EXIT
                  END IF
               END DO
               IF (.NOT.ldum) THEN
                 WRITE(kprint,*) 'WARNING intragroup pair not excluded:'
     &               ,ia,ib,' group ',ig
               END IF
            END DO
         END DO
      END DO

*----- distance check: every excluded pair must be within r0_excl.
*     Excluded pairs are intramolecular and never split, so no PBC.
*     xp0 is in the ORTHOGONAL frame (co-transformed), so the raw
*     difference is the physical displacement.
      r0_excl = rcutm - rtolm
      IF (r0_excl .LE. 0.d0) r0_excl = rcutm
      r0_excl = MIN(r0_excl, 4.5d0)
      DO k = 1, nexcl
         i1 = excl_i1(k)
         i2 = excl_i2(k)
         xa = xp0(i1) - xp0(i2)
         ya = yp0(i1) - yp0(i2)
         za = zp0(i1) - zp0(i2)
         rsq = xa*xa + ya*ya + za*za
         IF (rsq .GT. r0_excl**2) THEN
            WRITE(kprint,'(a,i7,a,i7,a,f10.4,a,f10.4)')
     &          ' ERROR: excluded pair (',i1,',',i2,') at r = ',
     &          SQRT(rsq),' exceeds r0 = ',r0_excl
            CALL xerror('Excluded contact beyond M-shell cutoff.'//
     &           ' Check topology or increase r_M.',80,1,2)
            jret = 1
            RETURN
          END IF
      END DO
      WRITE(kprint,'(5x,a,i8,a,f8.3)')
     &   'Excluded-contact distance check passed for ',nexcl,
     &   ' pairs, r0 = ',r0_excl

      RETURN
      END
      
