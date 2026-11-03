      subroutine FNDFMT(icase,x,fmt)

**********************************************************************
*                                                                    *
*     Find FORTRAN-format of a given input number in character       *
*     string.                                                        *
*                                                                    *
*     G. Kneller                                                     *
*     IBM Corp., Data Systems Division, Dept. 48B                    *
*     Kingston, NY, USA                                              *
*                                                                    *
*---- Last update : 04/06/89 ----------------------------------------*
*                                                                      *
*     Written by Massimo Marchi IBM Corp., Kingston NY,  1989          *
*                                                                    *
*     Arguments :                                                    *
*                                                                    *
*     ICASE   : Flag describing variable type expected.     (INPUT)  *
*               >> integer ICASE <<                                  *
*               If ICASE = 1, if an integer number is expected,      *
*               If ICASE = 2, if a real number is expected.          *
*               In all other cases a string variable is assumed.     *
*     X       : Input number in character string.           (INPUT)  *
*               >> character*80 X <<                                 *
*     FMT     : FORTRAN-format of X.                        (OUTPUT) *
*               >> character*8 FMT <<                                *
*                                                                    *
*                                                                    *
*     Written by Gerald Kneller Dept 48B, IBM Kingston 1989          *
*                                                                    *
**********************************************************************

*====================================================================*
*     DECLARATIONS :                                                 *
*====================================================================*

      IMPLICIT CHARACTER*80(a-z)

*     ARGUMENTS :

      integer      icase,indexd
      character*80 x
      character*8  fmt

*     CONSTANTS:
      integer npatti,npattr
      parameter(npatti=2,npattr=48)

*     LOCAL VARIABLES :

      integer      width,i,j,nt,digaft,
     .             multi(80)
      logical      known
      character*1  blank,
     .             xarray(80),digit(10),ctypes(0:80)
      character*80 match,
     .             patti(npatti),pattr(npattr)

*     DATA STATEMENTS :

      data digit /'1','2','3','4','5','6','7','8','9','0'/
      data blank /' '/
      data patti /'D','sD'/
      data pattr / 'D'      ,'pD'      ,'Dp'      ,'DpD',
     .             'DsD'    ,'pDsD'    ,'DpsD'    ,'DpDsD',
     .             'DeD'    ,'pDeD'    ,'DpeD'    ,'DpDeD',
     .             'DesD'   ,'pDesD'   ,'DpesD'   ,'DpDesD',
     .             'DeBD'   ,'pDeBD'   ,'DpeBD'   ,'DpDeBD',
     .             'DeBsD'  ,'pDeBsD'  ,'DpeBsD'  ,'DpDeBsD',
     .            'sD'     ,'spD'     ,'sDp'     ,'sDpD',
     .            'sDsD'   ,'spDsD'   ,'sDpsD'   ,'sDpDsD',
     .            'sDeD'   ,'spDeD'   ,'sDpeD'   ,'sDpDeD',
     .            'sDesD'  ,'spDesD'  ,'sDpesD'  ,'sDpDesD',
     .            'sDeBD'  ,'spDeBD'  ,'sDpeBD'  ,'sDpDeBD',
     .            'sDeBsD' ,'spDeBsD' ,'sDpeBsD' ,'sDpDeBsD'/

*====================================================================*
*     Executable statements :                                        *
*====================================================================*

*---- COPY INPUT STRING ON ARRAY: -----------------------------------*
      read(x,'(80a1)')(xarray(i),i=1,80)

*---- FIND FIELD WIDTH w: -------------------------------------------*
      width = 0
      do 100 i = 1,80
 100     if (xarray(i) .ne. blank) width = i
      if (width .eq. 0) then
	 IF(icase.EQ.1) THEN
            fmt = '(i80)'
	 ELSE IF(icase.EQ.2) THEN
            fmt = '(f80.10)'
         ELSE
	    fmt='(a80)'
         END IF
         return
      endif

*---- IF A CHARACTER VARIABLE IS EXPECTED: --------------------------*
      if (icase .ne. 1 .and. icase .ne. 2) then
         if (width .lt. 10) then
            write(fmt,'(a2,i1,a1)')'(A',width,')'
            return
         else
            write(fmt,'(a2,i2,a1)')'(A',width,')'
            return
         endif
      endif

*---- IF AN INTEGER OR A REAL NUMBER IS EXPECTED: -------------------*

*.....FIND PATTERN OF INPUT STRING:
      nt = 0
      ctypes(0) = blank
      do 200 i = 1,width
         known = .false.
         if (xarray(i) .eq. blank) then
            known = .true.
            nt = nt + 1
            if (nt .gt. 1 .and. ctypes(nt-1) .eq. 'B') then
               nt = nt - 1
               multi(nt) = multi(nt) + 1
            else
               multi(nt) = 1
               ctypes(nt) = 'B'
            GOTO 200
            endif
         endif
         do 210 j = 1,10
            if (xarray(i) .eq. digit(j)) then
               known = .true.
               nt = nt + 1
               if (nt .gt. 1 .and. ctypes(nt-1) .eq. 'D') then
                  nt = nt - 1
                  multi(nt) = multi(nt) + 1
               else
                  multi(nt) = 1
                  ctypes(nt) = 'D'
               GOTO 200
               endif
            endif
 210     continue
         if (xarray(i) .eq. '+' .or. xarray(i) .eq. '-') then
            known = .true.
            nt = nt + 1
            ctypes(nt) = 's'
            multi(nt) = 1
            GOTO 200
         endif
         if (xarray(i) .eq. '.') then
            known = .true.
            nt = nt + 1
            ctypes(nt) = 'p'
            multi(nt) = 1
            GOTO 200
         endif
         if (xarray(i) .eq. 'E' .or. xarray(i) .eq. 'D' .or.
     .       xarray(i) .eq. 'e' .or. xarray(i) .eq. 'd'     ) then
            known = .true.
            nt = nt + 1
            ctypes(nt) = 'e'
            multi(nt) = 1
            GOTO 200
         endif
         if (.not. known) then
            nt = nt + 1
            if (nt .gt. 1 .and. ctypes(nt-1) .eq. '?') then
               nt = nt - 1
               multi(nt) = multi(nt) + 1
            else
               multi(nt) = 1
               ctypes(nt) = '?'
            endif
         endif
 200  continue

* for testing purp.:
*     do 111 i = 1,nt
*        print *,ctypes(i),multi(i)
*111  continue
*

      if (ctypes(1) .ne. 'B') then
         write(match,'(80a1)')(ctypes(i),i=1,nt)
      else
         write(match,'(80a1)')(ctypes(i),i=2,nt)
      endif

*.....IN CASE THAT AN INTEGER NUMBER IS EXPECTED:
      if (icase .eq. 1) then
         do 300 i = 1,npatti
            if (match .eq. patti(i)) GOTO 1000
 300     continue
         fmt = '?'
         return
1000     if (width .lt. 10) then
            write(fmt,'(a2,i1,a1)')'(I',width,')'
            return
         else
            write(fmt,'(a2,i2,a1)')'(I',width,')'
            return
         endif
      endif

*.....IN CASE THAT A REAL NUMBER IS EXPECTED:
      if (icase .eq. 2) then
         do 310 i = 1,npattr
            if (match .eq. pattr(i)) GOTO 1100
 310     continue
         fmt = '?'
         return
1100     if (INDEXD(match,'pD') .ne. 0) then
            do 320 i = 1,nt
               if (ctypes(i) .eq. 'p') then
                  digaft = multi(i+1)
                  GOTO 1200
               endif
 320        continue
         else
            digaft = 0
            GOTO 1200
         endif
1200     if (width .lt. 10) then
            write(fmt,'(a2,i1,a1,i1,a1)')'(F',width,'.',digaft,')'
            return
         endif
         if (width .ge. 10 .and. digaft .lt. 10) then
            write(fmt,'(a2,i2,a1,i1,a1)')'(F',width,'.',digaft,')'
            return
         endif
         if (width .ge. 10 .and. digaft .ge. 10) then
            write(fmt,'(a2,i2,a1,i2,a1)')'(F',width,'.',digaft,')'
            return
         endif
      endif

*---- END OF SR -----------------------------------------------------*

      END

      INTEGER FUNCTION INDEXD(string,subs)
      character*80 string
      character*2  subs
      indexd=0
      do i=1,79
        if(subs.eq.string(I:I+1)) indexd=1
      end do
      return
      end 


      
