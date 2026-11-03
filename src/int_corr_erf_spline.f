c-------------------------------------------------------------------------
      subroutine int_corr_erf_spline(rmin,rmax,nbin,alpha,kcut,
     $       erf_arr_corr,tau)
c-------------------------------------------------------------------------
      implicit none
      integer nbin
      REAL*8 rmin,rmax,alpha,kcut,kcut_max,erf_arr_corr(4,*)
     &     ,tau(*),func,s,corr2 
      external func
      common /for_func/work,r

      REAL*8 del,r,pi,fac,erf,erft,work
      integer i,ibcbeg,ibcend
      work = alpha
      del = (rmax-rmin)/float(nbin) 
      pi = 3.14159265358979323846D0
      kcut_max=10.0D0*2.0D0*alpha
      write(6,2009) rmin,rmax,nbin
2009  FORMAT(6x,'doing erf table for intra term ( rmin =',
     &     f8.3, ' rmax=',f8.3,' nbin =',i6,')....' )
      do 100 i = 1,nbin+1
        r = del*(i-1) + rmin
        call QSIMP(FUNC,kcut,kcut_max,S) 
        tau(i) = r
        erf_arr_corr(1,i) = 2.0D0*s/pi
c$$$        write(6,110) r,erf_arr_corr(1,i)
c$$$110     format(2E15.5) 
100   continue
      ibcbeg = 0
      ibcend = 0
      write(6,*) '     erf table done ... '  
      call cubspl ( tau, erf_arr_corr,nbin+1, ibcbeg, ibcend )
      return
      end

      SUBROUTINE QSIMP(FUNC,A,B,S)                                              
      IMPLICIT NONE 
      INTEGER  j,JMAX
      REAL*8  FUNC,A,B,S,ST,OS,OST,EPS
      
      PARAMETER (EPS=1.D-6, JMAX=20)                                            
      OST=-1.D30                                                                
      OS= -1.D30                                                                
      DO 11 J=1,JMAX                                                            
        CALL TRAPZD(FUNC,A,B,ST,J)                                              
        S=(4.0D0*ST-OST)/3.0D0                                                        
        IF (ABS(S-OS).LT.EPS*ABS(OS)) RETURN                                    
        OS=S                                                                    
        OST=ST                                                                  
11    CONTINUE                                                                  
      RETURN 
      END                                                                       

      SUBROUTINE TRAPZD(FUNC,A,B,S,N)                                           
      IMPLICIT NONE 
      INTEGER  J,N,IT
      REAL*8  FUNC,A,B,S,TNM,DEL,X,sum
      SAVE it
      IF (N.EQ.1) THEN                                                          
        S=0.5*(B-A)*(FUNC(A)+FUNC(B))                                           
        IT=1                                                                    
      ELSE                                                                      
        TNM=IT                                                                  
        DEL=(B-A)/TNM                                                           
        X=A+0.5*DEL                                                             
        SUM=0.                                                                  
        DO 11 J=1,IT                                                            
          SUM=SUM+FUNC(X)                                                       
          X=X+DEL                                                               
11      CONTINUE                                                                
        S=0.5*(S+(B-A)*SUM/TNM)                                                 
        IT=2*IT                                                                 
      ENDIF                                                                     

      RETURN 
      END 

      FUNCTION FUNC(x) 
      implicit none 
      REAL*8  x,func,alpha,r
      common /for_func/alpha,r

      if(r.lt.1.d-8) THEN 
        func = dexp(-x**2.0D0/(4.0D0*alpha**2.0D0))
      ELSE
        func = dexp(-x**2.0D0/(4.D0*alpha**2.0D0)) * dsin(r*x)/(r*x) 
      END IF
      RETURN                                                                    
      END                                                                       
