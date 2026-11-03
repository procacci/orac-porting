subroutine cache_line(n,nthr,lcl,cl_size)
  implicit NONE
  integer n,nthr,cl_size,nchunks,i,j,k,l,nbalance(nthr),m
  logical*1 lcl(n,nthr)
  write(206,*) n,nthr
  nchunks=1+n/cl_size
  write(206,*) "lc",nchunks,n,nthr
  lcl=.false.
  do j=1,nthr
     nbalance(j)=0
     i=0
     do k=j,nchunks,nthr
        if(i.ge.n) exit
        do l=1,cl_size
           i=l+(k-1)*cl_size
           if(i.le.n) THEN 
              lcl(i,j)=.true.
           END IF
           write(206,*) i,l,k,j,(lcl(i,m),m=1,nthr)
           nbalance(j)=nbalance(j)+1
        end do
     end do
  end do
  do i=1,nthr
     write(206,*) nbalance(i)
  end do
  do i=1,n
     write(206,10) (lcl(i,j),j=1,nthr)
10   format(16(1x,L1,1x))
  end do
  return
end subroutine cache_line
