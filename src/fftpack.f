** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                  version 3  june 1979
*  a package of fortran subprograms for the fast fourier
*      transform of periodic and other symmetric sequences
*  paul n swarztrauber
*  national center for atmospheric research  boulder,colorado 80307
*   which is sponsored by the national science foundation
*  modified by P. Bjorstad
*
** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*1.   drffti    initialize  drfftf and drfftb
*2.   drfftf    forward transform of a real periodic sequence
*3.   drfftb    backward transform of a real coefficient array
*
*4.   defftf    a simplified real periodic forward transform
*5.   defftb    a simplified real periodic backward transform
*
*6.   dsinti    initialize dsint
*7.   dsint     sine transform of a real odd sequence
*
*8.   dcosti    initialize dcost
*9.   dcost     cosine transform of a real even sequence
*
*10.  dsinqi    initialize dsinqf and dsinqb
*11.  dsinqf    forward sine transform with odd wave numbers
*12.  dsinqb    unnormalized inverse of dsinqf
*
*13.  dcosqi    initialize dcosqf and dcosqb
*14.  dcosqf    forward cosine transform with odd wave numbers
*15.  dcosqb    unnormalized inverse of dcosqf
*
*16.  dcffti     initialize dcfftf and dcfftb
*17.  dcfftf     forward transform of a complex periodic sequence
*18.  dcfftb     unnormalized inverse of dcfftf
*
*For single precision, change initial d to s.
*Each subroutine is described below. the names used refer to
*the REAL*8 version, but the same description
*applies for the single precision version.
*
*****************************************************************
*
*subroutine drffti(n,wsave)
*
*****************************************************************
*
*subroutine drffti initializes the array wsave which is used in
*both drfftf and drfftb. the prime factorization of n together with
*a tabulation of the trigonometric functions are computed and
*stored in wsave.
*
*input parameter
*
*n       the length of the sequence to be transformed.
*
*output parameter
*
*wsave   a work array which must be dimensioned at least 2*n+15.
*        the same work array can be used for both drfftf and drfftb
*        as long as n remains unchanged. different wsave arrays
*        are required for different values of n. the contents of
*        wsave must not be changed between calls of drfftf or drfftb.
*
*******************************************************************
*
*subroutine drfftf(n,r,wsave)
*
*******************************************************************
*
*subroutine drfftf computes the fourier coefficients of a real
*perodic sequence (fourier analysis). the transform is defined
*below at output parameter r.
*
*input parameters
*
*n       the length of the array r to be transformed.  the method
*        is most efficient when n is a product of small primes.
*        n may change so long as different work arrays are provided
*
*r       a real array of length n which contains the sequence
*        to be transformed
*
*wsave   a work array which must be dimensioned at least 2*n+15.
*        in the program that calls drfftf. the wsave array must be
*        initialized by calling subroutine drffti(n,wsave) and a
*        different wsave array must be used for each different
*        value of n. this initialization does not have to be
*        repeated so long as n remains unchanged thus subsequent
*        transforms can be obtained faster than the first.
*        the same wsave array can be used by drfftf and drfftb.
*
*
*output parameters
*
*r       r(1) = the sum from i=1 to i=n of r(i)
*
*        if n is even set l =n/2   , if n is odd set l = (n+1)/2
*
*          then for k = 2,...,l
*
*             r(2*k-2) = the sum from i = 1 to i = n of
*
*                  r(i)*cos((k-1)*(i-1)*2*pi/n)
*
*             r(2*k-1) = the sum from i = 1 to i = n of
*
*                 -r(i)*sin((k-1)*(i-1)*2*pi/n)
*
*        if n is even
*
*             r(n) = the sum from i = 1 to i = n of
*
*                  (-1)**(i-1)*r(i)
*
* *****  note
*             this transform is unnormalized since a call of drfftf
*             followed by a call of drfftb will multiply the input
*             sequence by n.
*
*wsave   contains results which must not be destroyed between
*        calls of drfftf or drfftb.
*
*
*******************************************************************
*
*subroutine drfftb(n,r,wsave)
*
*******************************************************************
*
*subroutine drfftb computes the real perodic sequence from its
*fourier coefficients (fourier synthesis). the transform is defined
*below at output parameter r.
*
*input parameters
*
*n       the length of the array r to be transformed.  the method
*        is most efficient when n is a product of small primes.
*        n may change so long as different work arrays are provided
*
*r       a real array of length n which contains the sequence
*        to be transformed
*
*wsave   a work array which must be dimensioned at least 2*n+15.
*        in the program that calls drfftb. the wsave array must be
*        initialized by calling subroutine drffti(n,wsave) and a
*        different wsave array must be used for each different
*        value of n. this initialization does not have to be
*        repeated so long as n remains unchanged thus subsequent
*        transforms can be obtained faster than the first.
*        the same wsave array can be used by drfftf and drfftb.
*
*
*output parameters
*
*r       for n even and for i = 1,...,n
*
*             r(i) = r(1)+(-1)**(i-1)*r(n)
*
*                  plus the sum from k=2 to k=n/2 of
*
*                   2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)
*
*                  -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)
*
*        for n odd and for i = 1,...,n
*
*             r(i) = r(1) plus the sum from k=2 to k=(n+1)/2 of
*
*                  2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)
*
*                 -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)
*
* *****  note
*             this transform is unnormalized since a call of drfftf
*             followed by a call of drfftb will multiply the input
*             sequence by n.
*
*wsave   contains results which must not be destroyed between
*        calls of drfftb or drfftf.
*
*
*******************************************************************
*
*subroutine defftf(n,r,azero,a,b,wsave)
*
*******************************************************************
*
*subroutine defftf computes the fourier coefficients of a real
*perodic sequence (fourier analysis). the transform is defined
*below at output parameters azero,a and b. defftf is a simplified
*version of drfftf. it is not as fast as drfftf since scaling
*and initialization are computed for each transform. the repeated
*initialization can be suppressed by removeing the statment
*( call deffti(n,wsave) ) from both defftf and defftb and inserting
*it at the appropriate place in your program.
*
*input parameters
*
*n       the length of the array r to be transformed.  the method
*        is must efficient when n is the product of small primes.
*
*r       a real array of length n which contains the sequence
*        to be transformed. r is not destroyed.
*
*wsave   a work array with at least 3*n+15 locations.
*
*output parameters
*
*azero   the sum from i=1 to i=n of r(i)/n
*
*a,b     for n even b(n/2)=0. and a(n/2) is the sum from i=1 to
*        i=n of (-1)**(i-1)*r(i)/n
*
*        for n even define kmax=n/2-1
*        for n odd  define kmax=(n-1)/2
*
*        then for  k=1,...,kmax
*
*             a(k) equals the sum from i=1 to i=n of
*
*                  2./n*r(i)*cos(k*(i-1)*2*pi/n)
*
*             b(k) equals the sum from i=1 to i=n of
*
*                  2./n*r(i)*sin(k*(i-1)*2*pi/n)
*
*
*******************************************************************
*
*subroutine defftb(n,r,azero,a,b,wsave)
*
*******************************************************************
*
*subroutine defftb computes a real perodic sequence from its
*fourier coefficients (fourier synthesis). the transform is
*defined below at output parameter r. defftb is a simplified
*version of drfftb. it is not as fast as drfftb since scaling and
*initialization are computed for each transform. the repeated
*initialization can be suppressed by removeing the statment
*( call deffti(n,wsave) ) from both defftf and defftb and inserting
*( call deffti(n,wsave) ) from both defftf and defftb and inserting
*it at the appropriate place in your program.
*
*input parameters
*
*n       the length of the output array r.  the method is most
*        efficient when n is the product of small primes.
*
*azero   the constant fourier coefficient
*
*a,b     arrays which contain the remaining fourier coefficients
*        these arrays are not destroyed.
*
*        the length of these arrays depends on whether n is even or
*        odd.
*
*        if n is even n/2    locations are required
*        if n is odd (n-1)/2 locations are required
*
*wsave   a work array with at least 3*n+15 locations.
*
*
*output parameters
*
*r       if n is even define kmax=n/2
*        if n is odd  define kmax=(n-1)/2
*
*        then for i=1,...,n
*
*             r(i)=azero plus the sum from k=1 to k=kmax of
*
*             a(k)*cos(k*(i-1)*2*pi/n)+b(k)*sin(k*(i-1)*2*pi/n)
*
********************** complex notation **************************
*
*        for j=1,...,n
*
*        r(j) equals the sum from k=-kmax to k=kmax of
*
*             c(k)*exp(i*k*(j-1)*2*pi/n)
*
*        where
*
*             c(k) = .5*cmplx(a(k),-b(k))   for k=1,...,kmax
*
*             c(-k) = conjg(c(k))
*
*             c(0) = azero
*
*                  and i=sqrt(-1)
*
**************** amplitude - phase notation ***********************
*
*        for i=1,...,n
*
*        r(i) equals azero plus the sum from k=1 to k=kmax of
*
*             alpha(k)*cos(k*(i-1)*2*pi/n+beta(k))
*
*        where
*
*             alpha(k) = sqrt(a(k)*a(k)+b(k)*b(k))
*
*             cos(beta(k))=a(k)/alpha(k)
*
*             sin(beta(k))=-b(k)/alpha(k)
*
*******************************************************************
*
*subroutine dsinti(n,wsave)
*
*******************************************************************
*
*subroutine dsinti initializes the array wsave which is used in
*subroutine dsint. the prime factorization of n together with
*a tabulation of the trigonometric functions are computed and
*stored in wsave.
*
*input parameter
*
*n       the length of the sequence to be transformed.  the method
*        is most efficient when n+1 is a product of small primes.
*
*output parameter
*
*wsave   a work array with at least int(2.5*n+15) locations.
*        different wsave arrays are required for different values
*        of n. the contents of wsave must not be changed between
*        calls of dsint.
*
*******************************************************************
*
*subroutine dsint(n,x,wsave)
*
*******************************************************************
*
*subroutine dsint computes the discrete fourier sine transform
*of an odd sequence x(i). the transform is defined below at
*output parameter x.
*
*dsint is the unnormalized inverse of itself since a call of dsint
*followed by another call of dsint will multiply the input sequence
*x by 2*(n+1).
*
*the array wsave which is used by subroutine dsint must be
*initialized by calling subroutine dsinti(n,wsave).
*
*input parameters
*
*n       the length of the sequence to be transformed.  the method
*        is most efficient when n+1 is the product of small primes.
*
*x       an array which contains the sequence to be transformed
*
*              ************important*************
*
*              x must be dimensioned at least n+1
*
*wsave   a work array with dimension at least int(2.5*n+15)
*        in the program that calls dsint. the wsave array must be
*        initialized by calling subroutine dsinti(n,wsave) and a
*        different wsave array must be used for each different
*        value of n. this initialization does not have to be
*        repeated so long as n remains unchanged thus subsequent
*        transforms can be obtained faster than the first.
*
*output parameters
*
*x       for i=1,...,n
*
*             x(i)= the sum from k=1 to k=n
*
*                  2*x(k)*sin(k*i*pi/(n+1))
*
*             a call of dsint followed by another call of
*             dsint will multiply the sequence x by 2*(n+1).
*             hence dsint is the unnormalized inverse
*             of itself.
*
*wsave   contains initialization calculations which must not be
*        destroyed between calls of dsint.
*
*******************************************************************
*
*subroutine dcosti(n,wsave)
*
*******************************************************************
*
*subroutine dcosti initializes the array wsave which is used in
*subroutine dcost. the prime factorization of n together with
*a tabulation of the trigonometric functions are computed and
*stored in wsave.
*
*input parameter
*
*n       the length of the sequence to be transformed.  the method
*        is most efficient when n-1 is a product of small primes.
*
*output parameter
*
*wsave   a work array which must be dimensioned at least 3*n+15.
*        different wsave arrays are required for different values
*        of n. the contents of wsave must not be changed between
*        calls of dcost.
*
*******************************************************************
*
*subroutine dcost(n,x,wsave)
*
*******************************************************************
*
*subroutine dcost computes the discrete fourier cosine transform
*of an even sequence x(i). the transform is defined below at output
*parameter x.
*
*dcost is the unnormalized inverse of itself since a call of dcost
*followed by another call of dcost will multiply the input sequence
*x by 2*(n-1). the transform is defined below at output parameter x
*
*the array wsave which is used by subroutine dcost must be
*initialized by calling subroutine dcosti(n,wsave).
*
*input parameters
*
*n       the length of the sequence x. n must be greater than 1.
*        the method is most efficient when n-1 is a product of
*        small primes.
*
*x       an array which contains the sequence to be transformed
*
*wsave   a work array which must be dimensioned at least 3*n+15
*        in the program that calls dcost. the wsave array must be
*        initialized by calling subroutine dcosti(n,wsave) and a
*        different wsave array must be used for each different
*        value of n. this initialization does not have to be
*        repeated so long as n remains unchanged thus subsequent
*        transforms can be obtained faster than the first.
*
*output parameters
*
*x       for i=1,...,n
*
*           x(i) = x(1)+(-1)**(i-1)*x(n)
*
*             + the sum from k=2 to k=n-1
*
*               2*x(k)*cos((k-1)*(i-1)*pi/(n-1))
*
*             a call of dcost followed by another call of
*             dcost will multiply the sequence x by 2*(n-1)
*             hence dcost is the unnormalized inverse
*             of itself.
*
*wsave   contains initialization calculations which must not be
*        destroyed between calls of dcost.
*
*******************************************************************
*
*subroutine dsinqi(n,wsave)
*
*******************************************************************
*
*subroutine dsinqi initializes the array wsave which is used in
*both dsinqf and dsinqb. the prime factorization of n together with
*a tabulation of the trigonometric functions are computed and
*stored in wsave.
*
*input parameter
*
*n       the length of the sequence to be transformed. the method
*        is most efficient when n is a product of small primes.
*
*output parameter
*
*wsave   a work array which must be dimensioned at least 3*n+15.
*        the same work array can be used for both dsinqf and dsinqb
*        as long as n remains unchanged. different wsave arrays
*        are required for different values of n. the contents of
*        wsave must not be changed between calls of dsinqf or
*        dsinqb.
*
*******************************************************************
*
*subroutine dsinqf(n,x,wsave)
*
*******************************************************************
*
*subroutine dsinqf computes the fast fourier transform of quarter
*wave data. that is , dsinqf computes the coefficients in a sine
*series representation with only odd wave numbers. the transform
*is defined below at output parameter x.
*
*dsinqb is the unnormalized inverse of dsinqf -- a call of dsinqf
*followed by a call of dsinqb will multiply the input sequence x
*by 4*n.
*
*the array wsave which is used by subroutine dsinqf must be
*initialized by calling subroutine dsinqi(n,wsave).
*
*
*input parameters
*
*n       the length of the array x to be transformed.  the method
*        is most efficient when n is a product of small primes.
*
*x       an array which contains the sequence to be transformed
*
*wsave   a work array which must be dimensioned at least 3*n+15.
*        in the program that calls dsinqf. the wsave array must be
*        initialized by calling subroutine dsinqi(n,wsave) and a
*        different wsave array must be used for each different
*        value of n. this initialization does not have to be
*        repeated so long as n remains unchanged thus subsequent
*        transforms can be obtained faster than the first.
*
*output parameters
*
*x       for i=1,...,n
*
*             x(i) = (-1)**(i-1)*x(n)
*
*                + the sum from k=1 to k=n-1 of
*
*                2*x(k)*sin((2*i-1)*k*pi/(2*n))
*
*             a call of dsinqf followed by a call of
*             dsinqb will multiply the sequence x by 4*n.
*             therefore dsinqb is the unnormalized inverse
*             of dsinqf.
*
*wsave   contains initialization calculations which must not
*        be destroyed between calls of dsinqf or dsinqb.
*
*******************************************************************
*
*subroutine dsinqb(n,x,wsave)
*
*******************************************************************
*
*subroutine dsinqb computes the fast fourier transform of quarter
*wave data. that is , dsinqb computes a sequence from its
*representation in terms of a sine series with odd wave numbers.
*the transform is defined below at output parameter x.
*
*dsinqf is the unnormalized inverse of dsinqb -- a call of dsinqb
*followed by a call of dsinqf will multiply the input sequence x
*by 4*n.
*
*the array wsave which is used by subroutine dsinqb must be
*initialized by calling subroutine dsinqi(n,wsave).
*
*
*input parameters
*
*n       the length of the array x to be transformed.  the method
*        is most efficient when n is a product of small primes.
*
*x       an array which contains the sequence to be transformed
*
*wsave   a work array which must be dimensioned at least 3*n+15.
*        in the program that calls dsinqb. the wsave array must be
*        initialized by calling subroutine dsinqi(n,wsave) and a
*        different wsave array must be used for each different
*        value of n. this initialization does not have to be
*        repeated so long as n remains unchanged thus subsequent
*        transforms can be obtained faster than the first.
*
*output parameters
*
*x       for i=1,...,n
*
*             x(i)= the sum from k=1 to k=n of
*
*               4*x(k)*sin((2k-1)*i*pi/(2*n))
*
*             a call of dsinqb followed by a call of
*             dsinqf will multiply the sequence x by 4*n.
*             therefore dsinqf is the unnormalized inverse
*             of dsinqb.
*
*wsave   contains initialization calculations which must not
*        be destroyed between calls of dsinqb or dsinqf.
*
*******************************************************************
*
*subroutine dcosqi(n,wsave)
*
*******************************************************************
*
*subroutine dcosqi initializes the array wsave which is used in
*both dcosqf and dcosqb. the prime factorization of n together with
*a tabulation of the trigonometric functions are computed and
*stored in wsave.
*
*input parameter
*
*n       the length of the array to be transformed.  the method
*        is most efficient when n is a product of small primes.
*
*output parameter
*
*wsave   a work array which must be dimensioned at least 3*n+15.
*        the same work array can be used for both dcosqf and dcosqb
*        as long as n remains unchanged. different wsave arrays
*        are required for different values of n. the contents of
*        wsave must not be changed between calls of dcosqf or
*        dcosqb.
*
*******************************************************************
*
*subroutine dcosqf(n,x,wsave)
*
*******************************************************************
*
*subroutine dcosqf computes the fast fourier transform of quarter
*wave data. that is , dcosqf computes the coefficients in a cosine
*series representation with only odd wave numbers. the transform
*is defined below at output parameter x
*
*dcosqf is the unnormalized inverse of dcosqb -- a call of dcosqf
*followed by a call of dcosqb will multiply the input sequence x
*by 4*n.
*
*the array wsave which is used by subroutine dcosqf must be
*initialized by calling subroutine dcosqi(n,wsave).
*
*
*input parameters
*
*n       the length of the array x to be transformed.  the method
*        is most efficient when n is a product of small primes.
*
*x       an array which contains the sequence to be transformed
*
*wsave   a work array which must be dimensioned at least 3*n+15
*        in the program that calls dcosqf. the wsave array must be
*        initialized by calling subroutine dcosqi(n,wsave) and a
*        different wsave array must be used for each different
*        value of n. this initialization does not have to be
*        repeated so long as n remains unchanged thus subsequent
*        transforms can be obtained faster than the first.
*
*output parameters
*
*x       for i=1,...,n
*
*             x(i) = x(1) plus the sum from k=2 to k=n of
*
*                2*x(k)*cos((2*i-1)*(k-1)*pi/(2*n))
*
*             a call of dcosqf followed by a call of
*             dcosqb will multiply the sequence x by 4*n.
*             therefore dcosqb is the unnormalized inverse
*             of dcosqf.
*
*wsave   contains initialization calculations which must not
*        be destroyed between calls of dcosqf or dcosqb.
*
*******************************************************************
*
*subroutine dcosqb(n,x,wsave)
*
*******************************************************************
*
*subroutine dcosqb computes the fast fourier transform of quarter
*wave data. that is , dcosqb computes a sequence from its
*representation in terms of a cosine series with odd wave numbers.
*the transform is defined below at output parameter x.
*
*dcosqb is the unnormalized inverse of dcosqf -- a call of dcosqb
*followed by a call of dcosqf will multiply the input sequence x
*by 4*n.
*
*the array wsave which is used by subroutine dcosqb must be
*initialized by calling subroutine dcosqi(n,wsave).
*
*
*input parameters
*
*n       the length of the array x to be transformed.  the method
*        is most efficient when n is a product of small primes.
*
*x       an array which contains the sequence to be transformed
*
*wsave   a work array that must be dimensioned at least 3*n+15
*        in the program that calls dcosqb. the wsave array must be
*        initialized by calling subroutine dcosqi(n,wsave) and a
*        different wsave array must be used for each different
*        value of n. this initialization does not have to be
*        repeated so long as n remains unchanged thus subsequent
*        transforms can be obtained faster than the first.
*
*output parameters
*
*x       for i=1,...,n
*
*             x(i)= the sum from k=1 to k=n of
*
*               4*x(k)*cos((2*k-1)*(i-1)*pi/(2*n))
*
*             a call of dcosqb followed by a call of
*             dcosqf will multiply the sequence x by 4*n.
*             therefore dcosqf is the unnormalized inverse
*             of dcosqb.
*
*wsave   contains initialization calculations which must not
*        be destroyed between calls of dcosqb or dcosqf.
*
*******************************************************************
*
*subroutine dcffti(n,wsave)
*
*******************************************************************
*
*subroutine dcffti initializes the array wsave which is used in
*both dcfftf and dcfftb. the prime factorization of n together with
*a tabulation of the trigonometric functions are computed and
*stored in wsave.
*
*input parameter
*
*n       the length of the sequence to be transformed
*
*output parameter
*
*wsave   a work array which must be dimensioned at least 4*n+15
*        the same work array can be used for both dcfftf and dcfftb
*        as long as n remains unchanged. different wsave arrays
*        are required for different values of n. the contents of
*        wsave must not be changed between calls of dcfftf or dcfftb.
*
*******************************************************************
*
*subroutine dcfftf(n,c,wsave)
*
*******************************************************************
*
*subroutine dcfftf computes the forward complex discrete fourier
*transform (the fourier analysis). equivalently , dcfftf computes
*the fourier coefficients of a complex periodic sequence.
*the transform is defined below at output parameter c.
*
*the transform is not normalized. to obtain a normalized transform
*the output must be divided by n. otherwise a call of dcfftf
*followed by a call of dcfftb will multiply the sequence by n.
*
*the array wsave which is used by subroutine dcfftf must be
*initialized by calling subroutine dcffti(n,wsave).
*
*input parameters
*
*
*n      the length of the complex sequence c. the method is
*       more efficient when n is the product of small primes. n
*
*c      a complex array of length n which contains the sequence
*
*wsave   a real work array which must be dimensioned at least 4n+15
*        in the program that calls dcfftf. the wsave array must be
*        initialized by calling subroutine dcffti(n,wsave) and a
*        different wsave array must be used for each different
*        value of n. this initialization does not have to be
*        repeated so long as n remains unchanged thus subsequent
*        transforms can be obtained faster than the first.
*        the same wsave array can be used by dcfftf and dcfftb.
*
*output parameters
*
*      for j=1,...,n
*
*           c(j)=the sum from k=1,...,n of
*
*                 c(k)*exp(-i*j*k*2*pi/n)
*
*                       where i=sqrt(-1)
*
*wsave   contains initialization calculations which must not be
*        destroyed between calls of subroutine dcfftf or dcfftb
*
*******************************************************************
*
*subroutine dcfftb(n,c,wsave)
*
*******************************************************************
*
*subroutine dcfftb computes the backward complex discrete fourier
*transform (the fourier synthesis). equivalently , dcfftb computes
*a complex periodic sequence from its fourier coefficients.
*the transform is defined below at output parameter c.
*
*a call of dcfftf followed by a call of dcfftb will multiply the
*sequence by n.
*
*the array wsave which is used by subroutine dcfftb must be
*initialized by calling subroutine dcffti(n,wsave).
*
*input parameters
*
*
*n      the length of the complex sequence c. the method is
*       more efficient when n is the product of small primes.
*
*c      a complex array of length n which contains the sequence
*
*wsave   a real work array which must be dimensioned at least 4n+15
*        in the program that calls dcfftb. the wsave array must be
*        initialized by calling subroutine dcffti(n,wsave) and a
*        different wsave array must be used for each different
*        value of n. this initialization does not have to be
*        repeated so long as n remains unchanged thus subsequent
*        transforms can be obtained faster than the first.
*        the same wsave array can be used by dcfftf and dcfftb.
*
*output parameters
*
*      for j=1,...,n
*
*           c(j)=the sum from k=1,...,n of
*
*                 c(k)*exp(i*j*k*2*pi/n)
*
*                       where i=sqrt(-1)
*
*wsave   contains initialization calculations which must not be
*        destroyed between calls of subroutine dcfftf or dcfftb
*************************************************************************

      subroutine dcosqi (n,wsave)
      REAL*8 wsave(*), dt, fk, pih
      data pih /  1.570796326 7948966192 3132169163 975 d0 /
c
      dt = pih/dfloat(n)
      fk = 0.d0
      do 101 k=1,n
         fk = fk+1.d0
         wsave(k) = dcos(fk*dt)
  101 continue
c
      call drffti (n,wsave(n+1))
c
      return
      end
      subroutine dcosti (n,wsave)
      REAL*8 wsave(*), dt, fk, pi
      data pi /  3.141592653 5897932384 6264338327 950 d0 /
c
      if (n .le. 3) return
c
      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      dt = pi/dfloat(nm1)
      fk = 0.d0
      do 101 k=2,ns2
         kc = np1-k
         fk = fk+1.d0
         wsave(k) = 2.d0*dsin(fk*dt)
         wsave(kc) = 2.d0*dcos(fk*dt)
  101 continue
c
      call drffti (nm1,wsave(n+1))
c
      return
      end
      subroutine dcost (n,x,wsave)
      REAL*8 x(*), wsave(*), c1, t1, t2, tx2, x1h, x1p3, xi,
     1  xim2
c
      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      if (n-2) 106,101,102
  101 x1h = x(1)+x(2)
      x(2) = x(1)-x(2)
      x(1) = x1h
      return
c
  102 if (n .gt. 3) go to 103
      x1p3 = x(1)+x(3)
      tx2 = x(2)+x(2)
      x(2) = x(1)-x(3)
      x(1) = x1p3+tx2
      x(3) = x1p3-tx2
      return
c
  103 c1 = x(1)-x(n)
      x(1) = x(1)+x(n)
      do 104 k=2,ns2
         kc = np1-k
         t1 = x(k)+x(kc)
         t2 = x(k)-x(kc)
         c1 = c1+wsave(kc)*t2
         t2 = wsave(k)*t2
         x(k) = t1-t2
         x(kc) = t1+t2
  104 continue
      modn = mod(n,2)
      if (modn .ne. 0) x(ns2+1) = x(ns2+1)+x(ns2+1)
c
      call drfftf (nm1,x,wsave(n+1))
c
      xim2 = x(2)
      x(2) = c1
      do 105 i=4,n,2
         xi = x(i)
         x(i) = x(i-2)-x(i-1)
         x(i-1) = xim2
         xim2 = xi
  105 continue
      if (modn .ne. 0) x(n) = xim2
c
  106 return
      end
      subroutine dradb2 (ido,l1,cc,ch,wa1)
      REAL*8 cc(ido,2,l1), ch(ido,l1,2), wa1(*), ti2, tr2
c
      do 101 k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(ido,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(ido,2,k)
  101 continue
c
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            ch(i-1,k,1) = cc(i-1,1,k)+cc(ic-1,2,k)
            tr2 = cc(i-1,1,k)-cc(ic-1,2,k)
            ch(i,k,1) = cc(i,1,k)-cc(ic,2,k)
            ti2 = cc(i,1,k)+cc(ic,2,k)
            ch(i-1,k,2) = wa1(i-2)*tr2-wa1(i-1)*ti2
            ch(i,k,2) = wa1(i-2)*ti2+wa1(i-1)*tr2
  103    continue
  104 continue
c
      if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(ido,k,1) = cc(ido,1,k)+cc(ido,1,k)
         ch(ido,k,2) = -(cc(1,2,k)+cc(1,2,k))
  106 continue
c
  107 return
      end
      subroutine dradb3 (ido,l1,cc,ch,wa1,wa2)
      REAL*8 cc(ido,3,l1), ch(ido,l1,3), wa1(*), wa2(*),
     1  ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, taui, taur, ti2, tr2
      data taur / -0.5 d0 /
      data taui  /  0.8660254037 8443864676 3723170752 93618d0/
c
      do 101 k=1,l1
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         cr2 = cc(1,1,k)+taur*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ci3 = taui*(cc(1,3,k)+cc(1,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
  101 continue
c
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            cr2 = cc(i-1,1,k)+taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ci2 = cc(i,1,k)+taur*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = taui*(cc(i-1,3,k)-cc(ic-1,2,k))
            ci3 = taui*(cc(i,3,k)+cc(ic,2,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
  102    continue
  103 continue
c
      return
      end
      subroutine dradb4 (ido,l1,cc,ch,wa1,wa2,wa3)
      REAL*8 cc(ido,4,l1), ch(ido,l1,4), wa1(*), wa2(*),
     1  wa3(*), ci2, ci3, ci4, cr2, cr3, cr4, sqrt2, ti1, ti2, ti3, ti4,
     2  tr1, tr2, tr3, tr4
      data sqrt2 /  1.414213562 3730950488 0168872420 970 d0 /
c
      do 101 k=1,l1
         tr1 = cc(1,1,k)-cc(ido,4,k)
         tr2 = cc(1,1,k)+cc(ido,4,k)
         tr3 = cc(ido,2,k)+cc(ido,2,k)
         tr4 = cc(1,3,k)+cc(1,3,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,2) = tr1-tr4
         ch(1,k,3) = tr2-tr3
         ch(1,k,4) = tr1+tr4
  101 continue
c
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            ti1 = cc(i,1,k)+cc(ic,4,k)
            ti2 = cc(i,1,k)-cc(ic,4,k)
            ti3 = cc(i,3,k)-cc(ic,2,k)
            tr4 = cc(i,3,k)+cc(ic,2,k)
            tr1 = cc(i-1,1,k)-cc(ic-1,4,k)
            tr2 = cc(i-1,1,k)+cc(ic-1,4,k)
            ti4 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr3 = cc(i-1,3,k)+cc(ic-1,2,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1-tr4
            cr4 = tr1+tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-2)*cr2-wa1(i-1)*ci2
            ch(i,k,2) = wa1(i-2)*ci2+wa1(i-1)*cr2
            ch(i-1,k,3) = wa2(i-2)*cr3-wa2(i-1)*ci3
            ch(i,k,3) = wa2(i-2)*ci3+wa2(i-1)*cr3
            ch(i-1,k,4) = wa3(i-2)*cr4-wa3(i-1)*ci4
            ch(i,k,4) = wa3(i-2)*ci4+wa3(i-1)*cr4
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
c
  105 continue
      do 106 k=1,l1
         ti1 = cc(1,2,k)+cc(1,4,k)
         ti2 = cc(1,4,k)-cc(1,2,k)
         tr1 = cc(ido,1,k)-cc(ido,3,k)
         tr2 = cc(ido,1,k)+cc(ido,3,k)
         ch(ido,k,1) = tr2+tr2
         ch(ido,k,2) = sqrt2*(tr1-ti1)
         ch(ido,k,3) = ti2+ti2
         ch(ido,k,4) = -sqrt2*(tr1+ti1)
  106 continue
c
  107 return
      end
      subroutine dradb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
      REAL*8 cc(ido,5,l1), ch(ido,l1,5), wa1(*), wa2(*),
     1  wa3(*), wa4(*), ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5,
     2  di2, di3, di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2, ti3,
     3  ti4, ti5, tr11, tr12, tr2, tr3, tr4, tr5
      data tr11  /  0.3090169943 7494742410 2293417182 81906d0/
      data ti11  /  0.9510565162 9515357211 6439333379 38214d0/
      data tr12  / -0.8090169943 7494742410 2293417182 81906d0/
      data ti12  /  0.5877852522 9247312916 8705954639 07277d0/
c
      do 101 k=1,l1
         ti5 = cc(1,3,k)+cc(1,3,k)
         ti4 = cc(1,5,k)+cc(1,5,k)
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         tr3 = cc(ido,4,k)+cc(ido,4,k)
         ch(1,k,1) = cc(1,1,k)+tr2+tr3
         cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
         cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
         ci5 = ti11*ti5+ti12*ti4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2) = cr2-ci5
         ch(1,k,3) = cr3-ci4
         ch(1,k,4) = cr3+ci4
         ch(1,k,5) = cr2+ci5
  101 continue
      if (ido .eq. 1) return
c
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            ti5 = cc(i,3,k)+cc(ic,2,k)
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ti4 = cc(i,5,k)+cc(ic,4,k)
            ti3 = cc(i,5,k)-cc(ic,4,k)
            tr5 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            tr4 = cc(i-1,5,k)-cc(ic-1,4,k)
            tr3 = cc(i-1,5,k)+cc(ic-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
            ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
            cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
            ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
            ch(i-1,k,4) = wa3(i-2)*dr4-wa3(i-1)*di4
            ch(i,k,4) = wa3(i-2)*di4+wa3(i-1)*dr4
            ch(i-1,k,5) = wa4(i-2)*dr5-wa4(i-1)*di5
            ch(i,k,5) = wa4(i-2)*di5+wa4(i-1)*dr5
  102    continue
  103 continue
c
      return
      end
      subroutine dradbg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
      REAL*8 cc(ido,ip,l1), c1(ido,l1,ip), c2(idl1,ip),
     1  ch(ido,l1,ip), ch2(idl1,ip), wa(*), ai1, ai2, ar1, ar1h, ar2,
     2  ar2h, arg, dc2, dcp, ds2, dsp, tpi
      data tpi   /  6.2831853071 7958647692 5286766559 00577d0/
c
      arg = tpi/dfloat(ip)
      dcp = dcos(arg)
      dsp = dsin(arg)
      idp2 = ido+2
      nbd = (ido-1)/2
      ipp2 = ip+2
      ipph = (ip+1)/2
      if (ido .lt. l1) go to 103
      do 102 k=1,l1
         do 101 i=1,ido
            ch(i,k,1) = cc(i,1,k)
  101    continue
  102 continue
      go to 106
c
  103 do 105 i=1,ido
         do 104 k=1,l1
            ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
c
  106 do 108 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 107 k=1,l1
            ch(1,k,j) = cc(ido,j2-2,k)+cc(ido,j2-2,k)
            ch(1,k,jc) = cc(1,j2-1,k)+cc(1,j2-1,k)
  107    continue
  108 continue
c
      if (ido .eq. 1) go to 116
      if (nbd .lt. l1) go to 112
      do 111 j=2,ipph
         jc = ipp2-j
         do 110 k=1,l1
            do 109 i=3,ido,2
               ic = idp2-i
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
  109       continue
  110    continue
  111 continue
      go to 116
c
  112 do 115 j=2,ipph
         jc = ipp2-j
         do 114 i=3,ido,2
            ic = idp2-i
            do 113 k=1,l1
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
  113       continue
  114    continue
  115 continue
c
  116 ar1 = 1.
      ai1 = 0.
      do 120 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 117 ik=1,idl1
            c2(ik,l) = ch2(ik,1)+ar1*ch2(ik,2)
            c2(ik,lc) = ai1*ch2(ik,ip)
  117    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 119 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 118 ik=1,idl1
               c2(ik,l) = c2(ik,l)+ar2*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)+ai2*ch2(ik,jc)
  118       continue
  119    continue
  120 continue
c
      do 122 j=2,ipph
         do 121 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  121    continue
  122 continue
c
      do 124 j=2,ipph
         jc = ipp2-j
         do 123 k=1,l1
            ch(1,k,j) = c1(1,k,j)-c1(1,k,jc)
            ch(1,k,jc) = c1(1,k,j)+c1(1,k,jc)
  123    continue
  124 continue
c
      if (ido .eq. 1) go to 132
      if (nbd .lt. l1) go to 128
      do 127 j=2,ipph
         jc = ipp2-j
         do 126 k=1,l1
            do 125 i=3,ido,2
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
  125       continue
  126    continue
  127 continue
      go to 132
c
  128 do 131 j=2,ipph
         jc = ipp2-j
         do 130 i=3,ido,2
            do 129 k=1,l1
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
  129       continue
  130    continue
  131 continue
  132 continue
c
      if (ido .eq. 1) return
      do 133 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  133 continue
c
      do 135 j=2,ip
         do 134 k=1,l1
            c1(1,k,j) = ch(1,k,j)
  134    continue
  135 continue
c
      if (nbd .gt. l1) go to 139
      is = -ido
      do 138 j=2,ip
         is = is+ido
         idij = is
         do 137 i=3,ido,2
            idij = idij+2
            do 136 k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  136       continue
  137    continue
  138 continue
      go to 143
c
  139 is = -ido
      do 142 j=2,ip
         is = is+ido
         do 141 k=1,l1
            idij = is
            do 140 i=3,ido,2
               idij = idij+2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  140       continue
  141    continue
  142 continue
c
  143 return
      end
      subroutine dradf2 (ido,l1,cc,ch,wa1)
      REAL*8 cc(ido,l1,2), ch(ido,2,l1), wa1(*), ti2, tr2
c
      do 101 k=1,l1
         ch(1,1,k) = cc(1,k,1)+cc(1,k,2)
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,2)
  101 continue
c
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            tr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ti2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            ch(i,1,k) = cc(i,k,1)+ti2
            ch(ic,2,k) = ti2-cc(i,k,1)
            ch(i-1,1,k) = cc(i-1,k,1)+tr2
            ch(ic-1,2,k) = cc(i-1,k,1)-tr2
  103    continue
  104 continue
c
      if (mod(ido,2) .eq. 1) return
  105 do 106 k=1,l1
         ch(1,2,k) = -cc(ido,k,2)
         ch(ido,1,k) = cc(ido,k,1)
  106 continue
c
  107 return
      end
      subroutine dradf3 (ido,l1,cc,ch,wa1,wa2)
      REAL*8 cc(ido,l1,3), ch(ido,3,l1), wa1(*), wa2(*),
     1  ci2, cr2, di2, di3, dr2, dr3, taui, taur, ti2, ti3, tr2, tr3
      data taur / -0.5 d0 /
      data taui  /  0.8660254037 8443864676 3723170752 93618d0/
c
      do 101 k=1,l1
         cr2 = cc(1,k,2)+cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2
         ch(1,3,k) = taui*(cc(1,k,3)-cc(1,k,2))
         ch(ido,2,k) = cc(1,k,1)+taur*cr2
  101 continue
c
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr2 = dr2+dr3
            ci2 = di2+di3
            ch(i-1,1,k) = cc(i-1,k,1)+cr2
            ch(i,1,k) = cc(i,k,1)+ci2
            tr2 = cc(i-1,k,1)+taur*cr2
            ti2 = cc(i,k,1)+taur*ci2
            tr3 = taui*(di2-di3)
            ti3 = taui*(dr3-dr2)
            ch(i-1,3,k) = tr2+tr3
            ch(ic-1,2,k) = tr2-tr3
            ch(i,3,k) = ti2+ti3
            ch(ic,2,k) = ti3-ti2
  102    continue
  103 continue
c
      return
      end
      subroutine dradf4 (ido,l1,cc,ch,wa1,wa2,wa3)
      REAL*8 cc(ido,l1,4), ch(ido,4,l1), wa1(*), wa2(*),
     1  wa3(*), ci2, ci3, ci4, cr2, cr3, cr4, hsqt2, ti1, ti2, ti3,
     2  ti4, tr1, tr2, tr3, tr4
      data hsqt2 /   .7071067811 8654752440 0844362104 85 d0 /
c
      do 101 k=1,l1
         tr1 = cc(1,k,2)+cc(1,k,4)
         tr2 = cc(1,k,1)+cc(1,k,3)
         ch(1,1,k) = tr1+tr2
         ch(ido,4,k) = tr2-tr1
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,3)
         ch(1,3,k) = cc(1,k,4)-cc(1,k,2)
  101 continue
c
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      do 104 k=1,l1
         do 103 i=3,ido,2
            ic = idp2-i
            cr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ci2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            cr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            ci3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            ci4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            tr1 = cr2+cr4
            tr4 = cr4-cr2
            ti1 = ci2+ci4
            ti4 = ci2-ci4
            ti2 = cc(i,k,1)+ci3
            ti3 = cc(i,k,1)-ci3
            tr2 = cc(i-1,k,1)+cr3
            tr3 = cc(i-1,k,1)-cr3
            ch(i-1,1,k) = tr1+tr2
            ch(ic-1,4,k) = tr2-tr1
            ch(i,1,k) = ti1+ti2
            ch(ic,4,k) = ti1-ti2
            ch(i-1,3,k) = ti4+tr3
            ch(ic-1,2,k) = tr3-ti4
            ch(i,3,k) = tr4+ti3
            ch(ic,2,k) = tr4-ti3
  103    continue
  104 continue
      if (mod(ido,2) .eq. 1) return
  105 continue
c
      do 106 k=1,l1
         ti1 = -hsqt2*(cc(ido,k,2)+cc(ido,k,4))
         tr1 = hsqt2*(cc(ido,k,2)-cc(ido,k,4))
         ch(ido,1,k) = tr1+cc(ido,k,1)
         ch(ido,3,k) = cc(ido,k,1)-tr1
         ch(1,2,k) = ti1-cc(ido,k,3)
         ch(1,4,k) = ti1+cc(ido,k,3)
  106 continue
c
  107 return
      end
      subroutine dradf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
      REAL*8 cc(ido,l1,5), ch(ido,5,l1), wa1(*), wa2(*),
     1  wa3(*), wa4(*), ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2,
     2  di3, di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2, ti3, ti4,
     3  ti5, tr11, tr12, tr2, tr3, tr4, tr5
      data tr11  /  0.3090169943 7494742410 2293417182 81906d0/
      data ti11  /  0.9510565162 9515357211 6439333379 38214d0/
      data tr12  / -0.8090169943 7494742410 2293417182 81906d0/
      data ti12  /  0.5877852522 9247312916 8705954639 07277d0/
c
      do 101 k=1,l1
         cr2 = cc(1,k,5)+cc(1,k,2)
         ci5 = cc(1,k,5)-cc(1,k,2)
         cr3 = cc(1,k,4)+cc(1,k,3)
         ci4 = cc(1,k,4)-cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2+cr3
         ch(ido,2,k) = cc(1,k,1)+tr11*cr2+tr12*cr3
         ch(1,3,k) = ti11*ci5+ti12*ci4
         ch(ido,4,k) = cc(1,k,1)+tr12*cr2+tr11*cr3
         ch(1,5,k) = ti12*ci5-ti11*ci4
  101 continue
c
      if (ido .eq. 1) return
      idp2 = ido+2
      do 103 k=1,l1
         do 102 i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            dr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            di4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            dr5 = wa4(i-2)*cc(i-1,k,5)+wa4(i-1)*cc(i,k,5)
            di5 = wa4(i-2)*cc(i,k,5)-wa4(i-1)*cc(i-1,k,5)
            cr2 = dr2+dr5
            ci5 = dr5-dr2
            cr5 = di2-di5
            ci2 = di2+di5
            cr3 = dr3+dr4
            ci4 = dr4-dr3
            cr4 = di3-di4
            ci3 = di3+di4
            ch(i-1,1,k) = cc(i-1,k,1)+cr2+cr3
            ch(i,1,k) = cc(i,k,1)+ci2+ci3
            tr2 = cc(i-1,k,1)+tr11*cr2+tr12*cr3
            ti2 = cc(i,k,1)+tr11*ci2+tr12*ci3
            tr3 = cc(i-1,k,1)+tr12*cr2+tr11*cr3
            ti3 = cc(i,k,1)+tr12*ci2+tr11*ci3
            tr5 = ti11*cr5+ti12*cr4
            ti5 = ti11*ci5+ti12*ci4
            tr4 = ti12*cr5-ti11*cr4
            ti4 = ti12*ci5-ti11*ci4
            ch(i-1,3,k) = tr2+tr5
            ch(ic-1,2,k) = tr2-tr5
            ch(i,3,k) = ti2+ti5
            ch(ic,2,k) = ti5-ti2
            ch(i-1,5,k) = tr3+tr4
            ch(ic-1,4,k) = tr3-tr4
            ch(i,5,k) = ti3+ti4
            ch(ic,4,k) = ti4-ti3
  102    continue
  103 continue
c
      return
      end
      subroutine dradfg (ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
      REAL*8 cc(ido,ip,l1), c1(ido,l1,ip), c2(idl1,ip),
     1  ch(ido,l1,ip), ch2(idl1,ip), wa(*), ai1, ai2, ar1, ar1h, ar2,
     2  ar2h, arg, dc2, dcp, ds2, dsp, tpi
      data tpi   /  6.2831853071 7958647692 5286766559 00577d0/
c
      arg = tpi/dfloat(ip)
      dcp = dcos(arg)
      dsp = dsin(arg)
      ipph = (ip+1)/2
      ipp2 = ip+2
      idp2 = ido+2
      nbd = (ido-1)/2
      if (ido .eq. 1) go to 119
      do 101 ik=1,idl1
         ch2(ik,1) = c2(ik,1)
  101 continue
      do 103 j=2,ip
         do 102 k=1,l1
            ch(1,k,j) = c1(1,k,j)
  102    continue
  103 continue
c
      if (nbd .gt. l1) go to 107
      is = -ido
      do 106 j=2,ip
         is = is+ido
         idij = is
         do 105 i=3,ido,2
            idij = idij+2
            do 104 k=1,l1
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
  104       continue
  105    continue
  106 continue
      go to 111
c
  107 is = -ido
      do 110 j=2,ip
         is = is+ido
         do 109 k=1,l1
            idij = is
            do 108 i=3,ido,2
               idij = idij+2
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
  108       continue
  109    continue
  110 continue
c
  111 if (nbd .lt. l1) go to 115
      do 114 j=2,ipph
         jc = ipp2-j
         do 113 k=1,l1
            do 112 i=3,ido,2
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
  112       continue
  113    continue
  114 continue
      go to 121
c
  115 do 118 j=2,ipph
         jc = ipp2-j
         do 117 i=3,ido,2
            do 116 k=1,l1
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
  116       continue
  117    continue
  118 continue
      go to 121
c
  119 do 120 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  120 continue
c
  121 do 123 j=2,ipph
         jc = ipp2-j
         do 122 k=1,l1
            c1(1,k,j) = ch(1,k,j)+ch(1,k,jc)
            c1(1,k,jc) = ch(1,k,jc)-ch(1,k,j)
  122    continue
  123 continue
c
      ar1 = 1.d0
      ai1 = 0.d0
      do 127 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 124 ik=1,idl1
            ch2(ik,l) = c2(ik,1)+ar1*c2(ik,2)
            ch2(ik,lc) = ai1*c2(ik,ip)
  124    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 126 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 125 ik=1,idl1
               ch2(ik,l) = ch2(ik,l)+ar2*c2(ik,j)
               ch2(ik,lc) = ch2(ik,lc)+ai2*c2(ik,jc)
  125       continue
  126    continue
  127 continue
c
      do 129 j=2,ipph
         do 128 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+c2(ik,j)
  128    continue
  129 continue
c
      if (ido .lt. l1) go to 132
      do 131 k=1,l1
         do 130 i=1,ido
            cc(i,1,k) = ch(i,k,1)
  130    continue
  131 continue
      go to 135
c
  132 do 134 i=1,ido
         do 133 k=1,l1
            cc(i,1,k) = ch(i,k,1)
  133    continue
  134 continue
c
  135 do 137 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 136 k=1,l1
            cc(ido,j2-2,k) = ch(1,k,j)
            cc(1,j2-1,k) = ch(1,k,jc)
  136    continue
  137 continue
c
      if (ido .eq. 1) return
      if (nbd .lt. l1) go to 141
      do 140 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 139 k=1,l1
            do 138 i=3,ido,2
               ic = idp2-i
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
  138       continue
  139    continue
  140 continue
      return
c
  141 do 144 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 143 i=3,ido,2
            ic = idp2-i
            do 142 k=1,l1
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
  142       continue
  143    continue
  144 continue
c
      return
      end
      subroutine drfftb (n,r,wsave)
      REAL*8 r(*), wsave(*)
c
      if (n .eq. 1) return
c
      call drftb1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
c
      return
      end
      subroutine drfftf (n,r,wsave)
      REAL*8 r(*), wsave(*)
c
      if (n .eq. 1) return
c
      call drftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
c
      return
      end
      subroutine drffti (n,wsave)
      REAL*8 wsave(*)
c
      if (n .eq. 1) return
c
      call drfti1 (n,wsave(n+1),wsave(2*n+1))
c
      return
      end
      subroutine drftb1 (n,c,ch,wa,ifac)
      REAL*8 c(*), ch(*), wa(*)
      integer ifac(*)
c
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call dradb4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call dradb4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
c
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call dradb2 (ido,l1,c,ch,wa(iw))
         go to 105
  104    call dradb2 (ido,l1,ch,c,wa(iw))
  105    na = 1-na
         go to 115
c
  106    if (ip .ne. 3) go to 109
         ix2 = iw+ido
         if (na .ne. 0) go to 107
         call dradb3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 108
  107    call dradb3 (ido,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
c
  109    if (ip .ne. 5) go to 112
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 110
         call dradb5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call dradb5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
c
  112    if (na .ne. 0) go to 113
         call dradbg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         go to 114
  113    call dradbg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (ido .eq. 1) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*ido
  116 continue
c
      if (na .eq. 0) return
      do 117 i=1,n
         c(i) = ch(i)
  117 continue
c
      return
      end
      subroutine drftf1 (n,c,ch,wa,ifac)
      REAL*8 c(*), ch(*), wa(*)
      integer ifac(*)
c
      nf = ifac(2)
      na = 1
      l2 = n
      iw = n
      do 111 k1=1,nf
         kh = nf-k1
         ip = ifac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
         if (ip .ne. 4) go to 102
c
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call dradf4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 110
  101    call dradf4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
         go to 110
c
  102    if (ip .ne. 2) go to 104
         if (na .ne. 0) go to 103
         call dradf2 (ido,l1,c,ch,wa(iw))
         go to 110
  103    call dradf2 (ido,l1,ch,c,wa(iw))
         go to 110
c
  104    if (ip .ne. 3) go to 106
         ix2 = iw+ido
         if (na .ne. 0) go to 105
         call dradf3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 110
  105    call dradf3 (ido,l1,ch,c,wa(iw),wa(ix2))
         go to 110
c
  106    if (ip .ne. 5) go to 108
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 107
         call dradf5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  107    call dradf5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
c
  108    if (ido .eq. 1) na = 1-na
         if (na .ne. 0) go to 109
         call dradfg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         na = 1
         go to 110
  109    call dradfg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
         na = 0
c
  110    l2 = l1
  111 continue
c
      if (na .eq. 1) return
      do 112 i=1,n
         c(i) = ch(i)
  112 continue
c
      return
      end
      subroutine drfti1 (n,wa,ifac)
      REAL*8 wa(*), arg, argh, argld, fi, tpi
      integer ifac(*), ntryh(4)
      data ntryh(1), ntryh(2), ntryh(3), ntryh(4) /4, 2, 3, 5/
      data tpi   /  6.2831853071 7958647692 5286766559 00577d0/
c
      nl = n
      nf = 0
      j = 0
c
  101 j = j+1
      if (j.le.4) ntry = ntryh(j)
      if (j.gt.4) ntry = ntry + 2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr.ne.0) go to 101
c
  105 nf = nf+1
      ifac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         ifac(ib+2) = ifac(ib+1)
  106 continue
      ifac(3) = 2
  107 if (nl .ne. 1) go to 104
      ifac(1) = n
      ifac(2) = nf
c
      argh = tpi/dfloat(n)
      is = 0
      nfm1 = nf-1
      l1 = 1
      if (nfm1 .eq. 0) return
      do 110 k1=1,nfm1
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         do 109 j=1,ipm
            ld = ld+l1
            i = is
            argld = dfloat(ld)*argh
            fi = 0.d0
            do 108 ii=3,ido,2
               i = i+2
               fi = fi+1.d0
               arg = fi*argld
               wa(i-1) = dcos(arg)
               wa(i) = dsin(arg)
  108       continue
            is = is+ido
  109    continue
c
         l1 = l2
  110 continue
c
      return
      end

      subroutine dcffti (n,wsave)
      REAL*8 wsave(*)
c
      if (n .eq. 1) return
c
      iw1 = n+n+1
      iw2 = iw1+n+n
      call dcfti1 (n,wsave(iw1),wsave(iw2))
c
      return
      end
      subroutine dcfti1 (n,wa,ifac)
      REAL*8 wa(*), arg, argh, argld, fi, tpi
      integer ifac(*), ntryh(4)
      data ntryh(1), ntryh(2), ntryh(3), ntryh(4) /3, 4, 2, 5/
      data tpi   /  6.2831853071 7958647692 5286766559 00577d0/
c
      nl = n
      nf = 0
      j = 0
c
  101 j = j+1
      if (j.le.4) ntry = ntryh(j)
      if (j.gt.4) ntry = ntry + 2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr.ne.0) go to 101
c
  105 nf = nf+1
      ifac(nf+2) = ntry
      nl = nq
      if (ntry .ne. 2) go to 107
      if (nf .eq. 1) go to 107
      do 106 i=2,nf
         ib = nf-i+2
         ifac(ib+2) = ifac(ib+1)
  106 continue
      ifac(3) = 2
c
  107 if (nl .ne. 1) go to 104
c
      ifac(1) = n
      ifac(2) = nf
c
      argh = tpi/dfloat(n)
      i = 2
      l1 = 1
      do 110 k1=1,nf
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         idot = ido+ido+2
         ipm = ip-1
c
         do 109 j=1,ipm
            i1 = i
            wa(i-1) = 1.d0
            wa(i) = 0.d0
            ld = ld+l1
            fi = 0.d0
            argld = dfloat(ld)*argh
            do 108 ii=4,idot,2
               i = i+2
               fi = fi+1.d0
               arg = fi*argld
               wa(i-1) = dcos(arg)
               wa(i) = dsin(arg)
  108       continue
            if (ip .le. 5) go to 109
            wa(i1-1) = wa(i-1)
            wa(i1) = wa(i)
  109    continue
c
         l1 = l2
  110 continue
c
      return
      end

      subroutine dcfftf (n,c,wsave)
      REAL*8 c(*), wsave(*)
c
      if (n .eq. 1) return
c
      iw1 = n+n+1
      iw2 = iw1+n+n
      call dcftf1 (n,c,wsave,wsave(iw1),wsave(iw2))
c
      return
      end
      subroutine dcftf1 (n,c,ch,wa,ifac)
      REAL*8 c(*), ch(*), wa(*)
      integer ifac(*)
c
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido+ido
         idl1 = idot*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+idot
         ix3 = ix2+idot
         if (na .ne. 0) go to 101
         call dpssf4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call dpssf4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
c
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call dpssf2 (idot,l1,c,ch,wa(iw))
         go to 105
  104    call dpssf2 (idot,l1,ch,c,wa(iw))
  105    na = 1-na
         go to 115
c
  106    if (ip .ne. 3) go to 109
         ix2 = iw+idot
         if (na .ne. 0) go to 107
         call dpssf3 (idot,l1,c,ch,wa(iw),wa(ix2))
         go to 108
  107    call dpssf3 (idot,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
c
  109    if (ip .ne. 5) go to 112
         ix2 = iw+idot
         ix3 = ix2+idot
         ix4 = ix3+idot
         if (na .ne. 0) go to 110
         call dpssf5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call dpssf5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
c
  112    if (na .ne. 0) go to 113
         call dpssf (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         go to 114
  113    call dpssf (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (nac .ne. 0) na = 1-na
c
  115    l1 = l2
         iw = iw+(ip-1)*idot
  116 continue
      if (na .eq. 0) return
c
      n2 = n+n
      do 117 i=1,n2
         c(i) = ch(i)
  117 continue
c
      return
      end
      subroutine dpssf (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
      REAL*8 cc(ido,ip,l1), c1(ido,l1,ip), c2(idl1,ip),
     1  ch(ido,l1,ip), ch2(idl1,ip), wa(*), wai, war
c
      idot = ido/2
      nt = ip*idl1
      ipp2 = ip+2
      ipph = (ip+1)/2
      idp = ip*ido
c
      if (ido .lt. l1) go to 106
      do 103 j=2,ipph
         jc = ipp2-j
         do 102 k=1,l1
            do 101 i=1,ido
               ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  101       continue
  102    continue
  103 continue
c
      do 105 k=1,l1
         do 104 i=1,ido
            ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
      go to 112
c
  106 do 109 j=2,ipph
         jc = ipp2-j
         do 108 i=1,ido
            do 107 k=1,l1
               ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  107       continue
  108    continue
  109 continue
c
      do 111 i=1,ido
         do 110 k=1,l1
            ch(i,k,1) = cc(i,1,k)
  110    continue
  111 continue
c
  112 idl = 2-ido
      inc = 0
      do 116 l=2,ipph
         lc = ipp2-l
         idl = idl+ido
         do 113 ik=1,idl1
            c2(ik,l) = ch2(ik,1)+wa(idl-1)*ch2(ik,2)
            c2(ik,lc) = -wa(idl)*ch2(ik,ip)
  113    continue
         idlj = idl
         inc = inc+ido
         do 115 j=3,ipph
            jc = ipp2-j
            idlj = idlj+inc
            if (idlj .gt. idp) idlj = idlj-idp
            war = wa(idlj-1)
            wai = wa(idlj)
            do 114 ik=1,idl1
               c2(ik,l) = c2(ik,l)+war*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)-wai*ch2(ik,jc)
  114       continue
  115    continue
  116 continue
c
      do 118 j=2,ipph
         do 117 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  117    continue
  118 continue
c
      do 120 j=2,ipph
         jc = ipp2-j
         do 119 ik=2,idl1,2
            ch2(ik-1,j) = c2(ik-1,j)-c2(ik,jc)
            ch2(ik-1,jc) = c2(ik-1,j)+c2(ik,jc)
            ch2(ik,j) = c2(ik,j)+c2(ik-1,jc)
            ch2(ik,jc) = c2(ik,j)-c2(ik-1,jc)
  119    continue
  120 continue
c
      nac = 1
      if (ido .eq. 2) return
      nac = 0
c
      do 121 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  121 continue
c
      do 123 j=2,ip
         do 122 k=1,l1
            c1(1,k,j) = ch(1,k,j)
            c1(2,k,j) = ch(2,k,j)
  122    continue
  123 continue
c
      if (idot .gt. l1) go to 127
      idij = 0
      do 126 j=2,ip
         idij = idij+2
         do 125 i=4,ido,2
            idij = idij+2
            do 124 k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)+wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)-wa(idij)*ch(i-1,k,j)
  124       continue
  125    continue
  126 continue
      return
c
  127 idj = 2-ido
      do 130 j=2,ip
         idj = idj+ido
         do 129 k=1,l1
            idij = idj
            do 128 i=4,ido,2
               idij = idij+2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)+wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)-wa(idij)*ch(i-1,k,j)
  128       continue
  129    continue
  130 continue
c
      return
      end
      subroutine dpssf2 (ido,l1,cc,ch,wa1)
      REAL*8 cc(ido,2,l1), ch(ido,l1,2), wa1(*), ti2, tr2
c
      if (ido .gt. 2) go to 102
      do 101 k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(1,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(1,2,k)
         ch(2,k,1) = cc(2,1,k)+cc(2,2,k)
         ch(2,k,2) = cc(2,1,k)-cc(2,2,k)
  101 continue
      return
c
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ch(i-1,k,1) = cc(i-1,1,k)+cc(i-1,2,k)
            tr2 = cc(i-1,1,k)-cc(i-1,2,k)
            ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
            ti2 = cc(i,1,k)-cc(i,2,k)
            ch(i,k,2) = wa1(i-1)*ti2-wa1(i)*tr2
            ch(i-1,k,2) = wa1(i-1)*tr2+wa1(i)*ti2
  103    continue
  104 continue
c
      return
      end
      subroutine dpssf3 (ido,l1,cc,ch,wa1,wa2)
      REAL*8 cc(ido,3,l1), ch(ido,l1,3), wa1(*), wa2(*),
     1  ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, taui, taur, ti2, tr2
      data taur / -0.5 d0 /
      data taui  / -0.8660254037 8443864676 3723170752 93618d0/
c
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         tr2 = cc(1,2,k)+cc(1,3,k)
         cr2 = cc(1,1,k)+taur*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ti2 = cc(2,2,k)+cc(2,3,k)
         ci2 = cc(2,1,k)+taur*ti2
         ch(2,k,1) = cc(2,1,k)+ti2
         cr3 = taui*(cc(1,2,k)-cc(1,3,k))
         ci3 = taui*(cc(2,2,k)-cc(2,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
         ch(2,k,2) = ci2+cr3
         ch(2,k,3) = ci2-cr3
  101 continue
      return
c
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            tr2 = cc(i-1,2,k)+cc(i-1,3,k)
            cr2 = cc(i-1,1,k)+taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,2,k)+cc(i,3,k)
            ci2 = cc(i,1,k)+taur*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = taui*(cc(i-1,2,k)-cc(i-1,3,k))
            ci3 = taui*(cc(i,2,k)-cc(i,3,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i,k,2) = wa1(i-1)*di2-wa1(i)*dr2
            ch(i-1,k,2) = wa1(i-1)*dr2+wa1(i)*di2
            ch(i,k,3) = wa2(i-1)*di3-wa2(i)*dr3
            ch(i-1,k,3) = wa2(i-1)*dr3+wa2(i)*di3
  103    continue
  104 continue
c
      return
      end
      subroutine dpssf4 (ido,l1,cc,ch,wa1,wa2,wa3)
      REAL*8 cc(ido,4,l1), ch(ido,l1,4), wa1(*), wa2(*),
     1  wa3(*), ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4,
     2  tr1, tr2, tr3, tr4
c
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         ti1 = cc(2,1,k)-cc(2,3,k)
         ti2 = cc(2,1,k)+cc(2,3,k)
         tr4 = cc(2,2,k)-cc(2,4,k)
         ti3 = cc(2,2,k)+cc(2,4,k)
         tr1 = cc(1,1,k)-cc(1,3,k)
         tr2 = cc(1,1,k)+cc(1,3,k)
         ti4 = cc(1,4,k)-cc(1,2,k)
         tr3 = cc(1,2,k)+cc(1,4,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,3) = tr2-tr3
         ch(2,k,1) = ti2+ti3
         ch(2,k,3) = ti2-ti3
         ch(1,k,2) = tr1+tr4
         ch(1,k,4) = tr1-tr4
         ch(2,k,2) = ti1+ti4
         ch(2,k,4) = ti1-ti4
  101 continue
      return
c
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ti1 = cc(i,1,k)-cc(i,3,k)
            ti2 = cc(i,1,k)+cc(i,3,k)
            ti3 = cc(i,2,k)+cc(i,4,k)
            tr4 = cc(i,2,k)-cc(i,4,k)
            tr1 = cc(i-1,1,k)-cc(i-1,3,k)
            tr2 = cc(i-1,1,k)+cc(i-1,3,k)
            ti4 = cc(i-1,4,k)-cc(i-1,2,k)
            tr3 = cc(i-1,2,k)+cc(i-1,4,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1+tr4
            cr4 = tr1-tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-1)*cr2+wa1(i)*ci2
            ch(i,k,2) = wa1(i-1)*ci2-wa1(i)*cr2
            ch(i-1,k,3) = wa2(i-1)*cr3+wa2(i)*ci3
            ch(i,k,3) = wa2(i-1)*ci3-wa2(i)*cr3
            ch(i-1,k,4) = wa3(i-1)*cr4+wa3(i)*ci4
            ch(i,k,4) = wa3(i-1)*ci4-wa3(i)*cr4
  103    continue
  104 continue
c
      return
      end
      subroutine dpssf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
      REAL*8 cc(ido,5,l1), ch(ido,l1,5), wa1(*), wa2(*),
     1  wa3(*), wa4(*), ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2,
     2  di3, di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2, ti3, ti4,
     3  ti5, tr11, tr12, tr2, tr3, tr4, tr5
      data tr11  /  0.3090169943 7494742410 2293417182 81906d0/
      data ti11  / -0.9510565162 9515357211 6439333379 38214d0/
      data tr12  / -0.8090169943 7494742410 2293417182 81906d0/
      data ti12  / -0.5877852522 9247312916 8705954639 07277d0/
c
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         ti5 = cc(2,2,k)-cc(2,5,k)
         ti2 = cc(2,2,k)+cc(2,5,k)
         ti4 = cc(2,3,k)-cc(2,4,k)
         ti3 = cc(2,3,k)+cc(2,4,k)
         tr5 = cc(1,2,k)-cc(1,5,k)
         tr2 = cc(1,2,k)+cc(1,5,k)
         tr4 = cc(1,3,k)-cc(1,4,k)
         tr3 = cc(1,3,k)+cc(1,4,k)
         ch(1,k,1) = cc(1,1,k)+tr2+tr3
         ch(2,k,1) = cc(2,1,k)+ti2+ti3
         cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
         ci2 = cc(2,1,k)+tr11*ti2+tr12*ti3
         cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
         ci3 = cc(2,1,k)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2) = cr2-ci5
         ch(1,k,5) = cr2+ci5
         ch(2,k,2) = ci2+cr5
         ch(2,k,3) = ci3+cr4
         ch(1,k,3) = cr3-ci4
         ch(1,k,4) = cr3+ci4
         ch(2,k,4) = ci3-cr4
         ch(2,k,5) = ci2-cr5
  101 continue
      return
c
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ti5 = cc(i,2,k)-cc(i,5,k)
            ti2 = cc(i,2,k)+cc(i,5,k)
            ti4 = cc(i,3,k)-cc(i,4,k)
            ti3 = cc(i,3,k)+cc(i,4,k)
            tr5 = cc(i-1,2,k)-cc(i-1,5,k)
            tr2 = cc(i-1,2,k)+cc(i-1,5,k)
            tr4 = cc(i-1,3,k)-cc(i-1,4,k)
            tr3 = cc(i-1,3,k)+cc(i-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
            ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
            cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
            ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-1)*dr2+wa1(i)*di2
            ch(i,k,2) = wa1(i-1)*di2-wa1(i)*dr2
            ch(i-1,k,3) = wa2(i-1)*dr3+wa2(i)*di3
            ch(i,k,3) = wa2(i-1)*di3-wa2(i)*dr3
            ch(i-1,k,4) = wa3(i-1)*dr4+wa3(i)*di4
            ch(i,k,4) = wa3(i-1)*di4-wa3(i)*dr4
            ch(i-1,k,5) = wa4(i-1)*dr5+wa4(i)*di5
            ch(i,k,5) = wa4(i-1)*di5-wa4(i)*dr5
  103    continue
  104 continue
c
      return
      end

      subroutine dcfftb (n,c,wsave)
      REAL*8 c(*), wsave(*)
c
      if (n .eq. 1) return
c
      iw1 = n+n+1
      iw2 = iw1+n+n
      call dcftb1 (n,c,wsave,wsave(iw1),wsave(iw2))
c
      return
      end
      subroutine dcftb1 (n,c,ch,wa,ifac)
      REAL*8 c(*), ch(*), wa(*)
      integer ifac(*)
c
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idot = ido+ido
         idl1 = idot*l1
         if (ip .ne. 4) go to 103
         ix2 = iw+idot
         ix3 = ix2+idot
         if (na .ne. 0) go to 101
         call dpssb4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call dpssb4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
c
  103    if (ip .ne. 2) go to 106
         if (na .ne. 0) go to 104
         call dpssb2 (idot,l1,c,ch,wa(iw))
         go to 105
  104    call dpssb2 (idot,l1,ch,c,wa(iw))
  105    na = 1-na
         go to 115
c
  106    if (ip .ne. 3) go to 109
         ix2 = iw+idot
         if (na .ne. 0) go to 107
         call dpssb3 (idot,l1,c,ch,wa(iw),wa(ix2))
         go to 108
  107    call dpssb3 (idot,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
c
  109    if (ip .ne. 5) go to 112
         ix2 = iw+idot
         ix3 = ix2+idot
         ix4 = ix3+idot
         if (na .ne. 0) go to 110
         call dpssb5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call dpssb5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
c
  112    if (na .ne. 0) go to 113
         call dpssb (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         go to 114
  113    call dpssb (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (nac .ne. 0) na = 1-na
c
  115    l1 = l2
         iw = iw+(ip-1)*idot
  116 continue
      if (na .eq. 0) return
c
      n2 = n+n
      do 117 i=1,n2
         c(i) = ch(i)
  117 continue
c
      return
      end
      subroutine dpssb (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
      REAL*8 cc(ido,ip,l1), c1(ido,l1,ip), c2(idl1,ip),
     1  ch(ido,l1,ip), ch2(idl1,ip), wa(*), wai, war
c
      idot = ido/2
      nt = ip*idl1
      ipp2 = ip+2
      ipph = (ip+1)/2
      idp = ip*ido
c
      if (ido .lt. l1) go to 106
      do 103 j=2,ipph
         jc = ipp2-j
         do 102 k=1,l1
            do 101 i=1,ido
               ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  101       continue
  102    continue
  103 continue
c
      do 105 k=1,l1
         do 104 i=1,ido
            ch(i,k,1) = cc(i,1,k)
  104    continue
  105 continue
      go to 112
c
  106 do 109 j=2,ipph
         jc = ipp2-j
         do 108 i=1,ido
            do 107 k=1,l1
               ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
               ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
  107       continue
  108    continue
  109 continue
c
      do 111 i=1,ido
         do 110 k=1,l1
            ch(i,k,1) = cc(i,1,k)
  110    continue
  111 continue
c
  112 idl = 2-ido
      inc = 0
      do 116 l=2,ipph
         lc = ipp2-l
         idl = idl+ido
         do 113 ik=1,idl1
            c2(ik,l) = ch2(ik,1)+wa(idl-1)*ch2(ik,2)
            c2(ik,lc) = wa(idl)*ch2(ik,ip)
  113    continue
         idlj = idl
         inc = inc+ido
         do 115 j=3,ipph
            jc = ipp2-j
            idlj = idlj+inc
            if (idlj .gt. idp) idlj = idlj-idp
            war = wa(idlj-1)
            wai = wa(idlj)
            do 114 ik=1,idl1
               c2(ik,l) = c2(ik,l)+war*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)+wai*ch2(ik,jc)
  114       continue
  115    continue
  116 continue
c
      do 118 j=2,ipph
         do 117 ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
  117    continue
  118 continue
c
      do 120 j=2,ipph
         jc = ipp2-j
         do 119 ik=2,idl1,2
            ch2(ik-1,j) = c2(ik-1,j)-c2(ik,jc)
            ch2(ik-1,jc) = c2(ik-1,j)+c2(ik,jc)
            ch2(ik,j) = c2(ik,j)+c2(ik-1,jc)
            ch2(ik,jc) = c2(ik,j)-c2(ik-1,jc)
  119    continue
  120 continue
c
      nac = 1
      if (ido .eq. 2) return
      nac = 0
c
      do 121 ik=1,idl1
         c2(ik,1) = ch2(ik,1)
  121 continue
c
      do 123 j=2,ip
         do 122 k=1,l1
            c1(1,k,j) = ch(1,k,j)
            c1(2,k,j) = ch(2,k,j)
  122    continue
  123 continue
c
      if (idot .gt. l1) go to 127
      idij = 0
      do 126 j=2,ip
         idij = idij+2
         do 125 i=4,ido,2
            idij = idij+2
            do 124 k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  124       continue
  125    continue
  126 continue
      return
c
  127 idj = 2-ido
      do 130 j=2,ip
         idj = idj+ido
         do 129 k=1,l1
            idij = idj
            do 128 i=4,ido,2
               idij = idij+2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
  128       continue
  129    continue
  130 continue
c
      return
      end
      subroutine dpssb2 (ido,l1,cc,ch,wa1)
      REAL*8 cc(ido,2,l1), ch(ido,l1,2), wa1(*), ti2, tr2
c
      if (ido .gt. 2) go to 102
      do 101 k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(1,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(1,2,k)
         ch(2,k,1) = cc(2,1,k)+cc(2,2,k)
         ch(2,k,2) = cc(2,1,k)-cc(2,2,k)
  101 continue
      return
c
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ch(i-1,k,1) = cc(i-1,1,k)+cc(i-1,2,k)
            tr2 = cc(i-1,1,k)-cc(i-1,2,k)
            ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
            ti2 = cc(i,1,k)-cc(i,2,k)
            ch(i,k,2) = wa1(i-1)*ti2+wa1(i)*tr2
            ch(i-1,k,2) = wa1(i-1)*tr2-wa1(i)*ti2
  103    continue
  104 continue
c
      return
      end
      subroutine dpssb3 (ido,l1,cc,ch,wa1,wa2)
      REAL*8 cc(ido,3,l1), ch(ido,l1,3), wa1(*), wa2(*),
     1 ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, taui, taur, ti2, tr2
      data taur / -0.5 d0 /
      data taui  /  0.8660254037 8443864676 3723170752 93618d0/
c
c     one half sqrt(3) = .866025.....  .
c
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         tr2 = cc(1,2,k)+cc(1,3,k)
         cr2 = cc(1,1,k)+taur*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ti2 = cc(2,2,k)+cc(2,3,k)
         ci2 = cc(2,1,k)+taur*ti2
         ch(2,k,1) = cc(2,1,k)+ti2
         cr3 = taui*(cc(1,2,k)-cc(1,3,k))
         ci3 = taui*(cc(2,2,k)-cc(2,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
         ch(2,k,2) = ci2+cr3
         ch(2,k,3) = ci2-cr3
  101 continue
      return
c
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            tr2 = cc(i-1,2,k)+cc(i-1,3,k)
            cr2 = cc(i-1,1,k)+taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,2,k)+cc(i,3,k)
            ci2 = cc(i,1,k)+taur*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = taui*(cc(i-1,2,k)-cc(i-1,3,k))
            ci3 = taui*(cc(i,2,k)-cc(i,3,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i,k,2) = wa1(i-1)*di2+wa1(i)*dr2
            ch(i-1,k,2) = wa1(i-1)*dr2-wa1(i)*di2
            ch(i,k,3) = wa2(i-1)*di3+wa2(i)*dr3
            ch(i-1,k,3) = wa2(i-1)*dr3-wa2(i)*di3
  103    continue
  104 continue
c
      return
      end
      subroutine dpssb4 (ido,l1,cc,ch,wa1,wa2,wa3)
      REAL*8 cc(ido,4,l1), ch(ido,l1,4), wa1(*), wa2(*),
     1  wa3(*), ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1,
     2  tr2, tr3, tr4
c
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         ti1 = cc(2,1,k)-cc(2,3,k)
         ti2 = cc(2,1,k)+cc(2,3,k)
         tr4 = cc(2,4,k)-cc(2,2,k)
         ti3 = cc(2,2,k)+cc(2,4,k)
         tr1 = cc(1,1,k)-cc(1,3,k)
         tr2 = cc(1,1,k)+cc(1,3,k)
         ti4 = cc(1,2,k)-cc(1,4,k)
         tr3 = cc(1,2,k)+cc(1,4,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,3) = tr2-tr3
         ch(2,k,1) = ti2+ti3
         ch(2,k,3) = ti2-ti3
         ch(1,k,2) = tr1+tr4
         ch(1,k,4) = tr1-tr4
         ch(2,k,2) = ti1+ti4
         ch(2,k,4) = ti1-ti4
  101 continue
      return
c
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ti1 = cc(i,1,k)-cc(i,3,k)
            ti2 = cc(i,1,k)+cc(i,3,k)
            ti3 = cc(i,2,k)+cc(i,4,k)
            tr4 = cc(i,4,k)-cc(i,2,k)
            tr1 = cc(i-1,1,k)-cc(i-1,3,k)
            tr2 = cc(i-1,1,k)+cc(i-1,3,k)
            ti4 = cc(i-1,2,k)-cc(i-1,4,k)
            tr3 = cc(i-1,2,k)+cc(i-1,4,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1+tr4
            cr4 = tr1-tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-1)*cr2-wa1(i)*ci2
            ch(i,k,2) = wa1(i-1)*ci2+wa1(i)*cr2
            ch(i-1,k,3) = wa2(i-1)*cr3-wa2(i)*ci3
            ch(i,k,3) = wa2(i-1)*ci3+wa2(i)*cr3
            ch(i-1,k,4) = wa3(i-1)*cr4-wa3(i)*ci4
            ch(i,k,4) = wa3(i-1)*ci4+wa3(i)*cr4
  103    continue
  104 continue
c
      return
      end
      subroutine dpssb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
      REAL*8 cc(ido,5,l1), ch(ido,l1,5), wa1(*), wa2(*),
     1  wa3(*), wa4(*), ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5,
     2  di2, di3, di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2, ti3,
     3  ti4, ti5, tr11, tr12, tr2, tr3, tr4, tr5
      data tr11  /  0.3090169943 7494742410 2293417182 81906d0/
      data ti11  /  0.9510565162 9515357211 6439333379 38214d0/
      data tr12  / -0.8090169943 7494742410 2293417182 81906d0/
      data ti12  /  0.5877852522 9247312916 8705954639 07277d0/
c
c     sin(pi/10) = .30901699....    .
c     cos(pi/10) = .95105651....    .
c     sin(pi/5 ) = .58778525....    .
c     cos(pi/5 ) = .80901699....    .
c
      if (ido .ne. 2) go to 102
      do 101 k=1,l1
         ti5 = cc(2,2,k)-cc(2,5,k)
         ti2 = cc(2,2,k)+cc(2,5,k)
         ti4 = cc(2,3,k)-cc(2,4,k)
         ti3 = cc(2,3,k)+cc(2,4,k)
         tr5 = cc(1,2,k)-cc(1,5,k)
         tr2 = cc(1,2,k)+cc(1,5,k)
         tr4 = cc(1,3,k)-cc(1,4,k)
         tr3 = cc(1,3,k)+cc(1,4,k)
         ch(1,k,1) = cc(1,1,k)+tr2+tr3
         ch(2,k,1) = cc(2,1,k)+ti2+ti3
         cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
         ci2 = cc(2,1,k)+tr11*ti2+tr12*ti3
         cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
         ci3 = cc(2,1,k)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2) = cr2-ci5
         ch(1,k,5) = cr2+ci5
         ch(2,k,2) = ci2+cr5
         ch(2,k,3) = ci3+cr4
         ch(1,k,3) = cr3-ci4
         ch(1,k,4) = cr3+ci4
         ch(2,k,4) = ci3-cr4
         ch(2,k,5) = ci2-cr5
  101 continue
      return
c
  102 do 104 k=1,l1
         do 103 i=2,ido,2
            ti5 = cc(i,2,k)-cc(i,5,k)
            ti2 = cc(i,2,k)+cc(i,5,k)
            ti4 = cc(i,3,k)-cc(i,4,k)
            ti3 = cc(i,3,k)+cc(i,4,k)
            tr5 = cc(i-1,2,k)-cc(i-1,5,k)
            tr2 = cc(i-1,2,k)+cc(i-1,5,k)
            tr4 = cc(i-1,3,k)-cc(i-1,4,k)
            tr3 = cc(i-1,3,k)+cc(i-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
            ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
            cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
            ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-1)*dr2-wa1(i)*di2
            ch(i,k,2) = wa1(i-1)*di2+wa1(i)*dr2
            ch(i-1,k,3) = wa2(i-1)*dr3-wa2(i)*di3
            ch(i,k,3) = wa2(i-1)*di3+wa2(i)*dr3
            ch(i-1,k,4) = wa3(i-1)*dr4-wa3(i)*di4
            ch(i,k,4) = wa3(i-1)*di4+wa3(i)*dr4
            ch(i-1,k,5) = wa4(i-1)*dr5-wa4(i)*di5
            ch(i,k,5) = wa4(i-1)*di5+wa4(i)*dr5
  103    continue
  104 continue
c
      return
      end
      subroutine dcosqb (n,x,wsave)
      REAL*8 x(*), wsave(*), tsqrt2, x1
      data tsqrt2 /  2.828427124 7461900976 0337744841 94 d0 /
c
      if (n-2) 101,102,103
  101 x(1) = 4.d0*x(1)
      return
c
  102 x1 = 4.d0*(x(1)+x(2))
      x(2) = tsqrt2*(x(1)-x(2))
      x(1) = x1
      return
c
  103 call dcsqb1 (n,x,wsave,wsave(n+1))
c
      return
      end
      subroutine dcosqf (n,x,wsave)
      REAL*8 x(*), wsave(*), sqrt2, tsqx
      data sqrt2 /  1.414213562 3730950488 0168872420 970 d0 /
c
      if (n-2) 102,101,103
  101 tsqx = sqrt2*x(2)
      x(2) = x(1)-tsqx
      x(1) = x(1)+tsqx
  102 return
c
  103 call dcsqf1 (n,x,wsave,wsave(n+1))
c
      return
      end
      subroutine dcsqf1 (n,x,w,xh)
      REAL*8 x(*), w(*), xh(*), xim1
c
      ns2 = (n+1)/2
      np2 = n+2
      do 101 k=2,ns2
         kc = np2-k
         xh(k) = x(k)+x(kc)
         xh(kc) = x(k)-x(kc)
  101 continue
      modn = mod(n,2)
      if (modn .eq. 0) xh(ns2+1) = x(ns2+1)+x(ns2+1)
c
      do 102 k=2,ns2
         kc = np2-k
         x(k) = w(k-1)*xh(kc)+w(kc-1)*xh(k)
         x(kc) = w(k-1)*xh(k)-w(kc-1)*xh(kc)
  102 continue
      if (modn .eq. 0) x(ns2+1) = w(ns2)*xh(ns2+1)
c
      call drfftf (n,x,xh)
c
      do 103 i=3,n,2
         xim1 = x(i-1)-x(i)
         x(i) = x(i-1)+x(i)
         x(i-1) = xim1
  103 continue
c
      return
      end
      subroutine dcsqb1 (n,x,w,xh)
      REAL*8 x(*), w(*), xh(*), xim1
c
      ns2 = (n+1)/2
      np2 = n+2
      do 101 i=3,n,2
         xim1 = x(i-1)+x(i)
         x(i) = x(i)-x(i-1)
         x(i-1) = xim1
  101 continue
      x(1) = x(1)+x(1)
      modn = mod(n,2)
      if (modn .eq. 0) x(n) = x(n)+x(n)
c
      call drfftb (n,x,xh)
c
      do 102 k=2,ns2
         kc = np2-k
         xh(k) = w(k-1)*x(kc)+w(kc-1)*x(k)
         xh(kc) = w(k-1)*x(k)-w(kc-1)*x(kc)
  102 continue
c
      if (modn .eq. 0) x(ns2+1) = w(ns2)*(x(ns2+1)+x(ns2+1))
      do 103 k=2,ns2
         kc = np2-k
         x(k) = xh(k)+xh(kc)
         x(kc) = xh(k)-xh(kc)
  103 continue
      x(1) = x(1)+x(1)
c
      return
      end
