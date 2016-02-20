module module_fast_fourier_transform

    use, intrinsic :: iso_fortran_env, only: &
        wp     => REAL64, &
        ip     => INT32, &
        stderr => ERROR_UNIT

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: fft991
    public :: set99


contains

    subroutine fft99 (a,work,trigs,ifax,inc,jump,n,lot,isign)

        ! purpose      performs multiple fast fourier transforms.  this package
        !              will perform a number of simultaneous real/half-complex
        !              periodic fourier transforms or corresponding inverse
        !              transforms, i.e.  given a set of real data vectors, the
        !              package returns a set of 'half-complex' fourier
        !              coefficient vectors, or vice versa.  the length of the
        !              transforms must be an even number greater than 4 that has
        !              no other factors except possibly powers of 2, 3, and 5.
        !              this is an all-fortran version of a optimized routine
        !              fft99 written for xmp/ymps by dr. clive temperton of
        !              ecmwf.
        !
        !              the package fft99f contains several user-level routines:
        !
        !            subroutine set99
        !                an initialization routine that must be called once
        !                before a sequence of calls to the fft routines
        !                (provided that n is not changed).
        !
        !            subroutines fft99 and fft991
        !                two fft routines that return slightly different
        !                arrangements of the data in gridpoint space.
        !
        ! usage        let n be of the form 2**p * 3**q * 5**r, where p .ge. 1,
        !              q .ge. 0, and r .ge. 0.  then a typical sequence of
        !              calls to transform a given set of real vectors of length
        !              n to a set of 'half-complex' fourier coefficient vectors
        !              of length n is
        !
        !                   dimension ifax(13),trigs(3*n/2+1),a(m*(n+2)),
        !                  +          work(m*(n+1))
        !
        !                   call set99 (trigs, ifax, n)
        !                   call fft99 (a,work,trigs,ifax,inc,jump,n,m,isign)
        !
        !              see the individual write-ups for set99, fft99, and
        !              fft991 below, for a detailed description of the
        !              arguments.
        !
        ! history      the package was written by clive temperton at ecmwf in
        !              november, 1978.  it was modified, documented, and tested
        !              for ncar by russ rew in september, 1980.
        !
        !-----------------------------------------------------------------------
        !
        ! subroutine set99 (trigs, ifax, n)
        !
        ! purpose      a set-up routine for fft99 and fft991.  it need only be
        !              called once before a sequence of calls to the fft
        !              routines (provided that n is not changed).
        !
        ! argument     ifax(13),trigs(3*n/2+1)
        ! dimensions
        !
        ! arguments
        !
        ! on input     trigs
        !               a floating point array of dimension 3*n/2 if n/2 is
        !               even, or 3*n/2+1 if n/2 is odd.
        !
        !              ifax
        !               an integer array.  the number of elements actually used
        !               will depend on the factorization of n.  dimensioning
        !               ifax for 13 suffices for all n less than a million.
        !
        !              n
        !               an even number greater than 4 that has no prime factor
        !               greater than 5.  n is the length of the transforms (see
        !               the documentation for fft99 and fft991 for the
        !               definitions of the transforms).
        !
        ! on output    ifax
        !               contains the factorization of n/2.  ifax(1) is the
        !               number of factors, and the factors themselves are stored
        !               in ifax(2),ifax(3),...  if set99 is called with n odd,
        !               or if n has any prime factors greater than 5, ifax(1)
        !               is set to -99.
        !
        !              trigs
        !               an array of trigonometric function values subsequently
        !               used by the fft routines.
        !
        !-----------------------------------------------------------------------
        !
        ! subroutine fft991 (a,work,trigs,ifax,inc,jump,n,m,isign)
        !                       and
        ! subroutine fft99 (a,work,trigs,ifax,inc,jump,n,m,isign)
        !
        ! purpose      perform a number of simultaneous real/half-complex
        !              periodic fourier transforms or corresponding inverse
        !              transforms, using ordinary spatial order of gridpoint
        !              values (fft991) or explicit cyclic continuity in the
        !              gridpoint values (fft99).  given a set
        !              of real data vectors, the package returns a set of
        !              'half-complex' fourier coefficient vectors, or vice
        !              versa.  the length of the transforms must be an even
        !              number that has no other factors except possibly powers
        !              of 2, 3, and 5.  this is an all-fortran version of
        !              optimized routine fft991 written for xmp/ymps by
        !              dr. clive temperton of ecmwf.
        !
        ! argument     a(m*(n+2)), work(m*(n+1)), trigs(3*n/2+1), ifax(13)
        ! dimensions
        !
        ! arguments
        !
        ! on input     a
        !               an array of length m*(n+2) containing the input data
        !               or coefficient vectors.  this array is overwritten by
        !               the results.
        !
        !              work
        !               a work array of dimension m*(n+1)
        !
        !              trigs
        !               an array set up by set99, which must be called first.
        !
        !              ifax
        !               an array set up by set99, which must be called first.
        !
        !              inc
        !               the increment (in words) between successive elements of
        !               each data or coefficient vector (e.g.  inc=1 for
        !               consecutively stored data).
        !
        !              jump
        !               the increment (in words) between the first elements of
        !               successive data or coefficient vectors.  on crays,
        !               try to arrange data so that jump is not a multiple of 8
        !               (to avoid memory bank conflicts).  for clarification of
        !               inc and jump, see the examples below.
        !
        !              n
        !               the length of each transform (see definition of
        !               transforms, below).
        !
        !              m
        !               the number of transforms to be done simultaneously.
        !
        !              isign
        !               = +1 for a transform from fourier coefficients to
        !                    gridpoint values.
        !               = -1 for a transform from gridpoint values to fourier
        !                    coefficients.
        !
        ! on output    a
        !               if isign = +1, and m coefficient vectors are supplied
        !               each containing the sequence:
        !
        !               a(0),b(0),a(1),b(1),...,a(n/2),b(n/2)  (n+2 values)
        !
        !               then the result consists of m data vectors each
        !               containing the corresponding n+2 gridpoint values:
        !
        !               for fft991, x(0), x(1), x(2),...,x(n-1),0,0.
        !               for fft99, x(n-1),x(0),x(1),x(2),...,x(n-1),x(0).
        !                   (explicit cyclic continuity)
        !
        !               when isign = +1, the transform is defined by:
        !                 x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
        !                 where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
        !                 and i=sqrt (-1)
        !
        !               if isign = -1, and m data vectors are supplied each
        !               containing a sequence of gridpoint values x(j) as
        !               defined above, then the result consists of m vectors
        !               each containing the corresponding fourier cofficients
        !               a(k), b(k), 0 .le. k .le n/2.
        !
        !               when isign = -1, the inverse transform is defined by:
        !                 c(k)=(1/n)*sum(j=0,...,n-1)(x(j)*exp(-2*i*j*k*pi/n))
        !                 where c(k)=a(k)+i*b(k) and i=sqrt(-1)
        !
        !               a call with isign=+1 followed by a call with isign=-1
        !               (or vice versa) returns the original data.
        !
        !               note: the fact that the gridpoint values x(j) are real
        !               implies that b(0)=b(n/2)=0.  for a call with isign=+1,
        !               it is not actually necessary to supply these zeros.
        !
        ! examples      given 19 data vectors each of length 64 (+2 for explicit
        !               cyclic continuity), compute the corresponding vectors of
        !               fourier coefficients.  the data may, for example, be
        !               arranged like this:
        !
        ! first data   a(1)=    . . .                a(66)=             a(70)
        ! vector       x(63) x(0) x(1) x(2) ... x(63) x(0)  (4 empty locations)
        !
        ! second data  a(71)=   . . .                                  a(140)
        ! vector       x(63) x(0) x(1) x(2) ... x(63) x(0)  (4 empty locations)
        !
        !               and so on.  here inc=1, jump=70, n=64, m=19, isign=-1,
        !               and fft99 should be used (because of the explicit cyclic
        !               continuity).
        !
        !               alternatively the data may be arranged like this:
        !
        !                first         second                          last
        !                data          data                            data
        !                vector        vector                          vector
        !
        !                 a(1)=         a(2)=                           a(19)=
        !
        !                 x(63)         x(63)       . . .               x(63)
        !        a(20)=   x(0)          x(0)        . . .               x(0)
        !        a(39)=   x(1)          x(1)        . . .               x(1)
        !                  .             .                               .
        !                  .             .                               .
        !                  .             .                               .
        !
        !               in which case we have inc=19, jump=1, and the remaining
        !               parameters are the same as before.  in either case, each
        !               coefficient vector overwrites the corresponding input
        !               data vector.
        !
        !-----------------------------------------------------------------------
        integer (ip), intent(in)    :: inc
        integer (ip), intent(in)    ::jump
        integer (ip), intent(in)    ::n
        integer (ip), intent(in)    ::lot
        integer (ip), intent(in)    ::isign
        integer (ip), intent(in) :: ifax(:)
        real,    intent(in)    :: trigs(:)
        real,    intent(in out) :: a(*)
        real,    intent(in out) ::work(*)

        !     dimension a(n),work(n),trigs(n),ifax(1)
        !
        !     subroutine "fft99" - multiple fast real periodic transform
        !     corresponding to old scalar routine fft9
        !     procedure used to convert to half-length complex transform
        !     is given by cooley, lewis and welch (j. sound vib., vol. 12
        !     (1970), 315-337)
        !
        !     a is the array containing input and output data
        !     work is an area of size (n+1)*lot
        !     trigs is a previously prepared list of trig function values
        !     ifax is a previously prepared list of factors of n/2
        !     inc is the increment within each data 'vector'
        !         (e.g. inc=1 for consecutively stored data)
        !     jump is the increment between the start of each data vector
        !     n is the length of the data vectors
        !     lot is the number of data vectors
        !     isign = +1 for transform from spectral to gridpoint
        !           = -1 for transform from gridpoint to spectral
        !
        !     ordering of coefficients:
        !         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
        !         where b(0)=b(n/2)=0; (n+2) locations required
        !
        !     ordering of data:
        !         x(n-1),x(0),x(1),x(2),...,x(n),x(0)
        !         i.e. explicit cyclic continuity; (n+2) locations required
        !
        !     vectorization is achieved on cray by doing the transforms in
        !     parallel
        !
        !     *** n.b. n is assumed to be an even number
        !
        !     definition of transforms:
        !     -------------------------
        !
        !     isign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
        !         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
        !
        !     isign=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
        !               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
        !

        integer:: nfax
        integer:: nx
        integer:: nh
        integer:: ink
        integer:: igo
        integer:: ibase
        integer:: jbase
        integer:: i
        integer:: j
        integer:: k
        integer:: l
        integer:: m
        integer:: ia
        integer:: la
        integer:: ib

        nfax=ifax(1)
        nx=n+1
        nh=n/2
        ink=inc+inc
        if (isign == +1) go to 30

        !   if necessary, transfer data to work area
        igo=50
        if (mod(nfax,2) == 1) goto 40
        ibase=inc+1
        jbase=1
        do l=1,lot
            i=ibase
            j=jbase
            do m=1,n
                work(j)=a(i)
                i=i+inc
                j=j+1
            end do
            ibase=ibase+jump
            jbase=jbase+nx
        end do

        igo=60
        go to 40

    !   preprocessing (isign=+1)
    !   ------------------------

30  continue
    call fft99a(a,work,trigs,inc,jump,n,lot)
    igo=60

!   complex transform
!   -----------------

40 continue
   ia=inc+1
   la=1
   do 80 k=1,nfax
       if (igo == 60) go to 60
50 continue
   call vpassm (a(ia),a(ia+inc),work(1),work(2),trigs, &
       ink,2,jump,nx,lot,nh,ifax(k+1),la)
   igo=60
   go to 70
60 continue
   call vpassm (work(1),work(2),a(ia),a(ia+inc),trigs, &
       2,ink,nx,jump,lot,nh,ifax(k+1),la)
   igo=50
70 continue
   la=la*ifax(k+1)
80 continue

   if (isign == -1) go to 130

   ! if necessary, transfer data from work area

   if (mod(nfax,2)/=1) then
       ibase=1
       jbase=ia
       do l=1,lot
           i=ibase
           j=jbase
           do m=1,n
               a(j)=work(i)
               i=i+1
               j=j+inc
           end do
           ibase=ibase+nx
           jbase=jbase+jump
       end do
   end if

   !   fill in cyclic boundary points
   ia=1
   ib=n*inc+1
   do l=1,lot
       a(ia)=a(ib)
       a(ib+inc)=a(ia+inc)
       ia=ia+jump
       ib=ib+jump
   end do
   go to 140

   !   postprocessing (isign=-1):
   !   --------------------------

130 continue
    call fft99b(work,a,trigs,inc,jump,n,lot)

140 continue

end subroutine fft99

!##########################################################################

subroutine fft99a (a,work,trigs,inc,jump,n,lot)
    integer (ip), intent(in)    :: inc
    integer (ip), intent(in)    ::jump
    integer (ip), intent(in)    ::n
    integer (ip), intent(in)    ::lot
    real,    intent(in)    :: trigs(:)
    real,    intent(in out) :: a(*)
    real,    intent(in out) ::work(*)

    !     dimension a(n),work(n),trigs(n)
    !
    !     subroutine fft99a - preprocessing step for fft99, isign=+1
    !     (spectral to gridpoint transform)

    integer:: nh
    integer:: nx
    integer:: ink
    integer:: k
    integer:: l
    integer:: ia
    integer:: ib
    integer:: ja
    integer:: jb
    integer:: iabase
    integer:: ibbase
    integer:: jabase
    integer:: jbbase
    real:: c
    real:: s

    nh=n/2
    nx=n+1
    ink=inc+inc

    !   a(0) and a(n/2)
    ia=1
    ib=n*inc+1
    ja=1
    jb=2
    do l=1,lot
        work(ja)=a(ia)+a(ib)
        work(jb)=a(ia)-a(ib)
        ia=ia+jump
        ib=ib+jump
        ja=ja+nx
        jb=jb+nx
    end do

    !   remaining wavenumbers
    iabase=2*inc+1
    ibbase=(n-2)*inc+1
    jabase=3
    jbbase=n-1

    do k=3,nh,2
        ia=iabase
        ib=ibbase
        ja=jabase
        jb=jbbase
        c=trigs(n+k)
        s=trigs(n+k+1)
        do l=1,lot
            work(ja)=(a(ia)+a(ib))- &
                (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
            work(jb)=(a(ia)+a(ib))+ &
                (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
            work(ja+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))+ &
                (a(ia+inc)-a(ib+inc))
            work(jb+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))- &
                (a(ia+inc)-a(ib+inc))
            ia=ia+jump
            ib=ib+jump
            ja=ja+nx
            jb=jb+nx
        end do
        iabase=iabase+ink
        ibbase=ibbase-ink
        jabase=jabase+2
        jbbase=jbbase-2
    end do

    !   wavenumber n/4 (if it exists)
    if (iabase == ibbase) then
        ia=iabase
        ja=jabase
        do l=1,lot
            work(ja)=2.0*a(ia)
            work(ja+1)=-2.0*a(ia+inc)
            ia=ia+jump
            ja=ja+nx
        end do
    end if

end subroutine fft99a

!##########################################################################

subroutine fft99b (work,a,trigs,inc,jump,n,lot)
    integer (ip), intent(in)    :: inc
    integer (ip), intent(in)    ::jump
    integer (ip), intent(in)    ::n
    integer (ip), intent(in)    ::lot
    real,    intent(in)    :: trigs(:)
    real,    intent(in out) :: a(*)
    real,    intent(in out) ::work(*)

    !     dimension work(n),a(n),trigs(n)
    !
    !     subroutine fft99b - postprocessing step for fft99, isign=-1
    !     (gridpoint to spectral transform)

    integer:: nh
    integer:: nx
    integer:: ink
    integer:: k
    integer:: l
    integer:: ia
    integer:: ib
    integer:: ja
    integer:: jb
    integer:: iabase
    integer:: ibbase
    integer:: jabase
    integer:: jbbase
    real:: scale
    real:: c
    real:: s

    nh=n/2
    nx=n+1
    ink=inc+inc

    !   a(0) and a(n/2)
    scale=1.0/real(n)
    ia=1
    ib=2
    ja=1
    jb=n*inc+1
    do l=1,lot
        a(ja)=scale*(work(ia)+work(ib))
        a(jb)=scale*(work(ia)-work(ib))
        a(ja+inc)=0.0
        a(jb+inc)=0.0
        ia=ia+nx
        ib=ib+nx
        ja=ja+jump
        jb=jb+jump
    end do

    !   remaining wavenumbers
    scale=0.5*scale
    iabase=3
    ibbase=n-1
    jabase=2*inc+1
    jbbase=(n-2)*inc+1

    do k=3,nh,2
        ia=iabase
        ib=ibbase
        ja=jabase
        jb=jbbase
        c=trigs(n+k)
        s=trigs(n+k+1)
        do l=1,lot
            a(ja)=scale*((work(ia)+work(ib)) &
                +(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
            a(jb)=scale*((work(ia)+work(ib)) &
                -(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
            a(ja+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1))) &
                +(work(ib+1)-work(ia+1)))
            a(jb+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1))) &
                -(work(ib+1)-work(ia+1)))
            ia=ia+nx
            ib=ib+nx
            ja=ja+jump
            jb=jb+jump
        end do
        iabase=iabase+2
        ibbase=ibbase-2
        jabase=jabase+ink
        jbbase=jbbase-ink
    end do

    !   wavenumber n/4 (if it exists)
    if (iabase == ibbase) then
        ia=iabase
        ja=jabase
        scale=2.0*scale
        do l=1,lot
            a(ja)=scale*work(ia)
            a(ja+inc)=-scale*work(ia+1)
            ia=ia+nx
            ja=ja+jump
        end do
    end if

end subroutine fft99b

!##########################################################################

subroutine fft991(a,work,trigs,ifax,inc,jump,n,lot,isign)
    integer (ip), intent(in)    :: inc
    integer (ip), intent(in)    ::jump
    integer (ip), intent(in)    ::n
    integer (ip), intent(in)    ::lot
    integer (ip), intent(in)    ::isign
    integer (ip), intent(in) :: ifax(:)
    real,    intent(in)    :: trigs(:)
    real,    intent(in out) :: a(*)
    real,    intent(in out) ::work(*)

    !     dimension a(n),work(n),trigs(n),ifax(1)
    !
    !     subroutine "fft991" - multiple real/half-complex periodic
    !     fast fourier transform
    !
    !     same as fft99 except that ordering of data corresponds to
    !     that in mperform_multiple_real_fft2
    !
    !     procedure used to convert to half-length complex transform
    !     is given by cooley, lewis and welch (j. sound vib., vol. 12
    !     (1970), 315-337)
    !
    !     a is the array containing input and output data
    !     work is an area of size (n+1)*lot
    !     trigs is a previously prepared list of trig function values
    !     ifax is a previously prepared list of factors of n/2
    !     inc is the increment within each data 'vector'
    !         (e.g. inc=1 for consecutively stored data)
    !     jump is the increment between the start of each data vector
    !     n is the length of the data vectors
    !     lot is the number of data vectors
    !     isign = +1 for transform from spectral to gridpoint
    !           = -1 for transform from gridpoint to spectral
    !
    !     ordering of coefficients:
    !         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
    !         where b(0)=b(n/2)=0; (n+2) locations required
    !
    !     ordering of data:
    !         x(0),x(1),x(2),...,x(n-1)
    !
    !     vectorization is achieved on cray by doing the transforms in
    !     parallel
    !
    !     *** n.b. n is assumed to be an even number
    !
    !     definition of transforms:
    !     -------------------------
    !
    !     isign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
    !         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
    !
    !     isign=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
    !               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
    !

    integer:: nfax
    integer:: nx
    integer:: nh
    integer:: ink
    integer:: igo
    integer:: ibase
    integer:: jbase
    integer:: i
    integer:: j
    integer:: k
    integer:: l
    integer:: m
    integer:: ia
    integer:: la
    integer:: ib


    nfax=ifax(1)
    nx=n+1
    nh=n/2
    ink=inc+inc
    if (isign == +1) go to 30

    !     if necessary, transfer data to work area
    igo=50
    if (mod(nfax,2) == 1) goto 40
    ibase=1
    jbase=1
    do l=1,lot
        i=ibase
        j=jbase
        do m=1,n
            work(j)=a(i)
            i=i+inc
            j=j+1
        end do
        ibase=ibase+jump
        jbase=jbase+nx
    end do
    !
    igo=60
    go to 40
!
!   preprocessing (isign=+1)
!   ------------------------
!
30 continue
   call fft99a(a,work,trigs,inc,jump,n,lot)
   igo=60
   !
   !   complex transform
   !   -----------------
   !
40 continue
   ia=1
   la=1
   do k=1,nfax
       if (igo == 60) go to 60
50 continue
   call vpassm (a(ia),a(ia+inc),work(1),work(2),trigs, &
       ink,2,jump,nx,lot,nh,ifax(k+1),la)
   igo=60
   go to 70
60 continue
   call vpassm (work(1),work(2),a(ia),a(ia+inc),trigs, &
       2,ink,nx,jump,lot,nh,ifax(k+1),la)
   igo=50
70 continue
   la=la*ifax(k+1)
   end do

   if (isign == -1) go to 130

   ! if necessary, transfer data from work area
   if (mod(nfax,2)/=1) then
       ibase=1
       jbase=1
       do l=1,lot
           i=ibase
           j=jbase
           do m=1,n
               a(j)=work(i)
               i=i+1
               j=j+inc
           end do
           ibase=ibase+nx
           jbase=jbase+jump
       end do
   end if

   !   fill in zeros at end
   ib=n*inc+1
   do l=1,lot
       a(ib)=0.0
       a(ib+inc)=0.0
       ib=ib+jump
   end do
   go to 140

   !     postprocessing (isign=-1):
   !     --------------------------

130 continue
    call fft99b (work,a,trigs,inc,jump,n,lot)

140 continue

end subroutine fft991

!##########################################################################

subroutine set99 (trigs, ifax, n)
    integer (ip), intent(in)  :: n
    integer (ip), intent(out) :: ifax(:)
    real,    intent(out) :: trigs(:)

    !     dimension ifax(13),trigs(1)
    !
    ! mode 3 is used for real/half-complex transforms.  it is possible
    ! to do complex/complex transforms with other values of mode, but
    ! documentation of the details were not available when this routine
    ! was written.
    !
    integer:: mode = 3
    integer:: i

    call fax (ifax, n, mode)
    i = ifax(1)
    if (ifax(i+1) > 5 .or. n <= 4) ifax(1) = -99
    if (ifax(1) <= 0 ) then
        write(6,*) ' set99 -- invalid n'
        stop
    end if
    call fftrig (trigs, n, mode)

end subroutine set99

!##########################################################################

subroutine fax (ifax,n,mode)
    integer (ip), intent(out) :: ifax(:)
    integer (ip), intent(in)  :: n
    integer (ip), intent(in)  :: mode

    integer:: nn
    integer:: k
    integer:: l
    integer:: inc
    integer:: nfax
    integer:: ii
    integer:: istop
    integer:: i
    integer:: item

    nn=n
    if (iabs(mode) == 1) go to 10
    if (iabs(mode) == 8) go to 10
    nn=n/2
    if ((nn+nn) == n) go to 10
    ifax(1)=-99
    return
10  k=1
    !     test for factors of 4
20  if (mod(nn,4)/=0) go to 30
    k=k+1
    ifax(k)=4
    nn=nn/4
    if (nn == 1) go to 80
    go to 20
    !     test for extra factor of 2
30  if (mod(nn,2)/=0) go to 40
    k=k+1
    ifax(k)=2
    nn=nn/2
    if (nn == 1) go to 80
    !     test for factors of 3
40  if (mod(nn,3)/=0) go to 50
    k=k+1
    ifax(k)=3
    nn=nn/3
    if (nn == 1) go to 80
    go to 40
    !     now find remaining factors
50  l=5
    inc=2
    !     inc alternately takes on values 2 and 4
60  if (mod(nn,l)/=0) go to 70
    k=k+1
    ifax(k)=l
    nn=nn/l
    if (nn == 1) go to 80
    go to 60
70  l=l+inc
    inc=6-inc
    go to 60
80  ifax(1)=k-1
    !     ifax(1) contains number of factors
    nfax=ifax(1)
    !     sort factors into ascending order
    if (nfax == 1) go to 110
    do 100 ii=2,nfax
        istop=nfax+2-ii
        do 90 i=2,istop
            if (ifax(i+1)>=ifax(i)) go to 90
            item=ifax(i)
            ifax(i)=ifax(i+1)
            ifax(i+1)=item
90      continue
100 continue
110 continue

end subroutine fax

!##########################################################################

subroutine fftrig (trigs,n,mode)
    real,    intent(out) :: trigs(:)
    integer (ip), intent(in)  :: n
    integer (ip), intent(in)  :: mode

    real:: del
    real:: angle
    real::pi
    integer:: imode
    integer:: nn
    integer:: nh
    integer:: i
    integer:: l
    integer:: la

    pi = 4.*atan(1.0)
    imode=iabs(mode)
    nn=n
    if (imode>1.and.imode<6) nn=n/2
    del=(pi+pi)/real(nn)
    l=nn+nn
    do i=1,l,2
        angle=0.5*real(i-1)*del
        trigs(i)=cos(angle)
        trigs(i+1)=sin(angle)
    end do
    if (imode == 1) return
    if (imode == 8) return

    del=0.5*del
    nh=(nn+1)/2
    l=nh+nh
    la=nn+nn
    do i=1,l,2
        angle=0.5*real(i-1)*del
        trigs(la+i)=cos(angle)
        trigs(la+i+1)=sin(angle)
    end do
    if (imode<=3) return

    del=0.5*del
    la=la+nn
    if (mode/=5) then
        do i=2,nn
            angle=real(i-1)*del
            trigs(la+i)=2.0*sin(angle)
        end do
        return
    end if

    del=0.5*del
    do i=2,n
        angle=real(i-1)*del
        trigs(la+i)=sin(angle)
    end do

end subroutine fftrig

!##########################################################################

subroutine vpassm (a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,la)
    integer (ip), intent(in)  :: inc1
    integer (ip), intent(in)  :: inc2
    integer (ip), intent(in)  :: inc3
    integer (ip), intent(in)  :: inc4
    integer (ip), intent(in)  :: lot
    integer (ip), intent(in)  :: n
    integer (ip), intent(in)  :: ifac
    integer (ip), intent(in)  :: la
    real,    intent(in)  :: a(n)
    real,    intent(in)  ::b(n)
    real,    intent(in)  ::trigs(n)
    real,    intent(out) :: c(n)
    real,    intent(out) ::d(n)
    !
    !     subroutine "vpassm" - multiple version of "vpassa"
    !     performs one pass through data
    !     as part of multiple complex fft routine
    !     a is first real input vector
    !     b is first imaginary input vector
    !     c is first real output vector
    !     d is first imaginary output vector
    !     trigs is precalculated table of sines " cosines
    !     inc1 is addressing increment for a and b
    !     inc2 is addressing increment for c and d
    !     inc3 is addressing increment between a"s & b"s
    !     inc4 is addressing increment between c"s & d"s
    !     lot is the number of vectors
    !     n is length of vectors
    !     ifac is current factor of n
    !     la is product of previous factors
    !

    real:: sin36=0.587785252292473
    real:: cos36=0.809016994374947
    real:: sin72=0.951056516295154
    real:: cos72=0.309016994374947424102293417182819058860154589902881431067
    real:: sin60=0.866025403784438646763723170752936183471402626905190314027

    integer:: i
    integer:: j
    integer:: k
    integer:: l
    integer:: m
    integer:: iink
    integer:: jink
    integer:: jump
    integer:: ibase
    integer:: jbase
    integer:: igo
    integer:: ijk
    integer:: la1
    integer:: ia
    integer:: ja
    integer:: ib
    integer:: jb
    integer:: kb
    integer:: ic
    integer:: jc
    integer:: kc
    integer:: id
    integer:: jd
    integer:: kd
    integer:: ie
    integer:: je
    integer:: ke
    real:: c1
    real:: s1
    real:: c2
    real:: s2
    real:: c3
    real:: s3
    real:: c4
    real:: s4

    m=n/ifac
    iink=m*inc1
    jink=la*inc2
    jump=(ifac-1)*jink
    ibase=0
    jbase=0
    igo=ifac-1
    if (igo>4) return
    !del  go to (10,50,90,130),igo

    select case (igo)

        !   coding for factor 2

        case (1)
10          ia=1
            ja=1
            ib=ia+iink
            jb=ja+jink
            do 20 l=1,la
                i=ibase
                j=jbase
                do 15 ijk=1,lot
                    c(ja+j)=a(ia+i)+a(ib+i)
                    d(ja+j)=b(ia+i)+b(ib+i)
                    c(jb+j)=a(ia+i)-a(ib+i)
                    d(jb+j)=b(ia+i)-b(ib+i)
                    i=i+inc3
                    j=j+inc4
15              continue
                ibase=ibase+inc1
                jbase=jbase+inc2
20          continue
            if (la == m) return
            la1=la+1
            jbase=jbase+jump
            do 40 k=la1,m,la
                kb=k+k-2
                c1=trigs(kb+1)
                s1=trigs(kb+2)
                do 30 l=1,la
                    i=ibase
                    j=jbase
                    do 25 ijk=1,lot
                        c(ja+j)=a(ia+i)+a(ib+i)
                        d(ja+j)=b(ia+i)+b(ib+i)
                        c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
                        d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
                        i=i+inc3
                        j=j+inc4
25                  continue
                    ibase=ibase+inc1
                    jbase=jbase+inc2
30              continue
                jbase=jbase+jump
40          continue
        !     return

        !   coding for factor 3

        case (2)
50          ia=1
            ja=1
            ib=ia+iink
            jb=ja+jink
            ic=ib+iink
            jc=jb+jink
            do 60 l=1,la
                i=ibase
                j=jbase
                do 55 ijk=1,lot
                    c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
                    d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
                    c(jb+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))
                    c(jc+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))
                    d(jb+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i)))
                    d(jc+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))
                    i=i+inc3
                    j=j+inc4
55              continue
                ibase=ibase+inc1
                jbase=jbase+inc2
60          continue
            if (la == m) return
            la1=la+1
            jbase=jbase+jump
            do 80 k=la1,m,la
                kb=k+k-2
                kc=kb+kb
                c1=trigs(kb+1)
                s1=trigs(kb+2)
                c2=trigs(kc+1)
                s2=trigs(kc+2)
                do 70 l=1,la
                    i=ibase
                    j=jbase
                    do 65 ijk=1,lot
                        c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
                        d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
                        c(jb+j)=                                                           &
                            c1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))) &
                            -s1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
                        d(jb+j)=                                                           &
                            s1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))) &
                            +c1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
                        c(jc+j)=                                                           &
                            c2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))) &
                            -s2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
                        d(jc+j)=                                                           &
                            s2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))) &
                            +c2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
                        i=i+inc3
                        j=j+inc4
65                  continue
                    ibase=ibase+inc1
                    jbase=jbase+inc2
70              continue
                jbase=jbase+jump
80          continue
        !     return

        !   coding for factor 4

        case (3)
90          ia=1
            ja=1
            ib=ia+iink
            jb=ja+jink
            ic=ib+iink
            jc=jb+jink
            id=ic+iink
            jd=jc+jink
            do 100 l=1,la
                i=ibase
                j=jbase
                do 95 ijk=1,lot
                    c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
                    c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
                    d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
                    d(jc+j)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
                    c(jb+j)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
                    c(jd+j)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
                    d(jb+j)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
                    d(jd+j)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
                    i=i+inc3
                    j=j+inc4
95              continue
                ibase=ibase+inc1
                jbase=jbase+inc2
100         continue
            if (la == m) return
            la1=la+1
            jbase=jbase+jump
            do 120 k=la1,m,la
                kb=k+k-2
                kc=kb+kb
                kd=kc+kb
                c1=trigs(kb+1)
                s1=trigs(kb+2)
                c2=trigs(kc+1)
                s2=trigs(kc+2)
                c3=trigs(kd+1)
                s3=trigs(kd+2)
                do 110 l=1,la
                    i=ibase
                    j=jbase
                    do 105 ijk=1,lot
                        c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
                        d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
                        c(jc+j)=                                     &
                            c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
                            -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
                        d(jc+j)=                                     &
                            s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
                            +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
                        c(jb+j)=                                     &
                            c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
                            -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
                        d(jb+j)=                                     &
                            s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
                            +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
                        c(jd+j)=                                     &
                            c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
                            -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
                        d(jd+j)=                                     &
                            s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
                            +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
                        i=i+inc3
                        j=j+inc4
105                 continue
                    ibase=ibase+inc1
                    jbase=jbase+inc2
110             continue
                jbase=jbase+jump
120         continue
        !     return

        !   coding for factor 5

        case (4)
130         ia=1
            ja=1
            ib=ia+iink
            jb=ja+jink
            ic=ib+iink
            jc=jb+jink
            id=ic+iink
            jd=jc+jink
            ie=id+iink
            je=jd+jink
            do 140 l=1,la
                i=ibase
                j=jbase
                do 135 ijk=1,lot
                    c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
                    d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
                    c(jb+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
                        -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
                    c(je+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
                        +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
                    d(jb+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
                        +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
                    d(je+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
                        -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
                    c(jc+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
                        -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
                    c(jd+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
                        +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
                    d(jc+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
                        +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
                    d(jd+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
                        -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
                    i=i+inc3
                    j=j+inc4
135             continue
                ibase=ibase+inc1
                jbase=jbase+inc2
140         continue
            if (la == m) return
            la1=la+1
            jbase=jbase+jump
            do 160 k=la1,m,la
                kb=k+k-2
                kc=kb+kb
                kd=kc+kb
                ke=kd+kb
                c1=trigs(kb+1)
                s1=trigs(kb+2)
                c2=trigs(kc+1)
                s2=trigs(kc+2)
                c3=trigs(kd+1)
                s3=trigs(kd+2)
                c4=trigs(ke+1)
                s4=trigs(ke+2)
                do 150 l=1,la
                    i=ibase
                    j=jbase
                    do 145 ijk=1,lot
                        c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
                        d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
                        c(jb+j)=                                                          &
                            c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
                            -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
                            -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
                            +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
                        d(jb+j)=                                                          &
                            s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
                            -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
                            +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
                            +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
                        c(je+j)=                                                          &
                            c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
                            +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
                            -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
                            -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
                        d(je+j)=                                                          &
                            s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
                            +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
                            +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
                            -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
                        c(jc+j)=                                                          &
                            c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
                            -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
                            -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
                            +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
                        d(jc+j)=                                                          &
                            s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
                            -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
                            +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
                            +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
                        c(jd+j)=                                                          &
                            c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
                            +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
                            -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
                            -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
                        d(jd+j)=                                                          &
                            s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
                            +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
                            +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
                            -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
                        i=i+inc3
                        j=j+inc4
145                 continue
                    ibase=ibase+inc1
                    jbase=jbase+inc2
150             continue
                jbase=jbase+jump
160         continue

    end select

end subroutine vpassm

end module module_fast_fourier_transform
