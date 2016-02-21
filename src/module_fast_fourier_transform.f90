module module_fast_fourier_transform

    use, intrinsic :: iso_fortran_env, only: &
        wp     => REAL64, &
        ip     => INT32, &
        stderr => ERROR_UNIT

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: perform_fft991
    public :: initialize_fft99


contains
     !
     !*****************************************************************************************
     !
    subroutine perform_preprocessing_step_for_fft99 (a,work,trigs,inc,jump,n,lot)
        integer (ip), intent(in)    :: inc
        integer (ip), intent(in)    ::jump
        integer (ip), intent(in)    ::n
        integer (ip), intent(in)    ::lot
        real (wp),    intent(in)    :: trigs(:)
        real (wp),    intent(in out) :: a(*)
        real (wp),    intent(in out) ::work(*)

        !     dimension a(n),work(n),trigs(n)
        !
        !     subroutine fft99a - preprocessing step for fft99, isign=+1
        !     (spectral to gridpoint transform)

        integer (ip) :: nh
        integer (ip) :: nx
        integer (ip) :: ink
        integer (ip) :: k
        integer (ip) :: l
        integer (ip) :: ia
        integer (ip) :: ib
        integer (ip) :: ja
        integer (ip) :: jb
        integer (ip) :: iabase
        integer (ip) :: ibbase
        integer (ip) :: jabase
        integer (ip) :: jbbase
        real (wp) :: c
        real (wp) :: s

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
                work(ja)=2.0_wp*a(ia)
                work(ja+1)=-2.0_wp*a(ia+inc)
                ia=ia+jump
                ja=ja+nx
            end do
        end if

    end subroutine perform_preprocessing_step_for_fft99
        !
        !*****************************************************************************************
        !
    subroutine perform_postprocessing_step_for_fft99(work,a,trigs,inc,jump,n,lot)
        integer (ip), intent(in)    :: inc
        integer (ip), intent(in)    ::jump
        integer (ip), intent(in)    ::n
        integer (ip), intent(in)    ::lot
        real (wp),    intent(in)    :: trigs(:)
        real (wp),    intent(in out) :: a(*)
        real (wp),    intent(in out) ::work(*)

        !     dimension work(n),a(n),trigs(n)
        !
        !     subroutine fft99b - postprocessing step for fft99, isign=-1
        !     (gridpoint to spectral transform)

        integer (ip) :: nh
        integer (ip) :: nx
        integer (ip) :: ink
        integer (ip) :: k
        integer (ip) :: l
        integer (ip) :: ia
        integer (ip) :: ib
        integer (ip) :: ja
        integer (ip) :: jb
        integer (ip) :: iabase
        integer (ip) :: ibbase
        integer (ip) :: jabase
        integer (ip) :: jbbase
        real (wp) :: scale_constant
        real (wp) :: c
        real (wp) :: s

        nh=n/2
        nx=n+1
        ink=inc+inc

        !   a(0) and a(n/2)
        scale_constant = 1.0_wp/real(n, kind=wp)
        ia=1
        ib=2
        ja=1
        jb=n*inc+1
        do l=1,lot
            a(ja)=scale_constant*(work(ia)+work(ib))
            a(jb)=scale_constant*(work(ia)-work(ib))
            a(ja+inc)=0.0_wp
            a(jb+inc)=0.0_wp
            ia=ia+nx
            ib=ib+nx
            ja=ja+jump
            jb=jb+jump
        end do

        !   remaining wavenumbers
        scale_constant=0.5_wp*scale_constant
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
                a(ja)=scale_constant*((work(ia)+work(ib)) &
                    +(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
                a(jb)=scale_constant*((work(ia)+work(ib)) &
                    -(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
                a(ja+inc)=scale_constant*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1))) &
                    +(work(ib+1)-work(ia+1)))
                a(jb+inc)=scale_constant*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1))) &
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
            scale_constant=2.0_wp*scale_constant
            do l=1,lot
                a(ja)=scale_constant*work(ia)
                a(ja+inc)=-scale_constant*work(ia+1)
                ia=ia+nx
                ja=ja+jump
            end do
        end if

    end subroutine perform_postprocessing_step_for_fft99
    !
    !*****************************************************************************************
    !
    subroutine perform_fft991(a,work,trigs,ifax,inc,jump,n,lot,isign)
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (wp),    intent(in out) :: a(*)
        real (wp),    intent(in out) :: work(*)
        real (wp),    intent(in)     :: trigs(:)
        integer (ip), intent(in)     :: ifax(:)
        integer (ip), intent(in)     :: inc
        integer (ip), intent(in)     :: jump
        integer (ip), intent(in)     :: n
        integer (ip), intent(in)     :: lot
        integer (ip), intent(in)     :: isign

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

        integer (ip) :: nfax
        integer (ip) :: nx
        integer (ip) :: nh
        integer (ip) :: ink
        integer (ip) :: igo
        integer (ip) :: ibase
        integer (ip) :: jbase
        integer (ip) :: i
        integer (ip) :: j
        integer (ip) :: k
        integer (ip) :: l
        integer (ip) :: m
        integer (ip) :: ia
        integer (ip) :: la
        integer (ip) :: ib


        nfax=ifax(1)
        nx=n+1
        nh=n/2
        ink=inc+inc
        if (isign == +1) goto 30

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
        goto 40
    !
    !   preprocessing (isign=+1)
    !   ------------------------
    !
30  continue
    call perform_preprocessing_step_for_fft99(a,work,trigs,inc,jump,n,lot)
    igo=60
!
!   complex transform
!   -----------------
!
40 continue
   ia=1
   la=1
   do k=1,nfax
       if (igo == 60) goto 60
50 continue
   call vpassm (a(ia),a(ia+inc),work(1),work(2),trigs, &
       ink,2,jump,nx,lot,nh,ifax(k+1),la)
   igo=60
   goto 70
60 continue
   call vpassm (work(1),work(2),a(ia),a(ia+inc),trigs, &
       2,ink,nx,jump,lot,nh,ifax(k+1),la)
   igo=50
70 continue
   la=la*ifax(k+1)
   end do

   if (isign == -1) goto 130

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
       a(ib)=0.0_wp
       a(ib+inc)=0.0_wp
       ib=ib+jump
   end do
   goto 14011

   !     postprocessing (isign=-1):
   !     --------------------------

130 continue
    call perform_postprocessing_step_for_fft99 (work,a,trigs,inc,jump,n,lot)

14011 continue

  end subroutine perform_fft991
  !
  !*****************************************************************************************
  !
  subroutine initialize_fft99 (trigs, ifax, n)
      integer (ip), intent(in)  :: n
      integer (ip), intent(out) :: ifax(:)
      real (wp),    intent(out) :: trigs(:)

      !     dimension ifax(13),trigs(1)
      !
      ! mode 3 is used for real/half-complex transforms.  it is possible
      ! to do complex/complex transforms with other values of mode, but
      ! documentation of the details were not available when this routine
      ! was written.
      !
      integer (ip) :: mode = 3
      integer (ip) :: i

      call fax (ifax, n, mode)
      i = ifax(1)
      if (ifax(i+1) > 5 .or. n <= 4) ifax(1) = -99
      if (ifax(1) <= 0 ) then
          error stop ' initialize_fft99 -- invalid n'
      end if
      call fftrig (trigs, n, mode)

  end subroutine initialize_fft99
  !
  !*****************************************************************************************
  !
  subroutine fax (ifax,n,mode)
      integer (ip), intent(out) :: ifax(:)
      integer (ip), intent(in)  :: n
      integer (ip), intent(in)  :: mode

      integer (ip) :: nn
      integer (ip) :: k
      integer (ip) :: l
      integer (ip) :: inc
      integer (ip) :: nfax
      integer (ip) :: ii
      integer (ip) :: istop
      integer (ip) :: i
      integer (ip) :: item

      nn=n
      if (iabs(mode) == 1) goto 10
      if (iabs(mode) == 8) goto 10
      nn=n/2
      if ((nn+nn) == n) goto 10
      ifax(1)=-99
      return
10    k=1
      !     test for factors of 4
20    if (mod(nn,4)/=0) goto 30333
      k=k+1
      ifax(k)=4
      nn=nn/4
      if (nn == 1) goto 80
      goto 20
      !     test for extra factor of 2
30333 if (mod(nn,2)/=0) goto 40
      k=k+1
      ifax(k)=2
      nn=nn/2
      if (nn == 1) goto 80
      !     test for factors of 3
40    if (mod(nn,3)/=0) goto 50
      k=k+1
      ifax(k)=3
      nn=nn/3
      if (nn == 1) goto 80
      goto 40
      !     now find remaining factors
50    l=5
      inc=2
      !     inc alternately takes on values 2 and 4
60    if (mod(nn,l)/=0) goto 70
      k=k+1
      ifax(k)=l
      nn=nn/l
      if (nn == 1) goto 80
      goto 60
70    l=l+inc
      inc=6-inc
      goto 60
80    ifax(1)=k-1
      !     ifax(1) contains number of factors
      nfax=ifax(1)
      !     sort factors into ascending order
      if (nfax == 1) goto 11011
      do ii=2,nfax
          istop=nfax+2-ii
          do i=2,istop
              if (ifax(i+1)>=ifax(i)) goto 90
              item=ifax(i)
              ifax(i)=ifax(i+1)
              ifax(i+1)=item
90        continue
          end do
1001  continue
      end do
11011 continue

  end subroutine fax
      !
      !*****************************************************************************************
      !
  subroutine fftrig (trigs,n,mode)
      real (wp),    intent(out) :: trigs(:)
      integer (ip), intent(in)  :: n
      integer (ip), intent(in)  :: mode

      real (wp) :: del
      real (wp) :: angle
      real (wp), parameter :: PI = acos( -1.0_wp )
      integer (ip) :: imode
      integer (ip) :: nn
      integer (ip) :: nh
      integer (ip) :: i
      integer (ip) :: l
      integer (ip) :: la

      imode=iabs(mode)
      nn=n
      if (imode>1.and.imode<6) nn=n/2
      del=(2.0_wp * PI)/real(nn, kind=wp)
      l=nn+nn
      do i=1,l,2
          angle=0.5_wp*real(i-1, kind=wp)*del
          trigs(i)=cos(angle)
          trigs(i+1)=sin(angle)
      end do
      if (imode == 1) return
      if (imode == 8) return

      del=0.5_wp*del
      nh=(nn+1)/2
      l=nh+nh
      la=nn+nn
      do i=1,l,2
          angle=0.5_wp*real(i-1,kind=wp)*del
          trigs(la+i)=cos(angle)
          trigs(la+i+1)=sin(angle)
      end do
      if (imode<=3) return

      del=0.5_wp*del
      la=la+nn
      if (mode/=5) then
          do i=2,nn
              angle=real(i-1,kind=wp)*del
              trigs(la+i)=2.0_wp*sin(angle)
          end do
          return
      end if

      del=0.5_wp*del
      do i=2,n
          angle=real(i-1,kind=wp)*del
          trigs(la+i)=sin(angle)
      end do

  end subroutine fftrig
    !
    !*****************************************************************************************
    !
  subroutine vpassm (a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,la)
      integer (ip), intent(in)  :: inc1
      integer (ip), intent(in)  :: inc2
      integer (ip), intent(in)  :: inc3
      integer (ip), intent(in)  :: inc4
      integer (ip), intent(in)  :: lot
      integer (ip), intent(in)  :: n
      integer (ip), intent(in)  :: ifac
      integer (ip), intent(in)  :: la
      real (wp),    intent(in)  :: a(n)
      real (wp),    intent(in)  :: b(n)
      real (wp),    intent(in)  :: trigs(n)
      real (wp),    intent(out) :: c(n)
      real (wp),    intent(out) :: d(n)
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

      real (wp), parameter :: SIN_36 = 0.587785252292473129168705954639072768597652437643145991072_wp
      real (wp), parameter :: COS_36 = 0.809016994374947424102293417182819058860154589902881431067_wp
      real (wp), parameter :: SIN_72 = 0.951056516295153572116439333379382143405698634125750222447_wp
      real (wp), parameter :: COS_72 = 0.309016994374947424102293417182819058860154589902881431067_wp
      real (wp), parameter :: SIN_60 = 0.866025403784438646763723170752936183471402626905190314027_wp

      integer (ip) :: i
      integer (ip) :: j
      integer (ip) :: k
      integer (ip) :: l
      integer (ip) :: m
      integer (ip) :: iink
      integer (ip) :: jink
      integer (ip) :: jump
      integer (ip) :: ibase
      integer (ip) :: jbase
      integer (ip) :: igo
      integer (ip) :: ijk
      integer (ip) :: la1
      integer (ip) :: ia
      integer (ip) :: ja
      integer (ip) :: ib
      integer (ip) :: jb
      integer (ip) :: kb
      integer (ip) :: ic
      integer (ip) :: jc
      integer (ip) :: kc
      integer (ip) :: id
      integer (ip) :: jd
      integer (ip) :: kd
      integer (ip) :: ie
      integer (ip) :: je
      integer (ip) :: ke
      real (wp) :: c1
      real (wp) :: s1
      real (wp) :: c2
      real (wp) :: s2
      real (wp) :: c3
      real (wp) :: s3
      real (wp) :: c4
      real (wp) :: s4

      m=n/ifac
      iink=m*inc1
      jink=la*inc2
      jump=(ifac-1)*jink
      ibase=0
      jbase=0
      igo=ifac-1
      if (igo>4) return
      !del  goto (10,50,90,130),igo

      select case (igo)

          !   coding for factor 2

          case (1)
              ia=1
              ja=1
              ib=ia+iink
              jb=ja+jink
              do l=1,la
                  i=ibase
                  j=jbase
                  do ijk=1,lot
                      c(ja+j)=a(ia+i)+a(ib+i)
                      d(ja+j)=b(ia+i)+b(ib+i)
                      c(jb+j)=a(ia+i)-a(ib+i)
                      d(jb+j)=b(ia+i)-b(ib+i)
                      i=i+inc3
                      j=j+inc4

                  end do
                  ibase=ibase+inc1
                  jbase=jbase+inc2

              end do
              if (la == m) return
              la1=la+1
              jbase=jbase+jump
              do k=la1,m,la
                  kb=k+k-2
                  c1=trigs(kb+1)
                  s1=trigs(kb+2)
                  do l=1,la
                      i=ibase
                      j=jbase
                      do ijk=1,lot
                          c(ja+j)=a(ia+i)+a(ib+i)
                          d(ja+j)=b(ia+i)+b(ib+i)
                          c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
                          d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
                          i=i+inc3
                          j=j+inc4

                      end do
                      ibase=ibase+inc1
                      jbase=jbase+inc2

                  end do
                  jbase=jbase+jump

              end do
          !     return

          !   coding for factor 3

          case (2)
              ia=1
              ja=1
              ib=ia+iink
              jb=ja+jink
              ic=ib+iink
              jc=jb+jink
              do l=1,la
                  i=ibase
                  j=jbase
                  do ijk=1,lot
                      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
                      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
                      c(jb+j)=(a(ia+i)-0.5_wp*(a(ib+i)+a(ic+i)))-(SIN_60*(b(ib+i)-b(ic+i)))
                      c(jc+j)=(a(ia+i)-0.5_wp*(a(ib+i)+a(ic+i)))+(SIN_60*(b(ib+i)-b(ic+i)))
                      d(jb+j)=(b(ia+i)-0.5_wp*(b(ib+i)+b(ic+i)))+(SIN_60*(a(ib+i)-a(ic+i)))
                      d(jc+j)=(b(ia+i)-0.5_wp*(b(ib+i)+b(ic+i)))-(SIN_60*(a(ib+i)-a(ic+i)))
                      i=i+inc3
                      j=j+inc4
                  end do
                  ibase=ibase+inc1
                  jbase=jbase+inc2
              end do
              if (la == m) return
              la1=la+1
              jbase=jbase+jump
              do k=la1,m,la
                  kb=k+k-2
                  kc=kb+kb
                  c1=trigs(kb+1)
                  s1=trigs(kb+2)
                  c2=trigs(kc+1)
                  s2=trigs(kc+2)
                  do l=1,la
                      i=ibase
                      j=jbase
                      do ijk=1,lot
                          c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
                          d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
                          c(jb+j)=                                                           &
                              c1*((a(ia+i)-0.5_wp*(a(ib+i)+a(ic+i)))-(SIN_60*(b(ib+i)-b(ic+i)))) &
                              -s1*((b(ia+i)-0.5_wp*(b(ib+i)+b(ic+i)))+(SIN_60*(a(ib+i)-a(ic+i))))
                          d(jb+j)=                                                           &
                              s1*((a(ia+i)-0.5_wp*(a(ib+i)+a(ic+i)))-(SIN_60*(b(ib+i)-b(ic+i)))) &
                              +c1*((b(ia+i)-0.5_wp*(b(ib+i)+b(ic+i)))+(SIN_60*(a(ib+i)-a(ic+i))))
                          c(jc+j)=                                                           &
                              c2*((a(ia+i)-0.5_wp*(a(ib+i)+a(ic+i)))+(SIN_60*(b(ib+i)-b(ic+i)))) &
                              -s2*((b(ia+i)-0.5_wp*(b(ib+i)+b(ic+i)))-(SIN_60*(a(ib+i)-a(ic+i))))
                          d(jc+j)=                                                           &
                              s2*((a(ia+i)-0.5_wp*(a(ib+i)+a(ic+i)))+(SIN_60*(b(ib+i)-b(ic+i)))) &
                              +c2*((b(ia+i)-0.5_wp*(b(ib+i)+b(ic+i)))-(SIN_60*(a(ib+i)-a(ic+i))))
                          i=i+inc3
                          j=j+inc4
                      end do
                      ibase=ibase+inc1
                      jbase=jbase+inc2
                  end do
                  jbase=jbase+jump
              end do
          !     return

          !   coding for factor 4

          case (3)
              ia=1
              ja=1
              ib=ia+iink
              jb=ja+jink
              ic=ib+iink
              jc=jb+jink
              id=ic+iink
              jd=jc+jink
              do l=1,la
                  i=ibase
                  j=jbase
                  do ijk=1,lot
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
                  end do
                  ibase=ibase+inc1
                  jbase=jbase+inc2
              end do
              if (la == m) return
              la1=la+1
              jbase=jbase+jump
              do k=la1,m,la
                  kb=k+k-2
                  kc=kb+kb
                  kd=kc+kb
                  c1=trigs(kb+1)
                  s1=trigs(kb+2)
                  c2=trigs(kc+1)
                  s2=trigs(kc+2)
                  c3=trigs(kd+1)
                  s3=trigs(kd+2)
                  do l=1,la
                      i=ibase
                      j=jbase
                      do ijk=1,lot
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
                      end do
                      ibase=ibase+inc1
                      jbase=jbase+inc2
                  end do
                  jbase=jbase+jump
              end do
          !     return

          !   coding for factor 5

          case (4)
              ia=1
              ja=1
              ib=ia+iink
              jb=ja+jink
              ic=ib+iink
              jc=jb+jink
              id=ic+iink
              jd=jc+jink
              ie=id+iink
              je=jd+jink
              do l=1,la
                  i=ibase
                  j=jbase
                  do ijk=1,lot
                      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
                      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
                      c(jb+j)=(a(ia+i)+COS_72*(a(ib+i)+a(ie+i))-COS_36*(a(ic+i)+a(id+i))) &
                          -(SIN_72*(b(ib+i)-b(ie+i))+SIN_36*(b(ic+i)-b(id+i)))
                      c(je+j)=(a(ia+i)+COS_72*(a(ib+i)+a(ie+i))-COS_36*(a(ic+i)+a(id+i))) &
                          +(SIN_72*(b(ib+i)-b(ie+i))+SIN_36*(b(ic+i)-b(id+i)))
                      d(jb+j)=(b(ia+i)+COS_72*(b(ib+i)+b(ie+i))-COS_36*(b(ic+i)+b(id+i))) &
                          +(SIN_72*(a(ib+i)-a(ie+i))+SIN_36*(a(ic+i)-a(id+i)))
                      d(je+j)=(b(ia+i)+COS_72*(b(ib+i)+b(ie+i))-COS_36*(b(ic+i)+b(id+i))) &
                          -(SIN_72*(a(ib+i)-a(ie+i))+SIN_36*(a(ic+i)-a(id+i)))
                      c(jc+j)=(a(ia+i)-COS_36*(a(ib+i)+a(ie+i))+COS_72*(a(ic+i)+a(id+i))) &
                          -(SIN_36*(b(ib+i)-b(ie+i))-SIN_72*(b(ic+i)-b(id+i)))
                      c(jd+j)=(a(ia+i)-COS_36*(a(ib+i)+a(ie+i))+COS_72*(a(ic+i)+a(id+i))) &
                          +(SIN_36*(b(ib+i)-b(ie+i))-SIN_72*(b(ic+i)-b(id+i)))
                      d(jc+j)=(b(ia+i)-COS_36*(b(ib+i)+b(ie+i))+COS_72*(b(ic+i)+b(id+i))) &
                          +(SIN_36*(a(ib+i)-a(ie+i))-SIN_72*(a(ic+i)-a(id+i)))
                      d(jd+j)=(b(ia+i)-COS_36*(b(ib+i)+b(ie+i))+COS_72*(b(ic+i)+b(id+i))) &
                          -(SIN_36*(a(ib+i)-a(ie+i))-SIN_72*(a(ic+i)-a(id+i)))
                      i=i+inc3
                      j=j+inc4
                  end do
                  ibase=ibase+inc1
                  jbase=jbase+inc2
              end do
              if (la == m) return
              la1=la+1
              jbase=jbase+jump
              do k=la1,m,la
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
                  do l=1,la
                      i=ibase
                      j=jbase
                      do ijk=1,lot
                          c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
                          d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
                          c(jb+j)=                                                          &
                              c1*((a(ia+i)+COS_72*(a(ib+i)+a(ie+i))-COS_36*(a(ic+i)+a(id+i))) &
                              -(SIN_72*(b(ib+i)-b(ie+i))+SIN_36*(b(ic+i)-b(id+i))))         &
                              -s1*((b(ia+i)+COS_72*(b(ib+i)+b(ie+i))-COS_36*(b(ic+i)+b(id+i))) &
                              +(SIN_72*(a(ib+i)-a(ie+i))+SIN_36*(a(ic+i)-a(id+i))))
                          d(jb+j)=                                                          &
                              s1*((a(ia+i)+COS_72*(a(ib+i)+a(ie+i))-COS_36*(a(ic+i)+a(id+i))) &
                              -(SIN_72*(b(ib+i)-b(ie+i))+SIN_36*(b(ic+i)-b(id+i))))         &
                              +c1*((b(ia+i)+COS_72*(b(ib+i)+b(ie+i))-COS_36*(b(ic+i)+b(id+i))) &
                              +(SIN_72*(a(ib+i)-a(ie+i))+SIN_36*(a(ic+i)-a(id+i))))
                          c(je+j)=                                                          &
                              c4*((a(ia+i)+COS_72*(a(ib+i)+a(ie+i))-COS_36*(a(ic+i)+a(id+i))) &
                              +(SIN_72*(b(ib+i)-b(ie+i))+SIN_36*(b(ic+i)-b(id+i))))         &
                              -s4*((b(ia+i)+COS_72*(b(ib+i)+b(ie+i))-COS_36*(b(ic+i)+b(id+i))) &
                              -(SIN_72*(a(ib+i)-a(ie+i))+SIN_36*(a(ic+i)-a(id+i))))
                          d(je+j)=                                                          &
                              s4*((a(ia+i)+COS_72*(a(ib+i)+a(ie+i))-COS_36*(a(ic+i)+a(id+i))) &
                              +(SIN_72*(b(ib+i)-b(ie+i))+SIN_36*(b(ic+i)-b(id+i))))         &
                              +c4*((b(ia+i)+COS_72*(b(ib+i)+b(ie+i))-COS_36*(b(ic+i)+b(id+i))) &
                              -(SIN_72*(a(ib+i)-a(ie+i))+SIN_36*(a(ic+i)-a(id+i))))
                          c(jc+j)=                                                          &
                              c2*((a(ia+i)-COS_36*(a(ib+i)+a(ie+i))+COS_72*(a(ic+i)+a(id+i))) &
                              -(SIN_36*(b(ib+i)-b(ie+i))-SIN_72*(b(ic+i)-b(id+i))))         &
                              -s2*((b(ia+i)-COS_36*(b(ib+i)+b(ie+i))+COS_72*(b(ic+i)+b(id+i))) &
                              +(SIN_36*(a(ib+i)-a(ie+i))-SIN_72*(a(ic+i)-a(id+i))))
                          d(jc+j)=                                                          &
                              s2*((a(ia+i)-COS_36*(a(ib+i)+a(ie+i))+COS_72*(a(ic+i)+a(id+i))) &
                              -(SIN_36*(b(ib+i)-b(ie+i))-SIN_72*(b(ic+i)-b(id+i))))         &
                              +c2*((b(ia+i)-COS_36*(b(ib+i)+b(ie+i))+COS_72*(b(ic+i)+b(id+i))) &
                              +(SIN_36*(a(ib+i)-a(ie+i))-SIN_72*(a(ic+i)-a(id+i))))
                          c(jd+j)=                                                          &
                              c3*((a(ia+i)-COS_36*(a(ib+i)+a(ie+i))+COS_72*(a(ic+i)+a(id+i))) &
                              +(SIN_36*(b(ib+i)-b(ie+i))-SIN_72*(b(ic+i)-b(id+i))))         &
                              -s3*((b(ia+i)-COS_36*(b(ib+i)+b(ie+i))+COS_72*(b(ic+i)+b(id+i))) &
                              -(SIN_36*(a(ib+i)-a(ie+i))-SIN_72*(a(ic+i)-a(id+i))))
                          d(jd+j)=                                                          &
                              s3*((a(ia+i)-COS_36*(a(ib+i)+a(ie+i))+COS_72*(a(ic+i)+a(id+i))) &
                              +(SIN_36*(b(ib+i)-b(ie+i))-SIN_72*(b(ic+i)-b(id+i))))         &
                              +c3*((b(ia+i)-COS_36*(b(ib+i)+b(ie+i))+COS_72*(b(ic+i)+b(id+i))) &
                              -(SIN_36*(a(ib+i)-a(ie+i))-SIN_72*(a(ic+i)-a(id+i))))
                          i=i+inc3
                          j=j+inc4
                      end do
                      ibase=ibase+inc1
                      jbase=jbase+inc2
                  end do
                  jbase=jbase+jump
              end do

      end select

  end subroutine vpassm
      !
      !*****************************************************************************************
      !
  end module module_fast_fourier_transform
