module type_GaussianSphericalHarmonic

    use, intrinsic :: iso_fortran_env, only: &
        wp     => REAL64, &
        ip     => INT32, &
        stderr => ERROR_UNIT

    use module_fast_fourier_transform, only: &
        perform_fft991, &
        initialize_fft99

    ! Explicit typing
    implicit none
    ! everything private to module, unless otherwise specified.

    private
    public :: GaussianSphericalHarmonic
    public :: initialize_gaussian_spherical_harmonic
    public :: destroy_gaussian_spherical_harmonic
    public :: spharm
    public :: perform_multiple_real_fft
    public :: cosgrad
    public :: specsmooth
    public :: getuv
    public :: getvrtdiv
    public :: sumnm
    public :: compute_latitudes_and_gaussian_weights
    public :: compute_associated_legendre_functions

    type, public :: GaussianSphericalHarmonic
        !---------------------------------------------------------------------------------
        ! Class variables
        !---------------------------------------------------------------------------------
        integer (ip),               public :: NLONS = 0
        integer (ip),               public :: NLATS = 0
        integer (ip),               public :: NTRUNC = 0
        integer (ip), allocatable, private :: indxm(:)
        integer (ip), allocatable, private :: indxn(:)
        integer (ip), allocatable, private :: ifax(:)
        real (wp),                  public  :: RADIUS_OF_SPHERE = 0.0_wp
        real (wp), allocatable, private :: trigs(:)
        real (wp), allocatable, public  :: gaulats(:)
        real (wp), allocatable, private :: weights(:)
        real (wp), allocatable, private :: gwrc(:)
        real (wp), allocatable, public  :: lap(:)
        real (wp), allocatable, private :: ilap(:)
        real (wp), allocatable, private :: pnm(:,:)
        real (wp), allocatable, private :: hnm(:,:)
        logical,                 private :: initialized = .false.
    end type GaussianSphericalHarmonic

contains
    !
    !*****************************************************************************************
    !
    subroutine initialize_gaussian_spherical_harmonic(this,nlon,nlat,ntrunc,re)
        !
        ! Purpose:
        !
        ! To initialize an GaussianSphericalHarmonic object
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (GaussianSphericalHarmonic), intent(in out) :: this
        integer (ip),                      intent(in)     :: nlon
        integer (ip),                      intent(in)     :: nlat
        integer (ip),                      intent(in)     :: ntrunc
        real (wp),                         intent(in)     :: re
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) ::m, n, j !! Counters
        !--------------------------------------------------------------------------------

        ! Set contants
        this%nlons = nlon
        this%nlats = nlat
        this%ntrunc = ntrunc
        this%RADIUS_OF_SPHERE = re

        associate( nmdim => (ntrunc+1)*(ntrunc+2)/2 )

            !--------------------------------------------------------------------------------
            ! Allocate arrays
            !--------------------------------------------------------------------------------

            allocate(this%gwrc(nlat))
            allocate(this%gaulats(nlat))
            allocate(this%weights(nlat))
            allocate(this%indxm(nmdim))
            allocate(this%indxn(nmdim))
            allocate(this%lap(nmdim))
            allocate(this%ilap(nmdim))
            allocate(this%pnm(nmdim,nlat))
            allocate(this%hnm(nmdim,nlat))
            allocate(this%trigs((3*nlon/2)+1))
            allocate(this%ifax(13))

            associate( &
                gaulats => this%gaulats, &
                weights => this%weights, &
                gwrc => this%gwrc, &
                indxm => this%indxm, &
                indxn => this%indxn, &
                lap => this%lap, &
                ilap => this%ilap, &
                pnm => this%pnm, &
                hnm => this%hnm, &
                trigs => this%trigs, &
                ifax => this%ifax &
                )

                call compute_latitudes_and_gaussian_weights( gaulats, weights )

                gwrc = weights/(re * (1.0_wp - (gaulats**2)))
                indxm = [ ((m, n = m,ntrunc), m = 0, ntrunc) ]
                indxn = [ ((n, n = m,ntrunc), m = 0, ntrunc) ]
                lap(:)= -real(indxn(:), kind=wp) * real(indxn(:) + 1, kind=wp) / (re**2)
                ilap(1) = 0.0_wp
                ilap(2:nmdim) = 1.0_wp/lap(2:nmdim)

                do j = 1, nlat
                    call compute_associated_legendre_functions( gaulats(j), pnm(:,j), hnm(:,j) )
                end do

                call initialize_fft99( trigs, ifax, nlon)

            end associate
        end associate

        this%initialized = .true.

    end subroutine initialize_gaussian_spherical_harmonic
    !
    !*****************************************************************************************
    !
    subroutine destroy_gaussian_spherical_harmonic(this)

        ! deallocate allocatables in sphere object.

        class (GaussianSphericalHarmonic), intent(in out) :: this

        if (.not. this%initialized) return

        deallocate(this%gaulats)
        deallocate(this%weights)
        deallocate(this%gwrc)
        deallocate(this%pnm)
        deallocate(this%hnm)
        deallocate(this%indxm)
        deallocate(this%indxn)
        deallocate(this%lap)
        deallocate(this%ilap)
        deallocate(this%trigs)
        deallocate(this%ifax)

    end subroutine destroy_gaussian_spherical_harmonic
    !
    !*****************************************************************************************
    !
    subroutine compute_latitudes_and_gaussian_weights(sinlats, wts)
        !
        ! Purpose:
        !
        ! Compute sin of gaussian latitudes and gaussian weights.
        ! uses the iterative method presented in Numerical Recipes.
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (wp), intent (in out) :: sinlats(:)
        real (wp), intent (in out) :: wts(:)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: i
        integer (ip) :: iteration_counter
        integer (ip) :: j
        real (wp):: pp
        real (wp):: p1
        real (wp):: p2
        real (wp):: p3
        real (wp):: z
        real (wp):: z1
        logical :: success
        !--------------------------------------------------------------------------------

        if (size(sinlats) /= size(wts)) then
            error stop 'sinlats and wts must be same size in gaulw!'
        end if

        ! Initialize convergence flag
        success = .false.

        associate( &
            MAXIMUM_NUMBER_OF_ITERATIONS => 10, &
            pi                          => acos(-1.0_wp), &
            nlat                        => size(sinlats), &
            TOLERANCE                   => 0.1_wp**(precision(10.0_wp * epsilon(10.0_wp))) &
            )


            do i=1,nlat
                z = cos(pi*(real(i, kind=wp) - 0.25_wp)/(real(nlat, kind=wp) + 0.5_wp))

                do iteration_counter = 1, MAXIMUM_NUMBER_OF_ITERATIONS
                    p1 = 1.0_wp
                    p2 = 0.0_wp

                    do j=1,nlat
                        p3 = p2
                        p2 = p1
                        p1 = ((2.0_wp*j - 1.0_wp)*z*p2 - (j - 1.0_wp)*p3)/j
                    end do

                    pp = nlat*(z*p1 - p2)/(z*z - 1.0e+0_wp)
                    z1 = z
                    z  = z1 - p1/pp
                    if( abs(z - z1) < TOLERANCE ) then
                        sinlats (i) = z
                        wts (i) = 2.0_wp/((1.0_wp - z*z)*(pp**2))
                        success = .true.
                        exit
                    else
                        success = .false.
                    end if

                end do
            end do
        end associate

        if ( .not. success ) then
            error stop 'abscissas failed to converge in COMPUTE_LATITUDES_AND_GAUSSIAN_WEIGHTS'
        end if

    end subroutine compute_latitudes_and_gaussian_weights
    !
    !*****************************************************************************************
    !
    subroutine compute_associated_legendre_functions(x, pmn, hmn)
        !
        ! Purpose:
        !
        ! Computes associated legendre functions, pmn, and the derivative quantity
        ! hnn = -(1 - x*x) * d(pmn)/dx  at x = cos( colatitude )
        ! given some colatitude in the domain.
        !
        ! accuracy:
        !
        !     the recursion relation employed is less accurate
        !     near the poles for m + n >= roughly 75 because the recursion
        !     involves the difference of nearly equal numbers, so
        !     the associated legendre functions are computed using double
        !     precision accumulators.
        !     see belousov, tables of normalized associated legendre
        !     polynomials, c. 1962, for information on the accurate
        !     computation of associated legendre functions.
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        real (wp), intent(in)     :: x
        real (wp), intent(in out) ::  pmn(:)
        real (wp), intent(in out) ::  hmn(:)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: m
        integer (ip) ::n
        integer (ip) ::nm
        integer (ip) ::i
        integer (ip) ::nmax
        integer (ip) ::np1
        integer (ip) ::nmstrt
        integer (ip) ::j
        integer (ip) ::nmdim
        integer (ip) ::ntrunc
        real (wp):: a
        real (wp):: b
        real (wp):: prod
        real (wp):: sinsq
        real (wp):: eps
        real (wp):: epsnmp1
        real (wp):: epspmn
        real (wp):: pmnj
        real (wp):: pmnjm1
        real (wp):: pmnjm2


        !**** set parameters for entry into the recursive formulae.


        sinsq = 1.0_wp - x * x

        a     = 1.0_wp
        b     = 0.0_wp
        prod  = 1.0_wp

        nmdim = size(pmn)
        ntrunc = nint((-1.0_wp + sqrt(1.0_wp + 8.0_wp*real(nmdim, kind=wp)))/2.0_wp, kind=ip)-1

        if (size(pmn) /= size(hmn)) then
            error stop 'pnm and hnm must be same size in subroutine legend!'
        end if

        !**** loop for the 'm' index.

        nmstrt = 0
        do i = 1, ntrunc+1

            m = (i - 1)
            nmax = ntrunc+1-m

            !        when m=0 (i=1), standard legendre polynomials are
            !        generated.

            if (m /= 0) then
                a    = a + 2.0_wp
                b    = b + 2.0_wp
                prod = prod * sinsq * a / b
            end if

            !****    generate pmn and hmn for j = 1 and 2.

            pmnjm2   = sqrt(0.5_wp * prod)
            nm = nmstrt + 1
            pmn(nm) = pmnjm2

            pmnjm1   = sqrt( real(2 * m + 3, kind=wp) ) * x * pmnjm2
            if (nm /= nmdim) pmn(nm+1) = pmnjm1

            np1 = m + 1
            epsnmp1 = sqrt( real(np1*np1 - m*m, kind=wp) / real(4*np1*np1 - 1, kind=wp) )
            epspmn   = x * pmnjm1 - epsnmp1 * pmnjm2

            hmn(nm) = real(m, kind=wp) * epsnmp1 * pmnjm1
            if (nm /= nmdim) &
                hmn(nm+1) = real(m+1, kind=wp) * epspmn  -  &
                real(m+2, kind=wp) * epsnmp1 * pmnjm2

            !****    loop for the 'n' index.
            !        now apply the recursion formulae for j >= 3.

            do j = 3, nmax
                n = m + j - 1
                nm = nmstrt + j
                eps = sqrt( real(n*n - m*m, kind=wp) / real(4*n*n - 1, kind=wp) )
                pmnj = epspmn / eps
                pmn(nm) = pmnj

                !        compute eps * pmn for j+1.

                epspmn = x * pmnj - eps * pmnjm1

                !        compute the derivative.

                hmn(nm) = real( n, kind=wp) * epspmn -  &
                    real( n + 1, kind=wp) * eps * pmnjm1

                pmnjm2   = pmnjm1
                pmnjm1   = pmnj
            end do
            nmstrt = nmstrt + nmax
        end do

    end subroutine compute_associated_legendre_functions
    !
    !*****************************************************************************************
    !
    subroutine perform_multiple_real_fft(this, data, coeff, idir)

        ! real multiple fft (uses temperton fft991)

        class (GaussianSphericalHarmonic), intent(in) :: this

        real (wp), dimension(this%nlons,this%nlats), intent(in out) :: data
        real (wp), dimension(this%nlats*(this%nlons+2)) :: wrk1
        real (wp), dimension(this%nlats*(this%nlons+1)) :: wrk2
        complex (wp), intent(in out) :: coeff(this%ntrunc+1,this%nlats)
        integer (ip) ::  nlons
        integer (ip) ::nlats
        integer (ip) ::ntrunc
        integer (ip) ::mwaves
        integer (ip) ::i
        integer (ip) ::j
        integer (ip) ::m
        integer (ip) ::n
        integer (ip), intent(in) :: idir

        if (.not. this%initialized) then
            print *, 'uninitialized sphere object in perform_multiple_real_fft!'
            stop
        end if


        nlons = this%nlons
        nlats = this%nlats
        ntrunc = this%ntrunc
        if (ntrunc > nlons/2) then
            print *, 'ntrunc must be less than or equal to nlons in perform_multiple_real_fft'
            stop
        end if

        mwaves = ntrunc+1

        ! == > forward transform.

        select case (idir)
            case (+1)
    		
                ! == > copy the data into the work array.
                !    transforms are computed pairwise using a complex fft.
    		
                n = 0
                wrk1 = 0.0_wp
                do j=1,nlats
                    do i=1,nlons+2
                        n = n + 1
                        wrk1(n) = 0.0_wp
                        if (i <= nlons) then
                            wrk1(n) = data(i,j)
                        end if
                    end do
                end do
    		
                call perform_fft991(wrk1,wrk2,this%trigs,this%ifax,1,nlons+2,nlons,nlats,-1)
    		
                n = -1
                do j=1,nlats
                    do m=1,(nlons/2)+1
                        n = n + 2
                        if (m <= mwaves) then
                            coeff(m,j) = cmplx(wrk1(n),wrk1(n+1), kind=wp)
                        end if
                    end do
                end do
            case (-1)
    		
                wrk1 = 0.0_wp
                n = -1
                do j=1,nlats
                    do m=1,(nlons/2)+1
                        n = n + 2
                        if (m <= mwaves) then
                            wrk1(n) = real(coeff(m,j), kind=wp)
                            wrk1(n+1) = aimag(coeff(m,j))
                        end if
                    end do
                end do
    		
                call perform_fft991(wrk1,wrk2,this%trigs,this%ifax,1,nlons+2,nlons,nlats,1)
    		
                n = 0
                do j=1,nlats
                    do i=1,nlons+2
                        n = n + 1
                        if (i <= nlons) then
                            data(i,j) = wrk1(n)
                        end if
                    end do
                end do
            case default
                write(6,*) ' idir must be +1 or -1 in perform_multiple_real_fft!'
                write(6,*) ' execution terminated.'
                stop
        end select

    end subroutine perform_multiple_real_fft

    subroutine spharm(this, ugrid, anm, idir)

        ! spherical harmonic transform

        class (GaussianSphericalHarmonic), intent(in) :: this

        real (wp)    :: ugrid(this%nlons,this%nlats)
        complex (wp) :: anm((this%ntrunc+1)*(this%ntrunc+2)/2)
        complex (wp) :: am(this%ntrunc+1,this%nlats)
        integer (ip) :: nlats
        integer (ip) ::ntrunc
        integer (ip) ::mwaves
        integer (ip) ::nmstrt
        integer (ip) ::nm
        integer (ip) ::m
        integer (ip) ::n
        integer (ip) ::j
        integer (ip), intent(in) :: idir

        if (.not. this%initialized) then
            print *, 'uninitialized sphere object in spharm!'
            stop
        end if

        nlats = this%nlats
        ntrunc = this%ntrunc
        mwaves = ntrunc+1

        select case (idir)
            case (+1)
    		
                ! == >  GRID SPACE TO SPECTRAL SPACE TRANSFORMATION
                !     FIRST, INITIALIZE ARRAY.
    		
                anm = 0.
    		
                ! == > perform ffts on each latitude.
    		
                call perform_multiple_real_fft(this, ugrid, am, 1)
    		
                ! == >  SUM OVER ALL GAUSSIAN LATITUDES FOR EACH MODE AND EACH WAVE TO
                !     OBTAIN THE TRANSFORMED VARIABLE IN SPECTRAL SPACE.
    		
                do j=1,nlats
                    nmstrt = 0
                    do m = 1, mwaves
                        do n = 1, mwaves-m+1
                            nm = nmstrt + n
                            anm(nm)=anm(nm)+this%pnm(nm,j)*this%weights(j)*am(m,j)
                        end do
                        nmstrt = nmstrt + mwaves-m+1
                    end do
                end do
            case (-1)
    		
                do j = 1, nlats
    		
                    ! == >  INVERSE LEGendRE TRANSFORM TO GET VALUES OF THE ZONAL FOURIER
                    !     TRANSFORM AT LATITUDE j.
    		
                    ! == >  SUM THE VARIOUS MERIDIONAL MODES TO BUILD THE FOURIER SERIES
                    !     COEFFICIENT FOR ZONAL WAVENUMBER M=I-1 AT the GIVEN LATITUDE.
    		
                    nmstrt = 0
                    do m = 1, mwaves
                        am(m,j) = cmplx(0.0_wp, 0.0_wp, kind=wp)
                        do n = 1, mwaves-m+1
                            nm = nmstrt + n
                            am(m,j) = am(m,j)  +  anm(nm) * this%pnm(nm,j)
                        end do
                        nmstrt = nmstrt + mwaves-m+1
                    end do

                end do

                ! == >  FOURIER TRANSFORM TO COMPUTE THE VALUES OF THE VARIABLE IN GRID
                !     SPACE at THE J-TH LATITUDE.
    		
                call perform_multiple_real_fft(this, ugrid, am, -1)
            case default
                print *, 'error in spharm: idir must be -1 or +1!'
                print *, 'execution terminated in subroutine spharm'
        end select

    end subroutine spharm

    subroutine getuv(this,vrtnm,divnm,ug,vg)

        ! compute U,V (winds times cos(lat)) from vrtnm,divnm
        ! (spectral coeffs of vorticity and divergence).

        class (GaussianSphericalHarmonic), intent(in) :: this
        real (wp), dimension(this%nlons,this%nlats), intent(out) ::  ug
        real (wp), dimension(this%nlons,this%nlats), intent(out) ::vg
        complex (wp), intent(in) :: vrtnm((this%ntrunc+1)*(this%ntrunc+2)/2)
        complex (wp), intent(in) ::divnm((this%ntrunc+1)*(this%ntrunc+2)/2)
        complex (wp) :: um(this%ntrunc+1,this%nlats)
        complex (wp) ::vm(this%ntrunc+1,this%nlats)
        integer (ip) :: nlats
        integer (ip) ::ntrunc
        integer (ip) ::mwaves
        integer (ip) ::m
        integer (ip) ::j
        integer (ip) ::n
        integer (ip) ::nm
        integer (ip) ::nmstrt
        real (wp) ::  rm

        if (.not. this%initialized) then
            print *, 'uninitialized sphere object in getuv!'
            stop
        end if


        nlats = this%nlats
        ntrunc = this%ntrunc
        mwaves = ntrunc+1

        do j=1,nlats
            nmstrt = 0
            do m=1,mwaves
                rm = m-1
                um(m,j) = cmplx(0.0_wp,0.0_wp, kind=wp)
                vm(m,j) = cmplx(0.0_wp,0.0_wp, kind=wp)
                do n=1,mwaves-m+1
                    nm = nmstrt + n
                    um(m,j) = um(m,j) + (this%ilap(nm)/this%RADIUS_OF_SPHERE)*( &
                        cmplx(0.0_wp,rm, kind=wp)*divnm(nm)*this%pnm(nm,j) + &
                        vrtnm(nm)*this%hnm(nm,j) )
                    vm(m,j) = vm(m,j) + (this%ilap(nm)/this%RADIUS_OF_SPHERE)*( &
                        cmplx(0.0,rm, kind=wp)*vrtnm(nm)*this%pnm(nm,j) - &
                        divnm(nm)*this%hnm(nm,j) )
                end do
                nmstrt = nmstrt + mwaves-m+1
            end do
        end do

        call perform_multiple_real_fft(this, ug, um, -1)
        call perform_multiple_real_fft(this, vg, vm, -1)

    end subroutine getuv

    subroutine getvrtdiv(this,vrtnm,divnm,ug,vg)

        ! compute vrtnm,divnm (spectral coeffs of vorticity and
        ! divergence) from U,V (winds time cos(lat)).

        class (GaussianSphericalHarmonic), intent(in) :: this

        real (wp), dimension(this%nlons,this%nlats), intent(in out) ::  ug
        real (wp), dimension(this%nlons,this%nlats), intent(in out) ::vg
        complex (wp), dimension((this%ntrunc+1)*(this%ntrunc+2)/2), intent(in) :: vrtnm
        complex (wp), dimension((this%ntrunc+1)*(this%ntrunc+2)/2), intent(in) ::divnm
        complex (wp), dimension(this%ntrunc+1,this%nlats) :: um
        complex (wp), dimension(this%ntrunc+1,this%nlats) ::vm

        if (.not. this%initialized) then
            print *, 'uninitialized sphere object in getvrtdiv!'
            stop
        end if

        call perform_multiple_real_fft(this, ug, um, 1)
        call perform_multiple_real_fft(this, vg, vm, 1)

        call sumnm(this,um,vm,divnm,1,1)
        call sumnm(this,vm,um,vrtnm,1,-1)

    end subroutine getvrtdiv

    subroutine sumnm(this,am,bm,anm,isign1,isign2)
        !
        !  given the arrays of fourier coeffs, am and bm,
        !  compute the complex spherical harmonic coeffs of:
        !
        !  isign1*( (1./re*(1.-x**2))*d(ag)/d(lon) + (isign2/re)*d(bg)/dx )
        !
        !  where x = sin(lat), isign1 and isign2 are either +1 or -1,
        !  ag and bg are the physical space values of am,bm and re
        !  is the radius of the sphere.
        !
        !  the result is returned in anm.
        !
        !  for example on how to use this routine, see subroutine getvrtdiv.
        !
        !
        class (GaussianSphericalHarmonic), intent(in) :: this

        integer (ip), intent(in) :: isign1
        integer (ip), intent(in) ::isign2
        complex (wp), dimension((this%ntrunc+1)*(this%ntrunc+2)/2) :: anm
        complex (wp), dimension(this%ntrunc+1,this%nlats), intent(in) :: am
        complex (wp), dimension(this%ntrunc+1,this%nlats), intent(in) ::bm
        integer (ip) :: nlats
        integer (ip) ::ntrunc
        integer (ip) ::mwaves
        integer (ip) ::j
        integer (ip) ::m
        integer (ip) ::n
        integer (ip) ::nm
        integer (ip) ::nmstrt
        real (wp) ::  sign1
        real (wp) ::sign2
        real (wp) ::rm

        if (.not. this%initialized) then
            print *, 'uninitialized sphere object in sumnm!'
            stop
        end if


        sign1 = real(isign1, kind=wp)
        sign2 = real(isign2, kind=wp)
        if (isign2 /= 1 .and. isign2 /= -1) then
            print *, ' isign2 must be +1 or -1 in sumnm!'
            print *, ' execution terminated in sumnm'
        end if
        if (isign1 /= 1 .and. isign1 /= -1) then
            print *, ' isign1 must be +1 or -1 in sumnm!'
            print *, ' execution terminated in sumnm'
            stop
        end if
        nlats = this%nlats
        ntrunc = this%ntrunc
        mwaves = ntrunc+1

        anm = 0.
        do j=1,nlats
            nmstrt = 0
            do m = 1, mwaves
                rm = m-1
                do n   = 1, mwaves-m+1
                    nm = nmstrt + n
                    anm(nm) = anm(nm) + sign1*this%gwrc(j)*(cmplx(0.0_wp,rm,kind=wp) &
                        * this%pnm(nm,j) * am(m,j) &
                        + sign2 * this%hnm(nm,j) * bm(m,j))
                end do
                nmstrt = nmstrt + mwaves - m +1
            end do
        end do

    end subroutine sumnm

    subroutine cosgrad(this,divnm,ug,vg)

        ! compute coslat * gradient of spectral coefficients (divnm)
        ! vector gradient returned on grid as (ug,vg)

        class (GaussianSphericalHarmonic), intent(in) :: this

        real (wp), dimension(this%nlons,this%nlats), intent(out) ::  ug
        real (wp), dimension(this%nlons,this%nlats), intent(out) ::vg
        complex (wp), dimension((this%ntrunc+1)*(this%ntrunc+2)/2), intent(in) :: divnm
        complex (wp), dimension(this%ntrunc+1,this%nlats) :: um
        complex (wp), dimension(this%ntrunc+1,this%nlats) ::vm
        integer (ip) :: nlats
        integer (ip) ::ntrunc
        integer (ip) ::mwaves
        integer (ip) ::j
        integer (ip) ::m
        integer (ip) ::n
        integer (ip) ::nm
        integer (ip) ::nmstrt
        real (wp) :: rm

        if (.not. this%initialized) then
            print *, 'uninitialized sphere object in cosgrad!'
            stop
        end if


        nlats = this%nlats
        ntrunc = this%ntrunc
        mwaves = ntrunc+1

        do j=1,nlats
            nmstrt = 0
            do m=1,mwaves
                rm = (m-1)
                um(m,j) = cmplx(0.0_wp,0.0_wp, kind=wp)
                vm(m,j) = cmplx(0.0_wp,0.0_wp, kind=wp)
                do n=1,mwaves-m+1
                    nm = nmstrt + n
                    um(m,j) = um(m,j) + (1./this%RADIUS_OF_SPHERE)* &
                        cmplx(0.0_wp,rm, kind=wp)*divnm(nm)*this%pnm(nm,j)
                    vm(m,j) = vm(m,j) - (1./this%RADIUS_OF_SPHERE)* &
                        divnm(nm)*this%hnm(nm,j)
                end do
                nmstrt = nmstrt + mwaves - m +1
            end do
        end do

        call perform_multiple_real_fft(this, ug, um, -1)
        call perform_multiple_real_fft(this, vg, vm, -1)

    end subroutine cosgrad

    subroutine specsmooth(this,datagrid,smooth)

        ! isotropic spectral smoothing of datagrid.
        ! input: smooth(this%ntrunc+1) - smoothing factor as a
        ! function of degree (this%indxn).

        class (GaussianSphericalHarmonic), intent(in) :: this
        real (wp), dimension(this%nlons,this%nlats), intent(in out) :: datagrid
        real (wp), dimension(this%ntrunc+1), intent(in) :: smooth
        complex (wp), dimension((this%ntrunc+1)*(this%ntrunc+2)/2) :: dataspec
        integer (ip) :: n
        integer (ip) ::nm
        integer (ip) ::nmdim

        if (.not. this%initialized) then
            print *, 'uninitialized sphere object in specsmooth!'
            stop
        end if


        nmdim = (this%ntrunc+1)*(this%ntrunc+2)/2

        call spharm(this, datagrid, dataspec, 1)

        do nm=1,nmdim
            n = this%indxn(nm)
            dataspec(nm) = dataspec(nm)*smooth(n+1)
        end do

        call spharm(this, datagrid, dataspec, -1)

    end subroutine specsmooth

end module type_GaussianSphericalHarmonic
