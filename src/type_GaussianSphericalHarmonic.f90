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

    type, public :: GaussianSphericalHarmonic
        !---------------------------------------------------------------------------------
        ! Class variables
        !---------------------------------------------------------------------------------
        integer (ip),              public  :: NUMBER_OF_LONGITUDES = 0
        integer (ip),              public  :: NUMBER_OF_LATITUDES = 0
        integer (ip),              public  :: TRIANGULAR_TRUNCATION_LIMIT = 0
        integer (ip), allocatable, private :: indxm(:)
        integer (ip), allocatable, private :: indxn(:)
        integer (ip), allocatable, private :: ifax(:)
        real (wp),                 public  :: RADIUS_OF_SPHERE = 0.0_wp
        real (wp), allocatable,    private :: trigonometric_functions(:)
        real (wp), allocatable,    public  :: gaussian_latitudes(:)
        real (wp), allocatable,    private :: gaussian_weights(:)
        real (wp), allocatable,    private :: scaled_gaussian_weights(:)
        real (wp), allocatable,    public  :: laplacian(:)
        real (wp), allocatable,    private :: inverse_laplacian(:)
        real (wp), allocatable,    private :: associated_legendre_functions(:,:)
        real (wp), allocatable,    private :: legendre_derivative_quantity(:,:)
        logical,                   private :: initialized = .false.
        !---------------------------------------------------------------------------------
    contains
        !---------------------------------------------------------------------------------
        ! Class variables
        !---------------------------------------------------------------------------------
        procedure, public         :: create => initialize_gaussian_spherical_harmonic
        procedure, public         :: destroy => destroy_gaussian_spherical_harmonic
        procedure, public         :: perform_spherical_harmonic_transform
        procedure, public         :: perform_multiple_real_fft
        procedure, public         :: get_vector_gradient
        procedure, public         :: perform_isotropic_spectral_smoothing
        procedure, public         :: get_velocities_from_vorticity_and_divergence
        procedure, public         :: get_vorticity_and_divergence_from_velocities
        procedure, public         :: get_complex_spherical_harmonic_coefficients
        procedure, nopass, public :: get_latitudes_and_gaussian_weights
        procedure, nopass, public :: compute_associated_legendre_functions
        !---------------------------------------------------------------------------------
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
        this%NUMBER_OF_LONGITUDES = nlon
        this%NUMBER_OF_LATITUDES = nlat
        this%TRIANGULAR_TRUNCATION_LIMIT = ntrunc
        this%RADIUS_OF_SPHERE = re

        associate( nmdim => (ntrunc+1)*(ntrunc+2)/2 )

            !--------------------------------------------------------------------------------
            ! Allocate arrays
            !--------------------------------------------------------------------------------

            allocate(this%scaled_gaussian_weights(nlat))
            allocate(this%gaussian_latitudes(nlat))
            allocate(this%gaussian_weights(nlat))
            allocate(this%indxm(nmdim))
            allocate(this%indxn(nmdim))
            allocate(this%laplacian(nmdim))
            allocate(this%inverse_laplacian(nmdim))
            allocate(this%associated_legendre_functions(nmdim,nlat))
            allocate(this%legendre_derivative_quantity(nmdim,nlat))
            allocate(this%trigonometric_functions((3*nlon/2)+1))
            allocate(this%ifax(13))

            associate( &
                gaulats => this%gaussian_latitudes, &
                weights => this%gaussian_weights, &
                gwrc => this%scaled_gaussian_weights, &
                indxm => this%indxm, &
                indxn => this%indxn, &
                lap => this%laplacian, &
                ilap => this%inverse_laplacian, &
                pnm => this%associated_legendre_functions, &
                hnm => this%legendre_derivative_quantity, &
                trigs => this%trigonometric_functions, &
                ifax => this%ifax &
                )

                call get_latitudes_and_gaussian_weights( gaulats, weights )

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
        !
        ! Purpose:
        !
        ! Deallocate allocatables in object.
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (GaussianSphericalHarmonic), intent(in out) :: this
        !--------------------------------------------------------------------------------

        if (.not. this%initialized) return

        deallocate(this%gaussian_latitudes)
        deallocate(this%gaussian_weights)
        deallocate(this%scaled_gaussian_weights)
        deallocate(this%associated_legendre_functions)
        deallocate(this%legendre_derivative_quantity)
        deallocate(this%indxm)
        deallocate(this%indxn)
        deallocate(this%laplacian)
        deallocate(this%inverse_laplacian)
        deallocate(this%trigonometric_functions)
        deallocate(this%ifax)

        this%initialized = .false.

    end subroutine destroy_gaussian_spherical_harmonic
    !
    !*****************************************************************************************
    !
    subroutine get_latitudes_and_gaussian_weights(sinlats, wts)
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
        integer (ip) :: i, j !! Counters
        integer (ip) :: iteration_counter
        real (wp)    :: pp
        real (wp)    :: p1
        real (wp)    :: p2
        real (wp)    :: p3
        real (wp)    :: z
        real (wp)    :: z1
        logical      :: success = .false.
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

    end subroutine get_latitudes_and_gaussian_weights
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
        real (wp), intent(in out) :: pmn(:)
        real (wp), intent(in out) :: hmn(:)
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
        do i = 1, ntrunc + 1

            m = (i - 1)
            nmax = ntrunc + 1 - m

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
        !
        ! Purpose:
        !
        ! real multiple fft (uses temperton fft991)
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (GaussianSphericalHarmonic), intent(in)     :: this
        real (wp),                         intent(in out) :: data(:,:)
        complex (wp),                      intent(in out) :: coeff(:,:)
        integer (ip),                      intent(in)     :: idir
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        real (wp), allocatable :: input_output_data(:)
        real (wp), allocatable :: workspace(:)
        integer (ip)           :: i, j, m, n !! Counters
        !--------------------------------------------------------------------------------

        if (.not. this%initialized) then
            error stop 'uninitialized sphere object in perform_multiple_real_fft!'
        end if

        if ( this%TRIANGULAR_TRUNCATION_LIMIT > this%NUMBER_OF_LONGITUDES/2) then
            error stop 'ntrunc must be less than or equal to nlons in perform_multiple_real_fft'
        end if

        !--------------------------------------------------------------------------------
        ! Perform transform
        !--------------------------------------------------------------------------------

        associate( &
            nlons  => this%NUMBER_OF_LONGITUDES, &
            nlats  => this%NUMBER_OF_LATITUDES, &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            mwaves => this%TRIANGULAR_TRUNCATION_LIMIT + 1, &
            trigs  => this%trigonometric_functions, &
            ifax   => this%ifax &
            )

            !--------------------------------------------------------------------------------
            ! Allocate arrays
            !--------------------------------------------------------------------------------

            allocate( input_output_data( nlats * (nlons + 2) ) )
            allocate( workspace( nlats * (nlons + 1) ) )

            ! == > forward transform.

            select case (idir)
                case (+1)
    		
                    ! == > copy the data into the work array.
                    !    transforms are computed pairwise using a complex fft.
    		
                    n = 0
                    input_output_data = 0.0_wp
                    do j=1,nlats
                        do i=1,nlons+2
                            n = n + 1
                            input_output_data(n) = 0.0_wp
                            if (i <= nlons) then
                                input_output_data(n) = data(i,j)
                            end if
                        end do
                    end do
    		
                    call perform_fft991( &
                        input_output_data, workspace, trigs, ifax, 1, nlons+2, nlons, nlats, -1 )
    		
                    n = -1
                    do j=1,nlats
                        do m=1,(nlons/2)+1
                            n = n + 2
                            if (m <= mwaves) then
                                coeff(m,j) = &
                                    cmplx(input_output_data(n),input_output_data(n+1), kind=wp)
                            end if
                        end do
                    end do
                case (-1)
    		
                    input_output_data = 0.0_wp
                    n = -1
                    do j=1,nlats
                        do m=1,(nlons/2)+1
                            n = n + 2
                            if (m <= mwaves) then
                                input_output_data(n) = real(coeff(m,j), kind=wp)
                                input_output_data(n+1) = aimag(coeff(m,j))
                            end if
                        end do
                    end do
    		
                    call perform_fft991(input_output_data,workspace,trigs,ifax,1,nlons+2,nlons,nlats,1)
    		
                    n = 0
                    do j=1,nlats
                        do i=1,nlons+2
                            n = n + 1
                            if (i <= nlons) then
                                data(i,j) = input_output_data(n)
                            end if
                        end do
                    end do
                case default
                    error stop 'Calling argument idir must be +1 or -1 in perform_multiple_real_fft!'
            end select
        end associate

        !--------------------------------------------------------------------------------
        ! Free memory
        !--------------------------------------------------------------------------------

        deallocate( input_output_data )
        deallocate( workspace )

    end subroutine perform_multiple_real_fft
    !
    !*****************************************************************************************
    !
    subroutine perform_spherical_harmonic_transform(this, ugrid, anm, idir)

        ! Purpose:
        !
        ! Peforms the spherical harmonic transform
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (GaussianSphericalHarmonic), intent(in)      :: this
        real (wp),                         intent (in out) :: ugrid(:,:)
        complex (wp),                      intent (in out) :: anm(:)
        integer (ip),                      intent(in)      :: idir
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        complex (wp), allocatable  :: am(:,:)
        integer (ip)               :: nmstrt
        integer (ip)               :: nm !! index
        integer (ip)               :: m, n, j !! Counters
        !--------------------------------------------------------------------------------

        if (.not. this%initialized) then
            error stop 'uninitialized object in PERFORM_SPHERICAL_HARMONIC_TRANSFORM!'
        end if

        associate( &
            nlats   => this%NUMBER_OF_LATITUDES, &
            ntrunc  => this%TRIANGULAR_TRUNCATION_LIMIT, &
            mwaves  => this%TRIANGULAR_TRUNCATION_LIMIT + 1, &
            pnm     => this%associated_legendre_functions, &
            weights => this%gaussian_weights &
            )

            ! Allocate array
            allocate( am(mwaves, nlats) )

            select case (idir)
                case (+1)
    		
                    ! == >  grid space to spectral space transformation
                    !     first, initialize array.
    		
                    anm = 0.0_wp
    		
                    ! == > perform ffts on each latitude.
    		
                    call perform_multiple_real_fft(this, ugrid, am, 1)
    		
                    ! == >  sum over all gaussian latitudes for each mode and each wave to
                    !     obtain the transformed variable in spectral space.
    		
                    do j=1,nlats
                        nmstrt = 0
                        do m = 1, mwaves
                            do n = 1, mwaves-m+1
                                nm = nmstrt + n
                                anm(nm)=anm(nm)+pnm(nm,j)*weights(j)*am(m,j)
                            end do
                            nmstrt = nmstrt + mwaves-m+1
                        end do
                    end do
                case (-1)
    		
                    do j = 1, nlats
    		
                        ! == >  inverse legendre transform to get values of the zonal fourier
                        !     transform at latitude j.
    		
                        ! == >  sum the various meridional modes to build the fourier series
                        !     coefficient for zonal wavenumber m=i-1 at the given latitude.
    		
                        nmstrt = 0
                        do m = 1, mwaves
                            am(m,j) = cmplx(0.0_wp, 0.0_wp, kind=wp)
                            do n = 1, mwaves-m+1
                                nm = nmstrt + n
                                am(m,j) = am(m,j)  +  anm(nm) * pnm(nm,j)
                            end do
                            nmstrt = nmstrt + mwaves-m+1
                        end do
                    end do

                    ! == >  Fourier transform to compute the values of the variable in grid
                    !     space at the j-th latitude.
    		
                    call perform_multiple_real_fft(this, ugrid, am, -1)
                case default
                    error stop 'Calling argument idir must be -1 or +1!'
            end select
        end associate

        ! Free memory
        deallocate( am )

    end subroutine perform_spherical_harmonic_transform
    !
    !*****************************************************************************************
    !
    subroutine get_velocities_from_vorticity_and_divergence(this,vrtnm,divnm,ug,vg)
        !
        ! Purpose:
        !
        ! Computes U,V (winds times cos(lat)) from vrtnm,divnm
        ! (spectral coeffs of vorticity and divergence).
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (GaussianSphericalHarmonic), intent(in)  :: this
        real (wp),                         intent(out) :: ug(:,:)
        real (wp),                         intent(out) :: vg(:,:)
        complex (wp),                      intent(in)  :: vrtnm(:)
        complex (wp),                      intent(in)  :: divnm(:)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        complex (wp), allocatable :: um(:,:)
        complex (wp), allocatable :: vm(:,:)
        integer (ip)              :: n, m, j !! Counters
        integer (ip)              :: nm
        integer (ip)              :: nmstrt
        real (wp)                 :: rm
        !--------------------------------------------------------------------------------

        if (.not. this%initialized) then
            error stop 'uninitialized sphere object in getuv!'
        end if

        associate( &
            nlats  => this%NUMBER_OF_LATITUDES, &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            mwaves => this%TRIANGULAR_TRUNCATION_LIMIT + 1, &
            ilap   => this%inverse_laplacian, &
            re     => this%RADIUS_OF_SPHERE, &
            pnm    => this%associated_legendre_functions, &
            hnm    => this%legendre_derivative_quantity &
            )

            ! Allocate arrays
            allocate( um(mwaves, nlats) )
            allocate( vm(mwaves, nlats) )

            do j = 1, nlats
                nmstrt = 0
                do m = 1, mwaves
                    rm = real( m - 1, kind=wp)
                    um(m,j) = 0.0_wp
                    vm(m,j) = 0.0_wp
                    do n = 1, mwaves - m + 1

                        nm = nmstrt + n

                        um(m,j) = &
                            um(m,j) + (ilap(nm)/re)*( &
                            cmplx(0.0_wp,rm, kind=wp)*divnm(nm)*pnm(nm,j) + &
                            vrtnm(nm)*hnm(nm,j) )

                        vm(m,j) = &
                            vm(m,j) + (ilap(nm)/re)*( &
                            cmplx(0.0,rm, kind=wp) * vrtnm(nm)*pnm(nm,j) - &
                            divnm(nm)*hnm(nm,j) )

                    end do
                    nmstrt = nmstrt + mwaves - m + 1
                end do
            end do
        end associate

        call perform_multiple_real_fft(this, ug, um, -1)
        call perform_multiple_real_fft(this, vg, vm, -1)

        ! Free memory
        deallocate( um )
        deallocate( vm )

    end subroutine get_velocities_from_vorticity_and_divergence
    !
    !*****************************************************************************************
    !
    subroutine get_vorticity_and_divergence_from_velocities(this,vrtnm,divnm,ug,vg)
        !
        ! Purpose:
        !
        ! Computes vrtnm, divnm (spectral coeffs of vorticity and
        ! divergence) from U,V (winds time cos(lat)).
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (GaussianSphericalHarmonic), intent(in)     :: this
        real (wp),                         intent(in out) :: ug(:,:)
        real (wp),                         intent(in out) :: vg(:,:)
        complex (wp),                      intent(in)     :: vrtnm(:)
        complex (wp),                      intent(in)     :: divnm(:)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        complex (wp), allocatable :: um(:,:)
        complex (wp), allocatable :: vm(:,:)
        !--------------------------------------------------------------------------------

        if (.not. this%initialized) then
            error stop 'uninitialized sphere object in GET_VORTICITY_AND_DIVERGENCE_FROM_VELOCITIES!'
        end if

        ! Allocate arrays
        associate( &
            mwaves => this%TRIANGULAR_TRUNCATION_LIMIT + 1, &
            nlats  => this%NUMBER_OF_LATITUDES &
            )

            allocate( um( mwaves, nlats ) )
            allocate( vm( mwaves, nlats ) )

        end associate

        call perform_multiple_real_fft(this, ug, um, 1)
        call perform_multiple_real_fft(this, vg, vm, 1)

        call get_complex_spherical_harmonic_coefficients(this,um,vm,divnm,1,1)
        call get_complex_spherical_harmonic_coefficients(this,vm,um,vrtnm,1,-1)

        ! Free memory
        deallocate( um )
        deallocate( vm )

    end subroutine get_vorticity_and_divergence_from_velocities
    !
    !*****************************************************************************************
    !
    subroutine get_complex_spherical_harmonic_coefficients( this, &
        am, bm, anm, isign1, isign2 )
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
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (GaussianSphericalHarmonic), intent(in) :: this
        integer (ip),                      intent(in) :: isign1
        integer (ip),                      intent(in) :: isign2
        complex (wp)                                  :: anm(:)
        complex (wp),                      intent(in) :: am(:,:)
        complex (wp),                      intent(in) :: bm(:,:)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        integer (ip) :: j, n, m !! Counters
        integer (ip) :: nm
        integer (ip) :: nmstrt
        real (wp)    :: rm
        !--------------------------------------------------------------------------------

        if (.not. this%initialized) then
            error stop 'uninitialized sphere object in GET_COMPLEX_SPHERICAL_HARMONIC_COEFFICIENTS!'
        end if

        if (isign2 /= 1 .and. isign2 /= -1) then
            error stop 'isign2 must be +1 or -1 in GET_COMPLEX_SPHERICAL_HARMONIC_COEFFICIENTS!'
        end if

        if (isign1 /= 1 .and. isign1 /= -1) then
            error stop 'isign1 must be +1 or -1 in GET_COMPLEX_SPHERICAL_HARMONIC_COEFFICIENTS!'
        end if

        associate( &
            sign1  => real(isign1, kind=wp), &
            sign2  => real(isign2, kind=wp), &
            nlats  => this%NUMBER_OF_LATITUDES, &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            mwaves => this%TRIANGULAR_TRUNCATION_LIMIT+1, &
            gwrc   => this%scaled_gaussian_weights, &
            pnm    => this%associated_legendre_functions, &
            hnm    => this%legendre_derivative_quantity &
            )

            anm = 0.0_wp

            do j = 1,nlats
                nmstrt = 0
                do m = 1, mwaves
                    rm = real( m - 1, kind=wp )
                    do n   = 1, mwaves-m+1
                        nm = nmstrt + n
                        anm(nm) = &
                            anm(nm) &
                            + sign1 * gwrc(j) * (cmplx(0.0_wp,rm,kind=wp) &
                            * pnm(nm,j) * am(m,j) &
                            + sign2 * hnm(nm,j) * bm(m,j))
                    end do
                    nmstrt = nmstrt + mwaves - m +1
                end do
            end do

        end associate

    end subroutine get_complex_spherical_harmonic_coefficients
    !
    !*****************************************************************************************
    !
    subroutine get_vector_gradient( this, divnm, ug, vg )
        !
        ! Purpose:
        !
        ! Computes coslat * gradient of spectral coefficients (divnm)
        ! vector gradient returned on grid as (ug,vg)
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (GaussianSphericalHarmonic), intent(in)  :: this
        real (wp),                         intent(out) :: ug(:,:)
        real (wp),                         intent(out) :: vg(:,:)
        complex (wp),                      intent(in)  :: divnm(:)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        complex (wp), allocatable :: um(:,:)
        complex (wp), allocatable :: vm(:,:)
        integer (ip)              :: n, m, j !! Counters
        integer (ip)              :: nm
        integer (ip)              :: nmstrt
        real (wp)                 :: rm
        !--------------------------------------------------------------------------------

        if (.not. this%initialized) then
            error stop 'uninitialized sphere object in GET_VECTOR_GRADIENT!'
        end if

        associate( &
            nlats  => this%NUMBER_OF_LATITUDES, &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            mwaves => this%TRIANGULAR_TRUNCATION_LIMIT + 1, &
            re     => this%RADIUS_OF_SPHERE, &
            pnm    => this%associated_legendre_functions, &
            hnm    => this%legendre_derivative_quantity &
            )

            ! Allocate arrays
            allocate( um( mwaves, nlats ) )
            allocate( vm( mwaves, nlats ) )

            do j = 1, nlats
                nmstrt = 0
                do m = 1, mwaves
                    rm = real( m - 1, kind=wp)
                    um(m,j) = 0.0_wp
                    vm(m,j) = 0.0_wp
                    do n = 1, mwaves - m + 1
                        nm = nmstrt + n
                        um(m,j) = &
                            um(m,j) + (1.0_wp/re)* &
                            cmplx(0.0_wp,rm, kind=wp)*divnm(nm)*pnm(nm,j)

                        vm(m,j) = &
                            vm(m,j) - (1.0_wp/re)* &
                            divnm(nm)*hnm(nm,j)
                    end do
                    nmstrt = nmstrt + mwaves - m + 1
                end do
            end do

        end associate

        call perform_multiple_real_fft(this, ug, um, -1)
        call perform_multiple_real_fft(this, vg, vm, -1)

        ! Free memory
        deallocate( um )
        deallocate( vm )

    end subroutine get_vector_gradient
    !
    !*****************************************************************************************
    !
    subroutine perform_isotropic_spectral_smoothing(this,datagrid,smooth)
        !
        ! Purpose:
        !
        ! Performs isotropic spectral smoothing of datagrid.
        !
        ! Input: smooth(ntrunc+1) - smoothing factor as a function of degree (indxn).
        !
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (GaussianSphericalHarmonic), intent(in)     :: this
        real (wp),                         intent(in out) :: datagrid(:,:)
        real (wp),                         intent(in)     :: smooth(:)
        !--------------------------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------------------------
        complex (wp), allocatable :: dataspec(:)
        integer (ip)              :: n !! Counter
        integer (ip)              ::nm !! Index
        !--------------------------------------------------------------------------------

        if (.not. this%initialized) then
            error stop 'uninitialized object in PERFORM_ISOTROPIC_SPECTRAL_SMOOTHING!'
        end if

        associate( &
            ntrunc => this%TRIANGULAR_TRUNCATION_LIMIT, &
            nmdim  => (this%TRIANGULAR_TRUNCATION_LIMIT+1)*(this%TRIANGULAR_TRUNCATION_LIMIT+2)/2, &
            indxn  => this%indxn &
            )

            allocate( dataspec((ntrunc+1)*(ntrunc+2)/2) )

            call perform_spherical_harmonic_transform(this, datagrid, dataspec, 1)

            do nm =1, nmdim
                n = indxn(nm)
                dataspec(nm) = dataspec(nm)*smooth(n + 1)
            end do

            call perform_spherical_harmonic_transform(this, datagrid, dataspec, -1)

        end associate

        ! Free memory
        deallocate( dataspec )

    end subroutine perform_isotropic_spectral_smoothing
    !
    !*****************************************************************************************
    !
  end module type_GaussianSphericalHarmonic
