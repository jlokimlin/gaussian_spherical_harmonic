!
!  Test program for GaussianSphericalHarmonic
!
!  Non-linear steady-state geostropic flow in a shallow water model.
!
!  errors should be O(10E-5) or less in single-precision, O(10E-7) or less
!  in real (wp).
!
!
!     the nonlinear shallow-water equations on the sphere are
!     solved using a spectral method based on the spherical harmonics.
!     the method is described in the paper:
!
! [1] p. n. swarztrauber, spectral transform methods for solving
!     the shallow-water equations on the sphere, p.n. swarztrauber,
!     monthly weather review, vol. 124, no. 4, april 1996, pp. 730-744.
!
!     this program implements test case 3 (steady nonlinear rotated flow)
!     in the paper:
!
! [2] d.l. williamson, j.b. drake, j.j. hack, r. jakob, and
!     p.n. swarztrauber, j. comp. phys., a standard test set
!     for numerical approximations to the shallow-water
!     equations in spherical geometry, j. comp. phys.,
!     vol. 102, no. 1, sept. 1992, pp. 211-224.
!
! definitions:
!
!
!     nlat          number of latitudes
!     nlon          number of distinct longitudes
!     ntrunc        max wave number
!     omega         rotation rate of earth in radians per second
!     aa            radius of earth in meters
!     pzero         mean height of geopotential
!     uzero         maximum velocity
!     alpha         tilt angle of the rotated grid
!     ncycle        cycle number
!     time          model time in seconds
!     dt            time step
!     lambda        longitude
!     theta         colatitude
!
!   the second dimension of the following two dimensional arrays
!   corresponds to the latitude index with values j=1,...,nlat
!   going from north to south.
!   the second dimension is longitude with values i=1,...,nlon
!   where i=1 corresponds to zero longitude and j=nlon corresponds
!   to 2pi minus 2pi/nlon.
!
!     u(i,j)       east longitudinal velocity component at t=time
!     v(i,j)       latitudinal velocity component at t=time
!     p(i,j)       +pzero = geopotential at t=time
!
!     divg(i,j)    divergence (d/dtheta (cos(theta) v)
!                                          + du/dlambda)/cos(theta)
!     vrtg(i,j)    vorticity  (d/dtheta (cos(theta) u)
!                                          - dv/dlambda)/cos(theta)
!
!     uxact(i,j)   the "exact" longitudinal velocity component
!     vxact(i,j)   the "exact" latitudinal  velocity component
!     pxact(i,j)   the "exact" geopotential
!
program test

    use, intrinsic :: iso_fortran_env, only: &
        wp     => REAL64, &
        ip     => INT32, &
        stdout => OUTPUT_UNIT

    use type_GaussianSphericalHarmonic

    ! Explicit typing only
    implicit none

    integer (ip), parameter :: nlon=128
    integer (ip), parameter :: nlat=nlon/2
    integer (ip), parameter :: ntrunc=42
    integer (ip), parameter :: nmdim = (ntrunc+1)*(ntrunc+2)/2
    complex (wp), dimension(nmdim) :: vrtnm
    complex (wp), dimension(nmdim) ::divnm
    complex (wp), dimension(nmdim) ::pnm
    complex (wp), dimension(nmdim) ::scrnm
    complex (wp), dimension(nmdim,3) :: dvrtdtnm
    complex (wp), dimension(nmdim,3) ::ddivdtnm
    complex (wp), dimension(nmdim,3) ::dpdtnm
    complex (wp), dimension(ntrunc+1,nlat) :: scrm1
    complex (wp), dimension(ntrunc+1,nlat) ::scrm2
    real (wp), dimension(nlon,nlat) :: uxact
    real (wp), dimension(nlon,nlat) ::vxact
    real (wp), dimension(nlon,nlat) ::pxact
    real (wp), dimension(nlon,nlat) ::u
    real (wp), dimension(nlon,nlat) ::v
    real (wp), dimension(nlon,nlat) ::p
    real (wp), dimension(nlon,nlat) ::f
    real (wp), dimension(nlon,nlat) :: &
        coslat
    real (wp), dimension(nlon,nlat) ::ug
    real (wp), dimension(nlon,nlat) ::vg
    real (wp), dimension(nlon,nlat) ::pg
    real (wp), dimension(nlon,nlat) ::vrtg
    real (wp), dimension(nlon,nlat) ::divg
    real (wp), dimension(nlon,nlat) ::scrg1
    real (wp), dimension(nlon,nlat) ::scrg2
    integer (ip)::  MAXIMUM_NUMBER_OF_TIME_ITERATIONS
    integer (ip)::mprint
    integer (ip)::nl
    integer (ip)::nlm1
    integer (ip)::nlm2
    integer (ip)::i
    integer (ip)::j
    integer (ip)::cycle_number
    integer (ip)::&
        nsav1
    integer (ip)::nsav2
    integer (ip)::n_old
    integer (ip)::n_now
    integer (ip)::n_new
    real (wp)::lhat
    real (wp)::phlt(361)
    real (wp)::uhat
    real (wp), parameter :: RADIUS_OF_EARTH_IN_METERS = 6.37122e+6_wp
    real (wp)::uzero
    real (wp)::pzero
    real (wp), parameter :: PI = acos( -1.0_wp )
    real (wp), parameter :: HALF_PI = 0.5_wp * PI
    real (wp)::dtr
    real (wp), parameter :: ROTATION_RATE_OF_EARTH = 7.292e-5_wp
    real (wp)::alphad
    real (wp)::&
        alpha
    real (wp)::fzero
    real (wp)::dt
    real (wp)::cfn
    real (wp)::LATITUDINAL_MESH
    real (wp)::sint
    real (wp)::cost
    real (wp)::cos_a
    real (wp)::sin_a
    real (wp)::LONGITUDINAL_MESH
    real (wp)::cos_t
    real (wp)::sin_t
    real (wp)::cthclh
    real (wp)::cthslh
    real (wp)::&
        cos_lh
    real (wp)::time
    real (wp)::that
    real (wp)::sin_l
    real (wp)::sin_lh
    real (wp)::evmax
    real (wp)::epmax
    real (wp)::dvmax
    real (wp)::dpmax
    real (wp)::htime
    real (wp)::dvgm
    real (wp)::cos_l
    real (wp)::v2max
    real (wp)::p2max
    real (wp)::vmax
    real (wp)::pmax
    type (GaussianSphericalHarmonic):: this

    write( stdout, '(A)') 'Test program for GaussianSphericalHarmonic'
    write( stdout, '(A)') ' '
    write( stdout, '(A)') 'Non-linear steady-state geostropic flow in a shallow water model'
    write( stdout, '(A)') ' '
    write( stdout, '(A,I11)') 'Triangular trunction number  = ', ntrunc
    write( stdout, '(A,I11)') 'Number of gaussian latitudes = ', nlat
    write( stdout, '(A)') ' '

    dtr = PI/180.0_wp
    fzero = 2.0_wp * ROTATION_RATE_OF_EARTH
    uzero = 40.0_wp
    pzero = 2.94e+4_wp
    alphad = 60.0_wp
    alpha = dtr*alphad

    dt = 300.0_wp
    MAXIMUM_NUMBER_OF_TIME_ITERATIONS = nint(864.0e+2_wp * 5.0_wp/dt, kind=ip)
    mprint = MAXIMUM_NUMBER_OF_TIME_ITERATIONS/10

    !     compute the derivative of the unrotated geopotential
    !     p as a function of latitude

    nl = 91
    nlm1 = nl-1
    nlm2 = nl-2
    cfn = 1.0_wp/nlm1
    LATITUDINAL_MESH = PI/nlm1
    do i=1,nlm2
        associate( theta => real(i, kind=wp) * LATITUDINAL_MESH )
            sint = sin(theta)
            cost = cos(theta)
            uhat = compute_initial_unrotated_longitudinal_velocity(uzero,HALF_PI-theta)
            phlt(i) = cfn*cost*uhat*(uhat/sint+RADIUS_OF_EARTH_IN_METERS*fzero)
        end associate
    end do

    !     compute sine transform of the derivative of the geopotential
    !     for the purpose of computing the geopotential by integration
    !     see equation (3.9) in reference [1] above

    call compute_sine_transform(phlt(1:nlm2))

    !     compute the cosine coefficients of the unrotated geopotential
    !     by the formal integration of the sine series representation

    do i=1,nlm2
        phlt(i) = -phlt(i)/i
    end do

    !     phlt(i) contains the coefficients in the cosine series
    !     representation of the unrotated geopotential that are used
    !     below to compute the geopotential on the rotated grid.
    !
    !     compute the initial values of  east longitudinal
    !     and latitudinal velocities u and v as well as the
    !     geopotential p and coriolis f on the rotated grid.
    !
    cos_a = cos(alpha)
    sin_a = sin(alpha)
    LONGITUDINAL_MESH = (2.0_wp * PI)/nlon

    !  initialize sphere derived data type.

    call initialize_gaussian_spherical_harmonic(this,nlon,nlat,ntrunc,RADIUS_OF_EARTH_IN_METERS)

    do j=1,nlon
        associate( lambda => real(j - 1, kind=wp) * LONGITUDINAL_MESH )
            cos_l = cos(lambda)
            sin_l = sin(lambda)
        end associate
        associate( gaulats => this%gaulats )
            do i = 1, nlat

                !     lambda is longitude, theta is colatitude, and pi/2-theta is
                !     latitude on the rotated grid. lhat and that are longitude
                !     and colatitude on the unrotated grid. see text starting at
                !     equation (3.10)
                !
                associate( theta => HALF_PI-asin(gaulats(i)) )
                    cos_t = cos(theta)
                    sin_t = sin(theta)
                    sint = cos_a*cos_t+sin_a*sin_t*cos_l
                    cthclh = cos_a*sin_t*cos_l-sin_a*cos_t
                    cthslh = sin_t*sin_l
                    lhat = atanxy(cthclh,cthslh)
                    cos_lh = cos(lhat)
                    sin_lh = sin(lhat)
                    cost = cos_lh*cthclh+sin_lh*cthslh
                    that = atanxy(sint,cost)
                    uhat = compute_initial_unrotated_longitudinal_velocity(uzero,HALF_PI-that)
                    pxact(j,i) = compute_cosine_transform(that,phlt)
                    uxact(j,i) = uhat*(cos_a*sin_l*sin_lh+cos_l*cos_lh)
                    vxact(j,i) = uhat*(cos_a*cos_l*sin_lh*cos_t-cos_lh*sin_l*cos_t+sin_a*sin_lh*sin_t)
                    f(j,i) = fzero * sint
                    coslat(j,i) = sqrt(1.0_wp - gaulats(i)**2)
                end associate
            end do
        end associate
    end do

    vmax = 0.0_wp
    pmax = 0.0_wp
    v2max = 0.0_wp
    p2max = 0.0_wp
    do j=1,nlat
        do i=1,nlon
            v2max = v2max+uxact(i,j)**2+vxact(i,j)**2
            p2max = p2max+pxact(i,j)**2
            vmax = max(abs(uxact(i,j)),abs(vxact(i,j)),vmax)
            pmax = max(abs(pxact(i,j)),pmax)
        end do
    end do
    !
    !     initialize first time step
    !
    u = uxact
    v = vxact
    p = pxact
    ug = u*coslat
    vg = v*coslat
    pg = p

    !  compute spectral coeffs of initial vrt,div,p.

    call get_vorticity_and_divergence_from_velocities(this,vrtnm,divnm,ug,vg)
    call perform_spherical_harmonic_transform(this,p,pnm,1)


    !==> time step loop.

    n_new = 1
    n_now = 2
    n_old = 3

    do cycle_number = 0, MAXIMUM_NUMBER_OF_TIME_ITERATIONS

        time = real(cycle_number, kind=wp)*dt

        !==> INVERSE TRANSFORM TO GET VORT AND PHIG ON GRID.

        call perform_spherical_harmonic_transform(this,pg,pnm,-1)
        call perform_spherical_harmonic_transform(this,vrtg,vrtnm,-1)

        !==> compute u and v on grid from spectral coeffs. of vort and div.

        call get_velocities_from_vorticity_and_divergence(this,vrtnm,divnm,ug,vg)

        !==> compute error statistics.

        if (mod(cycle_number,mprint) == 0) then
            call perform_spherical_harmonic_transform(this,divg,divnm,-1)
            u = ug/coslat
            v = vg/coslat
            p = pg
            htime = time/3600.0_wp
            write(stdout,390) cycle_number,htime,dt,nlat,nlon,ntrunc,ROTATION_RATE_OF_EARTH,pzero,uzero,alphad

390         format(//' steady nonlinear rotated flow:',/    &
                ' cycle number              ',i10,     &
                ' model time in  hours      ',f10.2,/  &
                ' time step in seconds      ',f10.0,   &
                ' number of latitudes       ',i10,/    &
                ' number of longitudes      ',i10,     &
                ' max wave number           ',i10,/    &
                ' rotation rate        ',1pe15.6,      &
                ' mean height          ',1pe15.6,/     &
                ' maximum velocity     ',1pe15.6,      &
                ' tilt angle           ',1pe15.6)

            dvgm = 0.0_wp
            dvmax = 0.0_wp
            dpmax = 0.0_wp
            evmax = 0.0_wp
            epmax = 0.0_wp

            do j=1,nlat
                do i=1,nlon
                    dvgm = &
                        max(dvgm,abs(divg(i,j)))
                    dvmax = &
                        dvmax+(u(i,j)-uxact(i,j))**2+(v(i,j)-vxact(i,j))**2
                    dpmax = &
                        dpmax+(p(i,j)-pxact(i,j))**2
                    evmax = &
                        max(evmax,abs(v(i,j)-vxact(i,j)),abs(u(i,j)-uxact(i,j)))
                    epmax = &
                        max(epmax,abs(p(i,j)-pxact(i,j)))
                end do
            end do

            dvmax = sqrt(dvmax/v2max)
            dpmax = sqrt(dpmax/p2max)
            evmax = evmax/vmax
            epmax = epmax/pmax

            write(stdout,391) evmax,epmax,dvmax,dpmax,dvgm
391         format(&
                ' max error in velocity',1pe15.6,&
                ' max error in geopot. ',1pe15.6,/&
                ' l2 error in velocity ',1pe15.6,&
                ' l2 error in geopot.  ',1pe15.6,/&
                ' maximum divergence   ',1pe15.6)
        end if

        !==> COMPUTE RIGHT-HAND SIDES OF PROGNOSTIC EQNS.

        scrg1 = ug * ( vrtg + f )
        scrg2 = vg * ( vrtg + f )

        call perform_multiple_real_fft(this, scrg1, scrm1, 1)
        call perform_multiple_real_fft(this, scrg2, scrm2, 1)

        call get_complex_spherical_harmonic_coefficients(this,scrm1,scrm2,dvrtdtnm(:,n_new),-1,1)
        call get_complex_spherical_harmonic_coefficients(this,scrm2,scrm1,ddivdtnm(:,n_new),1,-1)

        scrg1 = ug * ( pg + pzero )
        scrg2 = vg * ( pg + pzero )

        call get_complex_spherical_harmonic_coefficients( &
            this, scrm1, scrm2, dpdtnm(:,n_new), -1, 1)

        scrg1 = pg + 0.5_wp * ( ( ug**2 + vg**2 ) / coslat**2 )

        call perform_spherical_harmonic_transform( this, scrg1, scrnm, 1)

        associate( lap => this%lap )
            ddivdtnm(:,n_new) = ddivdtnm(:,n_new) - lap * scrnm
        end associate

        !==> update vrt and div with third-order adams-bashforth.

        !==> forward euler, then 2nd-order adams-bashforth time steps to start.

        select case (cycle_number)
            case (0)
                dvrtdtnm(:,n_now) = dvrtdtnm(:,n_new)
                dvrtdtnm(:,n_old) = dvrtdtnm(:,n_new)
                ddivdtnm(:,n_now) = ddivdtnm(:,n_new)
                ddivdtnm(:,n_old) = ddivdtnm(:,n_new)
                dpdtnm(:,n_now) = dpdtnm(:,n_new)
                dpdtnm(:,n_old) = dpdtnm(:,n_new)
            case (1)
                dvrtdtnm(:,n_old) = dvrtdtnm(:,n_new)
                ddivdtnm(:,n_old) = ddivdtnm(:,n_new)
                dpdtnm(:,n_old) = dpdtnm(:,n_new)
        end select

        vrtnm = &
            vrtnm + dt * (&
            (23.0_wp/12.0_wp) * dvrtdtnm(:,n_new) &
            - (16.0_wp/12.0_wp) * dvrtdtnm(:,n_now) &
            + (5.0_wp/12.0_wp) * dvrtdtnm(:,n_old) )

        divnm = &
            divnm + dt *( &
            (23.0_wp/12.0_wp) * ddivdtnm(:,n_new) &
            - (16.0_wp/12.0_wp) * ddivdtnm(:,n_now) &
            + (5.0_wp/12.0_wp) * ddivdtnm(:,n_old) )

        pnm = &
            pnm + dt * (&
            (23.0_wp/12.0_wp) * dpdtnm(:,n_new) &
            - (16.0_wp/12.0_wp) * dpdtnm(:,n_now) &
            + (5.0_wp/12.0_wp) * dpdtnm(:,n_old) )

        !==> SWITCH INDICES.

        nsav1 = n_new
        nsav2 = n_now
        n_new = n_old
        n_now = nsav1
        n_old = nsav2

    !==> end time step loop

    end do

    !==> deallocate arrays in object.

    call destroy_gaussian_spherical_harmonic(this)

contains

    pure function compute_initial_unrotated_longitudinal_velocity(amp,thetad) result (return_value)
        !
        !     computes the initial unrotated longitudinal velocity
        !     see section 3.3.
        !
        real (wp), intent (in) :: amp
        real (wp), intent (in) :: thetad
        real (wp)              :: return_value

        real (wp), parameter :: ZERO = nearest( 1.0_wp, 1.0_wp) - nearest( 1.0_wp, -1.0_wp)
        real (wp), parameter :: PI = acos( -1.0_wp )
        real (wp)            :: x

        associate( &
            thetab => -pi/6.0_wp, &
            thetae => pi/2.0_wp, &
            xe => 3.0e-1_wp &
            )

            x =xe*(thetad-thetab)/(thetae-thetab)

            return_value = 0.0_wp

            if(x <= ZERO .or. x >= xe) return

            associate( arg => (-1.0_wp/x) - (1.0_wp/(xe-x)) + (4.0_wp/xe) )
                return_value = amp * exp( arg )
            end associate

        end associate

    end function compute_initial_unrotated_longitudinal_velocity

    real (wp) function atanxy(x,y)
        real (wp):: x
        real (wp)::y
        atanxy = 0.
        if(x==0. .and. y==0.) return
        atanxy = atan2(y,x)
    end function atanxy

    subroutine compute_sine_transform( x )
        !
        real (wp), intent (in out) :: x(:)

        integer (ip)           :: i,j !! Counters
        real (wp), allocatable ::  w(:)

        associate( n => size(x) )

            allocate( w(n) )

            associate( arg => acos(-1.0_wp)/(n+1) )
                do j=1,n
                    w(j) = 0.0_wp
                    do i=1,n
                        associate( sin_arg => real(i*j,kind=wp)*arg )
                            w(j) = w(j)+x(i)*sin(sin_arg)
                        end associate
                    end do
                end do
            end associate
        end associate

        x = 2.0_wp * w

        ! Free memory
        deallocate(w)

    end subroutine compute_sine_transform

    pure function compute_cosine_transform(theta, cf) result (return_value)

        real (wp), intent(in) :: theta
        real (wp), intent(in) :: cf(:)
        real (wp)             :: return_value

        integer (ip)          :: i !! Counter

        return_value = 0.0_wp

        associate( n => size(cf) )
            do i=1,n
                return_value = return_value + cf(i)*cos(i*theta)
            end do
        end associate

    end function compute_cosine_transform

end program test
