module type_GaussianSphericalHarmonic

    use, intrinsic :: iso_fortran_env, only: &
        wp     => REAL64, &
        ip     => INT32, &
        stderr => ERROR_UNIT

    use module_fast_fourier_transform, only: &
        fft991, &
        set99

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
    public :: sumnm,gaulw
    public :: legend

    type GaussianSphericalHarmonic

        ! rsphere is radius of sphere in m.
        real    :: rsphere = 0.

        ! nlons is number of longitudes (not including cyclic point).
        integer :: nlons = 0

        ! nlats is number of gaussian latitudes.
        integer :: nlats = 0

        ! for fft.
        real, dimension(:), allocatable :: trigs
        integer (ip), dimension(:), allocatable :: ifax

        ! ntrunc is triangular truncation limit (e.g. 42 for T42 truncation)
        integer :: ntrunc = 0

        ! gaulats is sin(gaussian lats), weights are gaussian weights,
        ! gwrc is gaussian weights divided by re*cos(lat)**2, lap is
        ! the laplacian operatore -n*(n+1)/re**2 (where n is the degree
        ! of the spherical harmonic),
        ! ilap is the inverse laplacian (1/lap, set to zero for n=0).
        real, dimension(:), allocatable :: gaulats, weights, gwrc, lap, ilap

        ! pnm is associated legendre polynomials, hnm is -(1-x**2)d(pnm)/dx,
        ! where x = gaulats.
        real, dimension(:,:), allocatable :: pnm, hnm

        ! indxm is zonal wavenumber index, indxn is spherical harmonic degree index.
        integer (ip), dimension(:), allocatable :: indxm, indxn
        logical :: initialized = .false.

    end type GaussianSphericalHarmonic

contains

    subroutine initialize_gaussian_spherical_harmonic(this,nlon,nlat,ntrunc,re)

        ! initialize a sphere object.

        integer (ip), intent(in) :: nlon
        integer (ip), intent(in) ::nlat
        integer (ip), intent(in) ::ntrunc
        real, intent(in) :: re
        class (GaussianSphericalHarmonic), intent(in out) :: this
        real, dimension(:), allocatable :: pnm_tmp
        real, dimension(:), allocatable ::hnm_tmp
        double precision, dimension(:), allocatable :: gaulats_tmp
        double precision, dimension(:), allocatable ::weights_tmp
        integer:: nmdim
        integer::m
        integer::n
        integer::j

        this%nlons = nlon
        this%nlats = nlat
        this%ntrunc = ntrunc
        this%rsphere = re

        nmdim = (ntrunc+1)*(ntrunc+2)/2

        allocate(gaulats_tmp(nlat))
        allocate(weights_tmp(nlat))
        call gaulw(gaulats_tmp,weights_tmp)
        allocate(this%gwrc(nlat))
        this%gwrc = weights_tmp/(dble(re)*(1.d0-gaulats_tmp**2))
        allocate(this%gaulats(nlat))
        allocate(this%weights(nlat))
        this%gaulats = gaulats_tmp
        this%weights = weights_tmp
        deallocate(weights_tmp)

        allocate(this%indxm(nmdim))
        allocate(this%indxn(nmdim))
        allocate(this%lap(nmdim))
        allocate(this%ilap(nmdim))

        this%indxm = (/((m,n=m,ntrunc),m=0,ntrunc)/)
        this%indxn = (/((n,n=m,ntrunc),m=0,ntrunc)/)
        this%lap(:)=-real(this%indxn(:))*real(this%indxn(:)+1)/re**2
        this%ilap(1) = 0.
        this%ilap(2:nmdim) = 1./this%lap(2:nmdim)

        allocate(this%pnm(nmdim,nlat))
        allocate(this%hnm(nmdim,nlat))
        allocate(pnm_tmp(nmdim))
        allocate(hnm_tmp(nmdim))
        do j=1,nlat
            call legend(gaulats_tmp(j),pnm_tmp,hnm_tmp)
            this%pnm(:,j) = pnm_tmp(:)
            this%hnm(:,j) = hnm_tmp(:)
        end do
        deallocate(gaulats_tmp)
        deallocate(pnm_tmp)
        deallocate(hnm_tmp)

        allocate(this%trigs((3*nlon/2)+1))
        allocate(this%ifax(13))
        call set99(this%trigs,this%ifax,nlon)

        this%initialized = .true.

    end subroutine initialize_gaussian_spherical_harmonic

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

    subroutine gaulw(sinlats, wts)

        ! compute sin of gaussian latitudes and gaussian weights.
        ! uses the iterative method presented in Numerical Recipes.

        double precision, intent (in out), dimension(:) :: sinlats
        double precision, intent (in out), dimension(:) :: wts

        integer:: itermax
        integer:: i
        integer:: iter
        integer:: j
        integer:: nlat
        integer:: nprec
        double precision:: pi
        double precision:: pp
        double precision:: p1
        double precision:: p2
        double precision:: p3
        double precision:: z
        double precision:: z1
        double precision:: converg
        double precision:: ten


        ten = 10.d0
        converg = ten * epsilon(ten)

        nprec = precision(converg)
        converg = .1**nprec

        itermax = 10

        pi = acos(-1.0)

        nlat = size(sinlats)
        if (size(sinlats) /= size(wts)) then
            print *, 'sinlats and wts must be same size in gaulw!'
            stop
        end if

        do i=1,nlat
            z = cos(pi*(i - 0.25)/(nlat + 0.5))
            do iter=1,itermax
                p1 = 1.0
                p2 = 0.0

                do j=1,nlat
                    p3 = p2
                    p2 = p1
                    p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j
                end do

                pp = nlat*(z*p1 - p2)/(z*z - 1.0e+00)
                z1 = z
                z  = z1 - p1/pp
                if(abs(z - z1) < converg) go to 10
            end do
            print *, 'abscissas failed to converge in itermax iterations'
            print *, 'stopping in gaulw!'
            stop

10      continue

        sinlats (i) = z
        wts (i) = 2.0/((1.0 - z*z)*pp*pp)

    end do

end subroutine gaulw

subroutine LEGend(x,pmn,hmn)
    !
    !     THIS SUBROUTINE COMPUTES ASSOCIATED LEGendRE
    !     FUNCTIONS, PMN, AND THE DERIVATIVE QUANTITY
    !     HMN = -(1 - X*X) * D(PMN)/DX  AT X = COS( COLATITUDE )
    !     GIVEN SOME COLATITUDE IN THE DOMAIN.
    !
    !     ACCURACY:
    !
    !     THE RECURSION RELATION EMPLOYED IS LESS ACCURATE
    !     NEAR THE POLES FOR M + N .GE. ROUGHLY 75 BECAUSE THE RECURSION
    !     INVOLVES THE DIFFERENCE OF NEARLY EQUAL NUMBERS, SO
    !     THE ASSOCIATED LEGendRE FUNCTIONS ARE COMPUTED USING DOUBLE
    !     PRECISION ACCUMULATORS.
    !     SEE BELOUSOV, TABLES OF NORMALIZED ASSOCIATED LEGendRE
    !     POLYNOMIALS, C. 1962, FOR INFORMATION ON THE ACCURATE
    !     COMPUTATION OF ASSOCIATED LEGendRE FUNCTIONS.

    real, dimension(:), intent(in out) ::  pmn
    real, dimension(:), intent(in out) ::  hmn
    double precision, intent(in) :: x
    integer:: m
    integer::n
    integer::nm
    integer::i
    integer::nmax
    integer::np1
    integer::nmstrt
    integer::j
    integer::nmdim
    integer::ntrunc
    double precision:: a
    double precision:: b
    double precision:: prod
    double precision:: sinsq
    double precision:: &
        eps
    double precision::epsnmp1
    double precision::epspmn
    double precision:: pmnj
    double precision:: pmnjm1
    double precision:: pmnjm2


    !**** SET PARAMETERS FOR ENTRY INTO THE RECURSIVE FORMULAE.


    sinsq = 1.D0 - x * x

    a     = 1.D0
    b     = 0.D0
    prod  = 1.D0

    nmdim = size(pmn)
    ntrunc = nint((-1.+sqrt(1+8*float(nmdim)))/2.)-1
    if (size(pmn) /= size(hmn)) then
        print *, 'pnm and hnm must be same size in subroutine legend!'
        stop
    end if

    !**** LOOP FOR THE 'M' INDEX.

    nmstrt = 0
    do i = 1, ntrunc+1

        m    = (i - 1)
        nmax = ntrunc+1-m

        !        WHEN M=0 (I=1), STANDARD LEGendRE POLYNOMIALS ARE
        !        GENERATED.

        if (m /= 0) then
            a    = a + 2.D0
            b    = b + 2.D0
            prod = prod * sinsq * a / b
        end if

        !****    GENERATE PMN AND HMN FOR J = 1 AND 2.

        pmnjm2   = SQRT(0.5D0 * prod)
        nm = nmstrt + 1
        pmn(nm) = pmnjm2

        pmnjm1   = SQRT( DBLE(2 * m + 3) ) * x * pmnjm2
        if (nm /= nmdim) pmn(nm+1) = pmnjm1

        np1 = m + 1
        epsnmp1 = SQRT( DBLE(np1*np1 - m*m) / DBLE(4*np1*np1 - 1) )
        epspmn   = x * pmnjm1 - epsnmp1 * pmnjm2

        hmn(nm) = DBLE(m) * epsnmp1 * pmnjm1
        if (nm /= nmdim) &
            hmn(nm+1) = DBLE(m+1) * epspmn  -  &
            DBLE(m+2) * epsnmp1 * pmnjm2

        !****    LOOP FOR THE 'N' INDEX.
        !        NOW APPLY THE RECURSION FORMULAE FOR J .GE. 3.

        do j    = 3, nmax
            n = m + j - 1
            nm = nmstrt + j
            eps = SQRT( DBLE(n*n - m*m) / DBLE(4*n*n - 1) )
            pmnj     = epspmn / eps
            pmn(nm) = pmnj

            !        COMPUTE EPS * PMN FOR J+1.

            epspmn   = x * pmnj - eps * pmnjm1

            !        COMPUTE THE DERIVATIVE.

            hmn(nm) = DBLE(n) * epspmn -  &
                DBLE(n+1) * eps * pmnjm1

            pmnjm2   = pmnjm1
            pmnjm1   = pmnj
        end do
        nmstrt = nmstrt + nmax
    end do

end subroutine LEGend

subroutine perform_multiple_real_fft(this, data, coeff, idir)

    ! real multiple fft (uses temperton fft991)

    class (GaussianSphericalHarmonic), intent(in) :: this

    real, dimension(this%nlons,this%nlats), intent(in out) :: data
    real, dimension(this%nlats*(this%nlons+2)) :: wrk1
    real, dimension(this%nlats*(this%nlons+1)) :: wrk2
    complex, dimension(this%ntrunc+1,this%nlats), intent(in out) :: coeff
    integer::  nlons
    integer::nlats
    integer::ntrunc
    integer::mwaves
    integer::i
    integer::j
    integer::m
    integer::n
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
            wrk1 = 0.
            do j=1,nlats
                do i=1,nlons+2
                    n = n + 1
                    wrk1(n) = 0.0
                    if (i <= nlons) then
                        wrk1(n) = data(i,j)
                    end if
                end do
            end do
    		
            call fft991(wrk1,wrk2,this%trigs,this%ifax,1,nlons+2,nlons,nlats,-1)
    		
            n = -1
            do j=1,nlats
                do m=1,(nlons/2)+1
                    n = n + 2
                    if (m <= mwaves) then
                        coeff(m,j) = cmplx(wrk1(n),wrk1(n+1))
                    end if
                end do
            end do
        case (-1)
    		
            wrk1 = 0.
            n = -1
            do j=1,nlats
                do m=1,(nlons/2)+1
                    n = n + 2
                    if (m <= mwaves) then
                        wrk1(n) = real(coeff(m,j))
                        wrk1(n+1) = aimag(coeff(m,j))
                    end if
                end do
            end do
    		
            call fft991(wrk1,wrk2,this%trigs,this%ifax,1,nlons+2,nlons,nlats,1)
    		
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

    real, dimension(this%nlons,this%nlats), intent(in out) :: ugrid
    complex, dimension((this%ntrunc+1)*(this%ntrunc+2)/2), intent(in out) :: anm
    complex, dimension(this%ntrunc+1,this%nlats) :: am
    integer::  nlats
    integer::ntrunc
    integer::mwaves
    integer::nmstrt
    integer::nm
    integer::m
    integer::n
    integer::j
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
                    am(m,j) = cmplx(0., 0.)
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
    real, dimension(this%nlons,this%nlats), intent(out) ::  ug
    real, dimension(this%nlons,this%nlats), intent(out) ::vg
    complex, dimension((this%ntrunc+1)*(this%ntrunc+2)/2), intent(in) :: vrtnm
    complex, dimension((this%ntrunc+1)*(this%ntrunc+2)/2), intent(in) ::divnm
    complex, dimension(this%ntrunc+1,this%nlats)  :: um
    complex, dimension(this%ntrunc+1,this%nlats)  ::vm
    integer:: nlats
    integer::ntrunc
    integer::mwaves
    integer::m
    integer::j
    integer::n
    integer::nm
    integer::nmstrt
    real::  rm

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
            um(m,j) = cmplx(0.,0.)
            vm(m,j) = cmplx(0.,0.)
            do n=1,mwaves-m+1
                nm = nmstrt + n
                um(m,j) = um(m,j) + (this%ilap(nm)/this%rsphere)*( &
                    cmplx(0.,rm)*divnm(nm)*this%pnm(nm,j) + &
                    vrtnm(nm)*this%hnm(nm,j) )
                vm(m,j) = vm(m,j) + (this%ilap(nm)/this%rsphere)*( &
                    cmplx(0.,rm)*vrtnm(nm)*this%pnm(nm,j) - &
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

    real, dimension(this%nlons,this%nlats), intent(in out) ::  ug
    real, dimension(this%nlons,this%nlats), intent(in out) ::vg
    complex, dimension((this%ntrunc+1)*(this%ntrunc+2)/2), intent(in) :: vrtnm
    complex, dimension((this%ntrunc+1)*(this%ntrunc+2)/2), intent(in) ::divnm
    complex, dimension(this%ntrunc+1,this%nlats) :: um
    complex, dimension(this%ntrunc+1,this%nlats) ::vm

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
    complex, dimension((this%ntrunc+1)*(this%ntrunc+2)/2) :: anm
    complex, dimension(this%ntrunc+1,this%nlats), intent(in) :: am
    complex, dimension(this%ntrunc+1,this%nlats), intent(in) ::bm
    integer:: nlats
    integer::ntrunc
    integer::mwaves
    integer::j
    integer::m
    integer::n
    integer::nm
    integer::nmstrt
    real::  sign1
    real::sign2
    real::rm

    if (.not. this%initialized) then
        print *, 'uninitialized sphere object in sumnm!'
        stop
    end if


    sign1 = float(isign1)
    sign2 = float(isign2)
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
                anm(nm) = aNM(nm) + sign1*this%GWrC(j)*(CMPLX(0.,rm) &
                    * this%PNM(nm,j) * am(m,j) &
                    + sign2 * this%HNM(nm,j) * bm(m,j))
            end do
            nmstrt = nmstrt + mwaves - m +1
        end do
    end do

end subroutine sumnm

subroutine cosgrad(this,divnm,ug,vg)

    ! compute coslat * gradient of spectral coefficients (divnm)
    ! vector gradient returned on grid as (ug,vg)

    class (GaussianSphericalHarmonic), intent(in) :: this

    real, dimension(this%nlons,this%nlats), intent(out) ::  ug
    real, dimension(this%nlons,this%nlats), intent(out) ::vg
    complex, dimension((this%ntrunc+1)*(this%ntrunc+2)/2), intent(in) :: divnm
    complex, dimension(this%ntrunc+1,this%nlats) :: um
    complex, dimension(this%ntrunc+1,this%nlats) ::vm
    integer:: nlats
    integer::ntrunc
    integer::mwaves
    integer::j
    integer::m
    integer::n
    integer::nm
    integer::nmstrt
    real:: rm

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
            um(m,j) = cmplx(0.,0.)
            vm(m,j) = cmplx(0.,0.)
            do n=1,mwaves-m+1
                nm = nmstrt + n
                um(m,j) = um(m,j) + (1./this%rsphere)* &
                    cmplx(0.,rm)*divnm(nm)*this%pnm(nm,j)
                vm(m,j) = vm(m,j) - (1./this%rsphere)* &
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
    real, dimension(this%nlons,this%nlats), intent(in out) :: datagrid
    real, dimension(this%ntrunc+1), intent(in) :: smooth
    complex, dimension((this%ntrunc+1)*(this%ntrunc+2)/2) :: dataspec
    integer:: n
    integer::nm
    integer::nmdim

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
