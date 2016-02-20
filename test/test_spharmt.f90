program test_spharmt

    use, intrinsic :: iso_fortran_env, only: &
        wp     => REAL64, &
        ip     => INT32, &
        stderr => ERROR_UNIT

    use type_GaussianSphericalHarmonic
    !
    !  Test program for spharmt module - non-linear steady-state geostropic
    !  flow in a shallow water model.
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
    implicit none
    integer (ip), parameter :: nlon=128
    integer (ip), parameter :: nlat=nlon/2
    integer (ip), parameter :: ntrunc=42
    integer (ip), parameter :: nmdim = (ntrunc+1)*(ntrunc+2)/2
    complex, dimension(nmdim) :: vrtnm,divnm,pnm,scrnm
    complex, dimension(nmdim,3) :: dvrtdtnm,ddivdtnm,dpdtnm
    complex, dimension(ntrunc+1,nlat) :: scrm1,scrm2
    real, dimension(nlon,nlat) :: uxact,vxact,pxact,u,v,p,f, &
        coslat,ug,vg,pg,vrtg,divg,scrg1,scrg2
    integer itmax,mprint,nl,nlm1,nlm2,i,j,ncycle,&
        nsav1,nsav2,nold,nnow,nnew
    real lambda,lhat,phlt(361),uhat,aa,uzero,pzero,pi,hpi,dtr,omega,alphad,&
        alpha,fzero,dt,cfn,dlath,theta,sth,cth,ca,sa,dlam,st,ct,cthclh,cthslh,&
        clh,time,that,sl,slh,evmax,epmax,dvmax,dpmax,htime,dvgm,cl,v2max,p2max,&
        vmax,pmax
    type (GaussianSphericalHarmonic) :: this

    print *,'triangular trunction T',ntrunc
    print *,nlat,' gaussian latitudes'
    pi = 4.*atan(1.)
    hpi = pi/2.
    dtr = pi/180.
    aa = 6.37122e6
    omega = 7.292e-5
    fzero = omega+omega
    uzero = 40.
    pzero = 2.94e4
    alphad = 60.
    alpha = dtr*alphad

    dt = 300.
    itmax = nint(86400.*5./dt)
    mprint = itmax/10


    !     compute the derivative of the unrotated geopotential
    !     p as a function of latitude

    nl = 91
    nlm1 = nl-1
    nlm2 = nl-2
    cfn = 1./nlm1
    dlath = pi/nlm1
    do i=1,nlm2
        theta = i*dlath
        sth = sin(theta)
        cth = cos(theta)
        uhat = ui(uzero,hpi-theta)
        phlt(i) = cfn*cth*uhat*(uhat/sth+aa*fzero)
    enddo

    !     compute sine transform of the derivative of the geopotential
    !     for the purpose of computing the geopotential by integration
    !     see equation (3.9) in reference [1] above

    call sine(nlm2,phlt)

    !     compute the cosine coefficients of the unrotated geopotential
    !     by the formal integration of the sine series representation

    do i=1,nlm2
        phlt(i) = -phlt(i)/i
    enddo

    !     phlt(i) contains the coefficients in the cosine series
    !     representation of the unrotated geopotential that are used
    !     below to compute the geopotential on the rotated grid.
    !
    !     compute the initial values of  east longitudinal
    !     and latitudinal velocities u and v as well as the
    !     geopotential p and coriolis f on the rotated grid.
    !
    ca = cos(alpha)
    sa = sin(alpha)
    dlam = (pi+pi)/nlon

    !  initialize sphere derived data type.

    call initialize_gaussian_spherical_harmonic(this,nlon,nlat,ntrunc,aa)

    do j=1,nlon
        lambda = (j-1)*dlam
        cl = cos(lambda)
        sl = sin(lambda)
        do i=1,nlat

            !     lambda is longitude, theta is colatitude, and pi/2-theta is
            !     latitude on the rotated grid. lhat and that are longitude
            !     and colatitude on the unrotated grid. see text starting at
            !     equation (3.10)
            !
            theta = hpi-asin(this%gaulats(i))
            st = cos(theta)
            ct = sin(theta)
            sth = ca*st+sa*ct*cl
            cthclh = ca*ct*cl-sa*st
            cthslh = ct*sl
            lhat = atanxy(cthclh,cthslh)
            clh = cos(lhat)
            slh = sin(lhat)
            cth = clh*cthclh+slh*cthslh
            that = atanxy(sth,cth)
            uhat = ui(uzero,hpi-that)
            pxact(j,i) = cosine(that,nlm2,phlt)
            uxact(j,i) = uhat*(ca*sl*slh+cl*clh)
            vxact(j,i) = uhat*(ca*cl*slh*st-clh*sl*st+sa*slh*ct)
            f(j,i) = fzero*sth
            coslat(j,i) = sqrt(1.-this%gaulats(i)**2)
        enddo
    enddo

    vmax = 0.
    pmax = 0.
    v2max = 0.
    p2max = 0.
    do j=1,nlat
        do i=1,nlon
            v2max = v2max+uxact(i,j)**2+vxact(i,j)**2
            p2max = p2max+pxact(i,j)**2
            vmax = amax1(abs(uxact(i,j)),abs(vxact(i,j)),vmax)
            pmax = amax1(abs(pxact(i,j)),pmax)
        enddo
    enddo
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

    call getvrtdiv(this,vrtnm,divnm,ug,vg)
    call spharm(this,p,pnm,1)


    !==> time step loop.

    nnew = 1
    nnow = 2
    nold = 3

    do ncycle=0,itmax

        time = float(ncycle)*dt

        !==> INVERSE TRANSFORM TO GET VORT AND PHIG ON GRID.

        call spharm(this,pg,pnm,-1)
        call spharm(this,vrtg,vrtnm,-1)

        !==> compute u and v on grid from spectral coeffs. of vort and div.

        call getuv(this,vrtnm,divnm,ug,vg)

        !==> compute error statistics.

        if (mod(ncycle,mprint) == 0) then
            call spharm(this,divg,divnm,-1)
            u = ug/coslat
            v = vg/coslat
            p = pg
            htime = time/3600.
            write(*,390) ncycle,htime,dt,nlat,nlon,ntrunc,omega,pzero,uzero,alphad

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

            dvgm = 0.
            dvmax = 0.
            dpmax = 0.
            evmax = 0.0
            epmax = 0.0
            do j=1,nlat
                do i=1,nlon
                    dvgm = amax1(dvgm,abs(divg(i,j)))
                    dvmax = dvmax+(u(i,j)-uxact(i,j))**2+(v(i,j)-vxact(i,j))**2
                    dpmax = dpmax+(p(i,j)-pxact(i,j))**2
                    evmax = &
                        amax1(evmax,abs(v(i,j)-vxact(i,j)),abs(u(i,j)-uxact(i,j)))
                    epmax = amax1(epmax,abs(p(i,j)-pxact(i,j)))
                enddo
            enddo
            dvmax = sqrt(dvmax/v2max)
            dpmax = sqrt(dpmax/p2max)
            evmax = evmax/vmax
            epmax = epmax/pmax
            write(*,391) evmax,epmax,dvmax,dpmax,dvgm
391         format(' max error in velocity',1pe15.6,&
                ' max error in geopot. ',1pe15.6,/&
                ' l2 error in velocity ',1pe15.6,&
                ' l2 error in geopot.  ',1pe15.6,/&
                ' maximum divergence   ',1pe15.6)
        end if

        !==> COMPUTE RIGHT-HAND SIDES OF PROGNOSTIC EQNS.

        scrg1(:,:) = ug(:,:)*(vrtg(:,:)+f(:,:))
        scrg2(:,:) = vg(:,:)*(vrtg(:,:)+f(:,:))
        call perform_multiple_real_fft(this, scrg1, scrm1, 1)
        call perform_multiple_real_fft(this, scrg2, scrm2, 1)
        call sumnm(this,scrm1,scrm2,dvrtdtnm(:,nnew),-1,1)
        call sumnm(this,scrm2,scrm1,ddivdtnm(:,nnew),1,-1)
        scrg1(:,:) = ug(:,:)*(pg(:,:)+pzero)
        scrg2(:,:) = vg(:,:)*(pg(:,:)+pzero)
        call sumnm(this,scrm1,scrm2,dpdtnm(:,nnew),-1,1)
        scrg1(:,:)=pg(:,:)+0.5*((ug(:,:)**2+vg(:,:)**2)/coslat(:,:)**2)
        call spharm(this,scrg1,scrnm,1)
        ddivdtnm(:,nnew)=ddivdtnm(:,nnew)-this%lap(:)*scrnm(:)

        !==> update vrt and div with third-order adams-bashforth.

        !==> forward euler, then 2nd-order adams-bashforth time steps to start.

        if (ncycle == 0) then
            dvrtdtnm(:,nnow) = dvrtdtnm(:,nnew)
            dvrtdtnm(:,nold) = dvrtdtnm(:,nnew)
            ddivdtnm(:,nnow) = ddivdtnm(:,nnew)
            ddivdtnm(:,nold) = ddivdtnm(:,nnew)
            dpdtnm(:,nnow) = dpdtnm(:,nnew)
            dpdtnm(:,nold) = dpdtnm(:,nnew)
        else if (ncycle == 1) then
            dvrtdtnm(:,nold) = dvrtdtnm(:,nnew)
            ddivdtnm(:,nold) = ddivdtnm(:,nnew)
            dpdtnm(:,nold) = dpdtnm(:,nnew)
        end if
        vrtnm(:) = vrtnm(:) + dt*(&
            (23./12.)*dvrtdtnm(:,nnew) - (16./12.)*dvrtdtnm(:,nnow)+&
            (5./12.)*dvrtdtnm(:,nold) )
        divnm(:) = divnm(:) + dt*(&
            (23./12.)*ddivdtnm(:,nnew) - (16./12.)*ddivdtnm(:,nnow)+&
            (5./12.)*ddivdtnm(:,nold) )
        pnm(:) = pnm(:) + dt*(&
            (23./12.)*dpdtnm(:,nnew) - (16./12.)*dpdtnm(:,nnow)+&
            (5./12.)*dpdtnm(:,nold) )

        !==> SWITCH INDICES.

        nsav1 = nnew
        nsav2 = nnow
        nnew = nold
        nnow = nsav1
        nold = nsav2

    !==> end time step loop

    enddo

    !==> deallocate pointers in sphere object.

    call destroy_gaussian_spherical_harmonic(this)

contains

    real function ui(amp,thetad)
        !
        !     computes the initial unrotated longitudinal velocity
        !     see section 3.3.
        !
        real amp, thetad, thetab, thetae, xe, x
        thetab=-pi/6.
        thetae= pi/2.
        xe=3.e-1
        x =xe*(thetad-thetab)/(thetae-thetab)
        ui = 0.
        if(x<=0. .or. x>=xe) return
        ui=amp*exp(-1./x-1./(xe-x)+4./xe)
    end function ui

    real function atanxy(x,y)
        real x,y
        atanxy = 0.
        if(x==0. .and. y==0.) return
        atanxy = atan2(y,x)
    end function atanxy

    subroutine sine(n,x)
        !     computes the sine transform
        integer n,i,j
        real x(n),w(n),arg
        arg = 4.*atan(1.)/(n+1)
        do j=1,n
            w(j) = 0.
            do i=1,n
                w(j) = w(j)+x(i)*sin(i*j*arg)
            enddo
        enddo
        do i=1,n
            x(i) = 2.*w(i)
        enddo
    end subroutine sine

    real function cosine(theta,n,cf)
        !     computes the cosine transform
        integer n,i
        real cf(n), theta
        cosine = 0.
        do i=1,n
            cosine = cosine+cf(i)*cos(i*theta)
        enddo
    end function cosine

end program test_spharmt
