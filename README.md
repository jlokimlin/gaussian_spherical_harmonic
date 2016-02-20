# Fortran 95 spherical harmonic module.
 
 For gaussian grids and triangular truncation only.

 version 1.1  9/30/2003 (second version - now all fortran using fft99)
 Jeff Whitaker <Jeffrey.S.Whitaker@noaa.gov>

-----------------------------------------------------------------------------
 To use this module, put "use spharmt" at top of main program
 (just after program statement).
 to initialize arrays needed to compute spherical harmonic transforms
 at T(ntrunc) resolution on a nlon x nlat gaussian grid,
 do something like this:

```fortran

 type (sphere) :: sphere_dat
 parameter (nlon=192,nlat=nlon/2)
 parameter (ntrunc=62)
 rsphere = 6.3712e6
 
 call spharmt_init(sphere_dat,nlon,nlat,ntrunc,re)


```

 derived data type "sphere" contains stuff needed for legendre
 transforms and fft.  By creating multiple instances of sphere
 data types, spherical harmonic transforms on multiple grids can
 be done easily in the same code.

### Components of sphere derived data type are:

 NLON (integer) - number of longitudes
 NLAT (integer) - number of latitudes
 NTRUNC (integer) - triangular truncation limit
 RE (real) - radius of sphere (m).
 PNM (real pointer array dimension ((ntrunc+1)*(ntrunc+2)/2, nlat) -
 Associated legendre polynomials.
 HNM (real pointer array, same size as pnm) = -(1 - x*x) * d(pnm)/dx
 at x = sin(latitude).
 GAULATS (real pointer array dimension nlat) - sin(gaussian latitudes).
 WEIGHTS (real pointer array dimension nlat) - gaussian weights.
 GWRC (real pointer array dimension nlat) - gaussian weights divided by
 rsphere*(1-x**2).
 INDXM (integer pointer array dimension (ntrunc+1)*(ntrunc+2)/2) - value of
 zonal wavenumber m corresponding to spectral array index
 nm=1,(ntrunc+1)*(ntrunc+2)/2)
 INDXN (integer pointer array same size as indxm) - value of
 spherical harmonic degree n corresponding to spectral array index
 nm=1,(ntrunc+1)*(ntrunc+2)/2)
 LAP (real pointer array dimension (ntrunc+1)*(ntrunc+2)/2) -
 lapacian operator in spectral space =
 -n*(n+1)/rsphere**2, where n is degree of spherical harmonic.
 ILAP (real pointer array same size as lap) - inverse laplacian operator in
 spectral space.
 TRIGS, IFAX: arrays needed by Temperton FFT.
 ISINITIALIZED (logical) - true if derived data type has been
 initialized by call to spharm_init, false if not initialized.

### Public routines:

 SPHARMT_INIT(sphere_dat,nlon,nlat,ntrunc,re):  initialize a sphere object
 (sphere_dat - derived data type "sphere"). inputs are nlon (number of unique
 longitudes), nlat (number of gaussian latitudes), and re
 (radius of sphere in m). Must be called before anything else.

 SPHARMT_DESTROY(sphere_dat):  cleans up pointer arrays allocated by
 spharmt_init.

 SPHARM(sphere_dat, ugrid, anm, idir):  spherical harmonic transform
 (forward, i.e. grid to spectral, for idir=1 and backward for idir=-1).
 Arguments are sphere derived data type (sphere_dat), gridded data (ugrid),
 complex spectral coefficients (anm), and flag specifying direction of
 transform (idir).  See "Import Details" below for information
 about indexing of grid and spectral arrays.

 GAULW(gaulats,weights): compute sin(gaussian latitudes) and gaussian
 weights.  Number of latitudes determined by size of input arrays.

 LEGEND(x,pnm,hnm): compute associated legendre functions (pnm) and
 their derivates hnm = (X**2-1)*D(PNM)/DX at x = sin(latitude). The
 input arrays pnm and hnm should have dimension (ntrunc+1)*(ntrunc+2)/2,
 where ntrunc is triangular truncation limit (see Import Details below
 for a description of the spectral indexing).

 GETUV(sphere_dat,vrtnm,divnm,ug,vg): get U,V (u*cos(lat),v*cos(lat) on grid)
 from spectral coefficients of vorticity and divergence.
 Input:  sphere_dat, spectral coeffs. of vort and div (vrtnm,divnm).
 Output: gridded U,V (ug,vg).

 GETVRTDIV(sphere_dat,vrtnm,divnm,ug,vg): get spectral coefficients of
 vorticity and divergence (vrtnm,divnm) from gridded U,V (ug,vg).

 COSGRAD(sphere_dat,chinm,uchig,vchig): compute cos(lat)*vector gradient.
 Inputs: sphere_dat, spectral coefficient array (chinm).
 Outputs: longitudinal and latitudinal components of gradient on the gaussian
 grid (uchig,vchig).

 SPECSMOOTH(sphere_dat,datagrid,smooth): isotropic spectral smoothing.
 Inputs: sphere_dat, smooth(ntrunc+1) (smoothing factor as a function
 of spherical harmonic degree), datagrid (gridded data to be smoothed).
 Outputs: datagrid (smoothed gridded data).

 SUMNM(sphere_dat,am,bm,anm,isign1,isign2):
 given the arrays of fourier coeffs, am and bm,
 compute the complex spherical harmonic coeffs (anm) of:
 isign1*( (1./rsphere*(1.-x**2))*d(ag)/d(lon) + (isign2/rsphere)*d(bg)/dx )
 where ag and bg are the grid pt counterparts of am,bm,
 isign1,isign2 are +1 or -1, rsphere is radius of sphere, x=sin(lat))
 am,bm can be computed from gridded data (ag,bg) using RFFT.

 RFFT(sphere_dat, data, coeff, idir): computes fourier harmonics in zonal
 direction of a gridded array. idir=+1 for forward (grid
 to fourier) and -1 for backward (fourier to grid) transform.
 data (nlon,nlat) contains gridded data, coeff (ntrunc+1,nlat) contains
 complex zonal fourier harmonics.

### Important Details:

 The gridded data is assumed to be oriented such that i=1 is the Greenwich
 meridian and j=1 is the northernmost point. Grid indices increase eastward
 and southward. The grid increment in longitude is 2*pi/nlon radians.
 For example nlon = 72 for a five degree grid. nlon must be greater than or
 equal to 4. The efficiency of the computation is improved when nlon is a
 product of small prime numbers.

 The spectral data is assumed to be in a complex array of dimension
 (NTRUNC+1)*(NTRUNC+2)/2. NTRUNC is the triangular truncation limit (NTRUNC=42
 for T42). NTRUNC must be <= nlon/2. Coefficients are ordered so that first
 (nm=1) is m=0,n=0, second is m=0,n=1, nm=mtrunc is m=0,n=mtrunc, nm=mtrunc+1
 is m=1,n=1, etc. In Fortran90 syntax, values of m (degree) and n (order)
 corresponding as a function of the index nm are:

 indxm = (/((m,n=m,mtrunc),m=0,mtrunc)/)
 indxn = (/((n,n=m,mtrunc),m=0,mtrunc)/)

 Conversely, the index nm as a function of m and n is:

 nm = sum((/(i,i=mtrunc+1,mtrunc-m+2,-1)/))+n-m+1

 The associated legendre polynomials are normalized so that the integral
 (pbar(n,m,theta)**2)*sin(theta) on the interval theta=0 to theta=pi is 1,
 where: pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
 *sin(theta)**m/(2**n*factorial(n)) times the (n+m)th derivative of
 (x**2-1)**n with respect to x=cos(theta).

 note: theta = 0.5*pi - phi, where phi is latitude and theta is colatitude,
 Therefore, cos(theta) = sin(phi) and sin(theta) = cos(phi).

 Note that pbar(0,0,theta)=sqrt(2)/2,
 and pbar(1,0,theta)=0.5*sqrt(6.)*sin(lat).
