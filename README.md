# **gaussian\_spherical\_harmonic**
 
A object-oriented library in Fortran to perform spherical harmonic transforms on gaussian grids using triangular truncation. 

The original work by Jeff Whitaker <Jeffrey.S.Whitaker@noaa.gov>; written in Fortran 95, was heavily refactored to incorporate features of modern Fortran (2008+). More specifically, the former object-based **sphere** type, is rebaptized as the fully object-oriented class **GausianSphericalHarmonic**. This class encapsulates all the relevant subroutines as type-bound procedures and all previous instances of pointers are replaced with allocatable arrays to circumvent potential memory leaks.

-----------------------------------------------------------------------------

## Requirements
* The GNU Make tool https://www.gnu.org/software/make/
* The GNU gfortran compiler https://gcc.gnu.org/wiki/GFortran

-----------------------------------------------------------------------------

## To build the project

Type the following command line arguments
```
git clone https://github.com/jlokimlin/gaussian_spherical_harmonic.git

cd gaussian_spherical_harmonic; make all
```

-----------------------------------------------------------------------------

## Usage

```fortran

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    use type_GaussianSphericalHarmonic, only: &
        GaussianSphericalHarmonic

    ! Explicit typing only
    implicit none
    
    type (GaussianSphericalHarmonic) :: sphere
    integer (ip), parameter          :: NLON = 192
    integer (ip), parameter          :: NLAT = NLON/2
    integer (ip), parameter          :: NTRUNC = 62
    real (wp),    parameter          :: RSPHERE = 6.3712e+6_wp
    
    call sphere%create( NLON, NLAT, NTRUNC, RSPHERE )

```

By creating multiple instances of **GausianSphericalHarmonic**, spherical harmonic transforms on multiple grids can be done easily in the same code.

-----------------------------------------------------------------------------


## Results

```
		*** Test program for TYPE(GaussianSphericalHarmonic) ***
		
		Non-linear steady-state geostropic flow in a shallow water model
		
		Triangular trunction number  =          42
		Number of gaussian latitudes =          65
		
		
		 steady nonlinear rotated flow:
		 cycle number                       0 model time in  hours            0.00
		 time step in seconds            600. number of latitudes               65
		 number of longitudes             128 max wave number                   42
		 rotation rate           7.292000E-05 mean height             2.940000E+04
		 maximum velocity        4.000000E+01 tilt angle              6.000000E+01
		 max error in velocity   1.297794E-07 max error in geopot.    2.061491E-08
		 l2 error in velocity    6.884086E-08
		 l2 error in geopot.     3.012380E-09 maximum divergence      9.622449E-17
		
		
		 steady nonlinear rotated flow:
		 cycle number                      72 model time in  hours           12.00
		 time step in seconds            600. number of latitudes               65
		 number of longitudes             128 max wave number                   42
		 rotation rate           7.292000E-05 mean height             2.940000E+04
		 maximum velocity        4.000000E+01 tilt angle              6.000000E+01
		 max error in velocity   1.077176E-07 max error in geopot.    2.061491E-08
		 l2 error in velocity    1.111029E-07
		 l2 error in geopot.     3.012380E-09 maximum divergence      1.542544E-11
		
		
		 steady nonlinear rotated flow:
		 cycle number                     144 model time in  hours           24.00
		 time step in seconds            600. number of latitudes               65
		 number of longitudes             128 max wave number                   42
		 rotation rate           7.292000E-05 mean height             2.940000E+04
		 maximum velocity        4.000000E+01 tilt angle              6.000000E+01
		 max error in velocity   1.514993E-07 max error in geopot.    2.061491E-08
		 l2 error in velocity    1.418913E-07
		 l2 error in geopot.     3.012380E-09 maximum divergence      2.722109E-11
		
		
		 steady nonlinear rotated flow:
		 cycle number                     216 model time in  hours           36.00
		 time step in seconds            600. number of latitudes               65
		 number of longitudes             128 max wave number                   42
		 rotation rate           7.292000E-05 mean height             2.940000E+04
		 maximum velocity        4.000000E+01 tilt angle              6.000000E+01
		 max error in velocity   2.016065E-07 max error in geopot.    2.061491E-08
		 l2 error in velocity    1.694504E-07
		 l2 error in geopot.     3.012380E-09 maximum divergence      3.806707E-11
		
		
		 steady nonlinear rotated flow:
		 cycle number                     288 model time in  hours           48.00
		 time step in seconds            600. number of latitudes               65
		 number of longitudes             128 max wave number                   42
		 rotation rate           7.292000E-05 mean height             2.940000E+04
		 maximum velocity        4.000000E+01 tilt angle              6.000000E+01
		 max error in velocity   2.601369E-07 max error in geopot.    2.061491E-08
		 l2 error in velocity    1.927176E-07
		 l2 error in geopot.     3.012380E-09 maximum divergence      4.700744E-11
		
		
		 steady nonlinear rotated flow:
		 cycle number                     360 model time in  hours           60.00
		 time step in seconds            600. number of latitudes               65
		 number of longitudes             128 max wave number                   42
		 rotation rate           7.292000E-05 mean height             2.940000E+04
		 maximum velocity        4.000000E+01 tilt angle              6.000000E+01
		 max error in velocity   3.167848E-07 max error in geopot.    2.061491E-08
		 l2 error in velocity    2.146840E-07
		 l2 error in geopot.     3.012380E-09 maximum divergence      5.396183E-11
		
		
		 steady nonlinear rotated flow:
		 cycle number                     432 model time in  hours           72.00
		 time step in seconds            600. number of latitudes               65
		 number of longitudes             128 max wave number                   42
		 rotation rate           7.292000E-05 mean height             2.940000E+04
		 maximum velocity        4.000000E+01 tilt angle              6.000000E+01
		 max error in velocity   3.711835E-07 max error in geopot.    2.061491E-08
		 l2 error in velocity    2.347049E-07
		 l2 error in geopot.     3.012380E-09 maximum divergence      5.892684E-11
		
		
		 steady nonlinear rotated flow:
		 cycle number                     504 model time in  hours           84.00
		 time step in seconds            600. number of latitudes               65
		 number of longitudes             128 max wave number                   42
		 rotation rate           7.292000E-05 mean height             2.940000E+04
		 maximum velocity        4.000000E+01 tilt angle              6.000000E+01
		 max error in velocity   4.221550E-07 max error in geopot.    2.061491E-08
		 l2 error in velocity    2.539370E-07
		 l2 error in geopot.     3.012380E-09 maximum divergence      6.195340E-11
		
		
		 steady nonlinear rotated flow:
		 cycle number                     576 model time in  hours           96.00
		 time step in seconds            600. number of latitudes               65
		 number of longitudes             128 max wave number                   42
		 rotation rate           7.292000E-05 mean height             2.940000E+04
		 maximum velocity        4.000000E+01 tilt angle              6.000000E+01
		 max error in velocity   4.699382E-07 max error in geopot.    2.061491E-08
		 l2 error in velocity    2.717727E-07
		 l2 error in geopot.     3.012380E-09 maximum divergence      6.301328E-11
		
		
		 steady nonlinear rotated flow:
		 cycle number                     648 model time in  hours          108.00
		 time step in seconds            600. number of latitudes               65
		 number of longitudes             128 max wave number                   42
		 rotation rate           7.292000E-05 mean height             2.940000E+04
		 maximum velocity        4.000000E+01 tilt angle              6.000000E+01
		 max error in velocity   5.125485E-07 max error in geopot.    2.061491E-08
		 l2 error in velocity    2.890043E-07
		 l2 error in geopot.     3.012380E-09 maximum divergence      6.270651E-11
		
		
		 steady nonlinear rotated flow:
		 cycle number                     720 model time in  hours          120.00
		 time step in seconds            600. number of latitudes               65
		 number of longitudes             128 max wave number                   42
		 rotation rate           7.292000E-05 mean height             2.940000E+04
		 maximum velocity        4.000000E+01 tilt angle              6.000000E+01
		 max error in velocity   5.517727E-07 max error in geopot.    2.061491E-08
		 l2 error in velocity    3.043037E-07
		 l2 error in geopot.     3.012380E-09 maximum divergence      5.983927E-11
		
		
		This file was compiled by GCC version 5.3.1 
```


## Description of Class variables, i.e., type components are:

**NUMBER\_OF\_LONGITUDES** (integer) - number of longitudinal points.
 
**NUMBER\_OF\_LATITUDES** (integer) - number of latitudinal points.
 
**TRIANGULAR\_TRUNCATION\_LIMIT** (integer) - triangular truncation limit.
 
**RADIUS\_OF\_SPHERE** (real) - radius of sphere in meters.
 
**associated\_legendre\_functions** (real allocatable array dimension ```((ntrunc+1)*(ntrunc+2)/2, nlat)``` ) - Associated legendre polynomials.
 
**legendre\_derivative\_quantity** (real allocatable array, same size as *associated_legendre_functions*) = ```(-(1 - x**2)) (pnm/dx)```
 at ```x = sin(theta).```
 
**gaussian\_latitudes** (real allocatable array dimension ```nlat```) - ```sin(theta).```
 
**gaussian\_weights** (real allocatable array dimension ```nlat```) - gaussian weights.
 
**scaled\_gaussian\_weights** (real allocatable array dimension ```nlat```) - gaussian weights divided by ```Re(1-x**2)```.
 
**INDEX\_ORDER\_M** (integer allocatable array dimension (```(ntrunc+1)(ntrunc+2)/2```) - value of zonal wavenumber ```m``` corresponding to spectral array index ```nm=1,(ntrunc+1)*(ntrunc+2)/2)```
 
**INDEX\_DEGREE\_N** (integer allocatable array same size as *INDEX\_ORDER\_M*) - value of spherical harmonic degree ```n``` corresponding to spectral array index ```nm = 1,(ntrunc+1)*(ntrunc+2)/2)```
 
**laplacian** (real allocatable array dimension (```ntrunc+1)*(ntrunc+2)/2``` ) - lapacian operator in spectral space = ```-n(n+1)/Re**2```, where ```n``` is degree of spherical harmonic.
 
**inverse\_laplacian** (real allocatable array same size as *laplacian*) - inverse laplacian operator in spectral space.
 
**trigonometric\_functions**, **ifax**: private arrays needed by Temperton FFT.
 
**initialized** (logical) - true if instance of object has been initialized by call to create, false if not initialized.

-----------------------------------------------------------------------------

### Class methods, i.e., type-bound procedures are:

```fortran

     call sphere%create( nlon, nlat, ntrunc, re ) 
```

Initializes an object instance of **GausianSphericalHarmonic**. Inputs are ```nlon``` (number of unique longitudes), ```nlat``` (number of gaussian latitudes), and ```re``` (radius of sphere in meters). Must be called before anything else.
    
```fortran

    call sphere%destroy():
```
Cleans up allocatable arrays allocated by *create*.
 
```fortran
    
    call sphere%perform_spherical_harmonic_transform( ugrid, anm, idir )
```

Spherical harmonic transform (forward, i.e. grid to spectral, for ```idir=1``` and backward for ```idir=-1```). Arguments are gridded data ```ugrid```, complex spectral coefficients ```anm```, and flag specifying direction of transform (```idir```).  See **Import Details** below for information about indexing of grid and spectral arrays.

```fortran
    
     call sphere%get_latitudes_and_gaussian_weights( gaulats, weights )
```

Computes ```sin(gaussian latitudes)``` and ```gaussian weights```. Number of latitudes determined by size of input arrays.

```fortran
    
    call sphere%compute_associated_legendre_functions( x, pnm, hnm )
```
    
Computes associated legendre functions ```pnm``` and their derivates ```hnm = (X**2-1)*(dpnm/dx)```at ```x = sin(latitude)```. The input arrays pnm and hnm should have dimension ```(ntrunc+1)*(ntrunc+2)/2,``` where ```ntrunc``` is triangular truncation limit (see **Import Details** below for a description of the spectral indexing).

```fortran
    
    call sphere%get_velocities_from_vorticity_and_divergence( vrtnm, divnm, ug, vg )
```
    
Computes U,V (```u*cos(lat), v*cos(lat)``` on grid) from spectral coefficients of vorticity and divergence.
Input:  spectral coeffs. of vort and div ```vrtnm```, ```divnm```.
Output: gridded U,V ```(ug,vg)```.

```fortran
    
    call sphere%get_vorticity_and_divergence_from_velocities( vrtnm, divnm, ug, vg )
```
    
Computes spectral coefficients of vorticity and divergence ```vrtnm```,```divnm``` from gridded U,V ```(ug,vg)```.

```fortran
    
    call sphere%get_vector_gradient( chinm, uchig, vchig )
```
    
Compute ```cos(lat)*vector``` gradient.
Inputs: sphere, spectral coefficient array ```chinm```.
Outputs: longitudinal and latitudinal components of gradient on the gaussian grid ```(uchig,vchig)```.

```fortran
    
    call sphere%perform_isotropic_spectral_smoothing( datagrid, smooth )
```
    
Isotropic spectral smoothing.
Inputs: ```smooth(ntrunc+1)``` (smoothing factor as a function of spherical harmonic degree), ```datagrid``` (gridded data to be smoothed).
Outputs: ```datagrid``` (smoothed gridded data).

```fortran

    call sphere%get_complex_spherical_harmonic_coefficients( am, bm, anm, isign1, isign2 )
```
    
Given the arrays of fourier coeffs, ```am``` and ```bm```, computes the complex spherical harmonic coeffs ```anm``` of:
 ```isign1*( (1./rsphere*(1.-x**2))*d(ag)/d(lon) + (isign2/rsphere)*d(bg)/dx )``` where ```ag``` and ```bg``` are the grid point counterparts of ```am```, ```bm,``` ```isign1```, ```isign2``` are +1 or -1, ```rsphere``` is radius of sphere, ```x=sin(lat)```) ```am```, ```bm``` can be computed from gridded data ```(ag,bg)``` using *perform\_multiple\_real\_fft*.

```fortran

    call sphere%perform_multiple_real_fft( data, coeff, idir )
```
        
Computes fourier harmonics in zonal direction of a gridded array.  ```idir=+1 ``` for forward (grid to fourier) and -1 for backward (fourier to grid) transform.  ```data(nlon,nlat) ``` contains gridded data,  ```coeff(ntrunc+1,nlat) ``` contains complex zonal fourier harmonics.

-----------------------------------------------------------------------------

### Important Details:

 The gridded data is assumed to be oriented such that  ```i=1 ``` is the Greenwich
 meridian and  ```j=1 ``` is the northernmost point. Grid indices increase eastward
 and southward. The grid increment in longitude is  ```2*pi/nlon ``` radians.
 
 
 For example  ```nlon = 72 ``` for a five degree grid.  ```nlon ``` must be greater than or
 equal to 4. The efficiency of the computation is improved when  ```nlon ``` is a
 product of small prime numbers.

The spectral data is assumed to be in a complex array of dimension ```(NTRUNC+1)*(NTRUNC+2)/2. ``` ```NTRUNC ``` is the triangular truncation limit ( ```NTRUNC=42 ``` for  ```T42 ```).  ```NTRUN ``` must be  ```<= nlon/2. ``` Coefficients are ordered so that first  ```(nm=1) ``` is  ```m=0,n=0, ``` second is  ```m=0,n=1, nm=mtrunc ``` is  ```m=0,n=mtrunc, nm=mtrunc+1 ``` is  ```m=1,n=1, ``` etc. In Fortran syntax, values of  ```m ``` (degree) and  ```n ``` (order) corresponding as a function of the index  ```nm ``` are:

```fortran
    
    INDEX_ORDER_M = [ ((m,n=m,mtrunc),m=0,mtrunc) ]
    INDEX_DEGREE_N = [ ((n,n=m,mtrunc),m=0,mtrunc) ]
```

 Conversely, the index nm as a function of m and n is:

```fortran

     nm = sum([(i,i=mtrunc+1,mtrunc-m+2,-1)])+n-m+1
```

 The associated legendre polynomials are normalized so that the integral  ```pbar(n,m,theta)**2)*sin(theta) ``` on the interval  ```theta=0 ``` to  ```theta=pi ``` is 1, where:  ```pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))*sin(theta)m/(2**n*factorial(n)) ``` times the  ```(n+m) ```th derivative of ```(x**2-1)**n ``` with respect to  ```x=cos(theta) ```.

 note:  ```theta = 0.5*pi - phi, ``` where  ```phi ``` is latitude and  ```theta ``` is colatitude; hence,  ```cos(theta) = sin(phi) ``` and  ```sin(theta) = cos(phi) ```.

 Note that  ```pbar(0,0,theta)=sqrt(2)/2 ```, and  ```pbar(1,0,theta)=0.5*sqrt(6.)*sin(lat). ```
