
!IDEAL:MODEL_LAYER:INITIALIZATION
!

!  This MODULE holds the routines which are used to perform various initializations
!  for the individual domains.

!  This MODULE CONTAINS the following routines:

!  initialize_field_test - 1. Set different fields to different constant
!                             values.  This is only a test.  If the correct
!                             domain is not found (based upon the "id")
!                             then a fatal error is issued.

!-----------------------------------------------------------------------

MODULE module_initialize_ideal

   USE module_domain
   USE module_io_domain
   USE module_state_description
   USE module_model_constants
   USE module_bc
   USE module_timing
   USE module_configure
   USE module_init_utilities
   USE module_soil_pre
#ifdef DM_PARALLEL
   USE module_dm
#endif


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine normal_01_cdf ( x, cdf )

!*****************************************************************************80
!
!! NORMAL_01_CDF evaluates the Normal 01 CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    AG Adams,
!    Algorithm 39,
!    Areas Under the Normal Curve,
!    Computer Journal,
!    Volume 12, pages 197-198, 1969.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ), parameter :: a1 = 0.398942280444D+00
  real ( kind = 8 ), parameter :: a2 = 0.399903438504D+00
  real ( kind = 8 ), parameter :: a3 = 5.75885480458D+00
  real ( kind = 8 ), parameter :: a4 = 29.8213557808D+00
  real ( kind = 8 ), parameter :: a5 = 2.62433121679D+00
  real ( kind = 8 ), parameter :: a6 = 48.6959930692D+00
  real ( kind = 8 ), parameter :: a7 = 5.92885724438D+00
  real ( kind = 8 ), parameter :: b0 = 0.398942280385D+00
  real ( kind = 8 ), parameter :: b1 = 3.8052D-08
  real ( kind = 8 ), parameter :: b2 = 1.00000615302D+00
  real ( kind = 8 ), parameter :: b3 = 3.98064794D-04
  real ( kind = 8 ), parameter :: b4 = 1.98615381364D+00
  real ( kind = 8 ), parameter :: b5 = 0.151679116635D+00
  real ( kind = 8 ), parameter :: b6 = 5.29330324926D+00
  real ( kind = 8 ), parameter :: b7 = 4.8385912808D+00
  real ( kind = 8 ), parameter :: b8 = 15.1508972451D+00
  real ( kind = 8 ), parameter :: b9 = 0.742380924027D+00
  real ( kind = 8 ), parameter :: b10 = 30.789933034D+00
  real ( kind = 8 ), parameter :: b11 = 3.99019417011D+00
  real ( kind = 8 ) cdf
  real ( kind = 8 ) q
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  |X| <= 1.28.
!
  if ( abs ( x ) <= 1.28D+00 ) then

    y = 0.5D+00 * x * x

    q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
      + a6 / ( y + a7 ) ) ) )
!
!  1.28 < |X| <= 12.7
!
  else if ( abs ( x ) <= 12.7D+00 ) then

    y = 0.5D+00 * x * x

    q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
      + b2 / ( abs ( x ) + b3 &
      + b4 / ( abs ( x ) - b5 &
      + b6 / ( abs ( x ) + b7 &
      - b8 / ( abs ( x ) + b9 &
      + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
!
!  12.7 < |X|
!
  else

    q = 0.0D+00

  end if
!
!  Take account of negative X.
!
  if ( x < 0.0D+00 ) then
    cdf = q
  else
    cdf = 1.0D+00 - q
  end if

  return
end subroutine
function r8_uniform_01 ( seed ) result(value)

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 8 ) value
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  value = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end function
function r8poly_value_horner ( m, c, x ) result(value)

!*****************************************************************************80
!
!! R8POLY_VALUE_HORNER evaluates a polynomial using Horner's method.
!
!  Discussion:
!
!    The polynomial
!
!      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
!
!    is to be evaluated at the value X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the degree.
!
!    Input, real ( kind = 8 ) C(0:M), the polynomial coefficients.
!    C(I) is the coefficient of X^I.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) R8POLY_VALUE_HORNER, the polynomial value.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) c(0:m)
  integer ( kind = 4 ) i
  !real ( kind = 8 ) r8poly_value_horner
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  value = c(m)
  do i = m - 1, 0, -1
    value = value * x + c(i)
  end do

!  r8poly_value_horner = value

  return
end function
subroutine normal_01_cdf_inv ( p, x )

!*****************************************************************************80
!
!! NORMAL_01_CDF_INV inverts the standard normal CDF.
!
!  Discussion:
!
!    The result is accurate to about 1 part in 10^16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2015
!
!  Author:
!
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    Algorithm AS241:
!    The Percentage Points of the Normal Distribution,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the value of the cumulative probability
!    densitity function.  0 < P < 1.  If P is outside this range, an
!    "infinite" value will be returned.
!
!    Output, real ( kind = 8 ) X, the normal deviate value
!    with the property that the probability of a standard normal deviate being
!    less than or equal to the value is P.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 8 ) :: a = (/ &
    3.3871328727963666080D+00, &
    1.3314166789178437745D+02, &
    1.9715909503065514427D+03, &
    1.3731693765509461125D+04, &
    4.5921953931549871457D+04, &
    6.7265770927008700853D+04, &
    3.3430575583588128105D+04, &
    2.5090809287301226727D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: b = (/ &
    1.0D+00, &
    4.2313330701600911252D+01, &
    6.8718700749205790830D+02, &
    5.3941960214247511077D+03, &
    2.1213794301586595867D+04, &
    3.9307895800092710610D+04, &
    2.8729085735721942674D+04, &
    5.2264952788528545610D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: c = (/ &
    1.42343711074968357734D+00, &
    4.63033784615654529590D+00, &
    5.76949722146069140550D+00, &
    3.64784832476320460504D+00, &
    1.27045825245236838258D+00, &
    2.41780725177450611770D-01, &
    2.27238449892691845833D-02, &
    7.74545014278341407640D-04 /)
  real ( kind = 8 ), parameter :: const1 = 0.180625D+00
  real ( kind = 8 ), parameter :: const2 = 1.6D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    1.0D+00, &
    2.05319162663775882187D+00, &
    1.67638483018380384940D+00, &
    6.89767334985100004550D-01, &
    1.48103976427480074590D-01, &
    1.51986665636164571966D-02, &
    5.47593808499534494600D-04, &
    1.05075007164441684324D-09 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: e = (/ &
    6.65790464350110377720D+00, &
    5.46378491116411436990D+00, &
    1.78482653991729133580D+00, &
    2.96560571828504891230D-01, &
    2.65321895265761230930D-02, &
    1.24266094738807843860D-03, &
    2.71155556874348757815D-05, &
    2.01033439929228813265D-07 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: f = (/ &
    1.0D+00, &
    5.99832206555887937690D-01, &
    1.36929880922735805310D-01, &
    1.48753612908506148525D-02, &
    7.86869131145613259100D-04, &
    1.84631831751005468180D-05, &
    1.42151175831644588870D-07, &
    2.04426310338993978564D-15 /)
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
!  real ( kind = 8 ) r8poly_value_horner
  real ( kind = 8 ), parameter :: split1 = 0.425D+00
  real ( kind = 8 ), parameter :: split2 = 5.0D+00
  real ( kind = 8 ) x

  if ( p <= 0.0D+00 ) then
    x = - huge ( x )
    return
  end if

  if ( 1.0D+00 <= p ) then
    x = huge ( x )
    return
  end if

  q = p - 0.5D+00

  if ( abs ( q ) <= split1 ) then

    r = const1 - q * q
    x = q * r8poly_value_horner ( 7, a, r ) &
          / r8poly_value_horner ( 7, b, r )

  else

    if ( q < 0.0D+00 ) then
      r = p
    else
      r = 1.0D+00 - p
    end if

    if ( r <= 0.0D+00 ) then

      x = huge ( x )

    else

      r = sqrt ( - log ( r ) )

      if ( r <= split2 ) then

        r = r - const2
        x = r8poly_value_horner ( 7, c, r ) &
          / r8poly_value_horner ( 7, d, r )

      else

        r = r - split2
        x = r8poly_value_horner ( 7, e, r ) &
          / r8poly_value_horner ( 7, f, r )

      end if

    end if

    if ( q < 0.0D+00 ) then
      x = -x
    end if

  end if

  return
end subroutine
subroutine truncated_normal_ab_sample ( mu, sigma, a, b, seed, x )

!*****************************************************************************80
!
!! TRUNCATED_NORMAL_AB_SAMPLE samples the truncated Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) MU, SIGMA, the mean and standard deviation of the
!    parent Normal distribution.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper truncation limits.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) alpha_cdf
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) beta_cdf
  real ( kind = 8 ) mu
!  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) sigma
  integer ( kind = 4 ) seed
  real ( kind = 8 ) u
  real ( kind = 8 ) x
  real ( kind = 8 ) xi
  real ( kind = 8 ) xi_cdf

  alpha = ( a - mu ) / sigma
  beta = ( b - mu ) / sigma

  call normal_01_cdf ( alpha, alpha_cdf )
  call normal_01_cdf ( beta, beta_cdf )

  u = r8_uniform_01 ( seed )
  xi_cdf = alpha_cdf + u * ( beta_cdf - alpha_cdf )
  call normal_01_cdf_inv ( xi_cdf, xi )

  x = mu + sigma * xi

  return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!-------------------------------------------------------------------
! this is a wrapper for the solver-specific init_domain routines.
! Also dereferences the grid variables and passes them down as arguments.
! This is crucial, since the lower level routines may do message passing
! and this will get fouled up on machines that insist on passing down
! copies of assumed-shape arrays (by passing down as arguments, the
! data are treated as assumed-size -- ie. f77 -- arrays and the copying
! business is avoided).  Fie on the F90 designers.  Fie and a pox.
! NOTE:  Modified to remove all but arrays of rank 4 or more from the
!        argument list.  Arrays with rank>3 are still problematic due to the
!        above-noted fie- and pox-ities.  TBH 20061129.

   SUBROUTINE init_domain ( grid )

   IMPLICIT NONE

   !  Input data.
   TYPE (domain), POINTER :: grid
   !  Local data.
   INTEGER :: idum1, idum2


   CALL set_scalar_indices_from_config ( head_grid%id , idum1, idum2 )

     CALL init_domain_rk( grid &
!
#include "actual_new_args.inc"
!
                        )

   END SUBROUTINE init_domain

!-------------------------------------------------------------------

   SUBROUTINE init_domain_rk ( grid &
!
# include "dummy_new_args.inc"
!
)
   IMPLICIT NONE

   !  Input data.
   TYPE (domain), POINTER :: grid

# include "dummy_new_decl.inc"

   TYPE (grid_config_rec_type)              :: config_flags

   !  Local data
   INTEGER                             ::                       &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte, &
                                  i, j, k, kk

   ! Local data

   INTEGER, PARAMETER :: nl_max = 1000
   REAL, DIMENSION(nl_max) :: zk, p_in, theta, rho, u, v, qv, pd_in
   INTEGER :: nl_in, icount


   INTEGER :: icm,jcm, ii, im1, jj, jm1, loop, error, fid, nxc, nyc, lm
   REAL    :: u_mean,v_mean, f0, p_surf, p_level, qvf, z_at_v, z_at_u
   REAL    :: z_scale, xrad, yrad, zrad, rad, delt, cof1, cof2
!   REAL, EXTERNAL :: interp_0
   REAL    :: hm, xa
   REAL (kind=8) :: pi, rnd, rnd2
!  for random normal number generator
   INTEGER (kind=4) :: SEED
   REAL  (kind=8) :: mu, sigma, a, b

!  stuff from original initialization that has been dropped from the Registry
   REAL    :: vnu, xnu, xnus, dinit0, cbh, p0_temp, t0_temp, zd, zt
   REAL    :: qvf1, qvf2, pd_surf, theta_surf
   INTEGER :: it
   real :: thtmp, ptmp, temp(3)
   real :: t_min, t_max, xpos, xposml, xpospl ! used by grav2d_x

   LOGICAL :: moisture_init
   LOGICAL :: stretch_grid, dry_sounding
   character (len=256) :: mminlu2

   REAL    :: xa1, xal1,pii,hm1  !  data for intercomparison setup from dale
!  space for initial jet in b_wave

   INTEGER, parameter :: nz_jet=64, ny_jet=80
   REAL, DIMENSION(nz_jet, ny_jet) :: u_jet, rho_jet, th_jet, z_jet

!  perturbation parameters for b_wave

   REAL, PARAMETER :: htbub=8000., radbub=2000000., radz=8000., tpbub=1.0
   REAL :: piov2, tp
   INTEGER :: icen, jcen

   REAL    :: B1, B2, B3, B4, B5, sin_arg

   REAL    :: Nsq, z, z1, z2
   INTEGER :: iter_loop
   INTEGER :: xs , xe , ys , ye
   REAL :: mtn_ht
   REAL :: randx       ! for LES
   INTEGER :: ks, ke, id

   LOGICAL, EXTERNAL :: wrf_dm_on_monitor

   SELECT CASE ( model_data_order )
         CASE ( DATA_ORDER_ZXY )
   kds = grid%sd31 ; kde = grid%ed31 ;
   ids = grid%sd32 ; ide = grid%ed32 ;
   jds = grid%sd33 ; jde = grid%ed33 ;

   kms = grid%sm31 ; kme = grid%em31 ;
   ims = grid%sm32 ; ime = grid%em32 ;
   jms = grid%sm33 ; jme = grid%em33 ;

   kts = grid%sp31 ; kte = grid%ep31 ;   ! note that tile is entire patch
   its = grid%sp32 ; ite = grid%ep32 ;   ! note that tile is entire patch
   jts = grid%sp33 ; jte = grid%ep33 ;   ! note that tile is entire patch
         CASE ( DATA_ORDER_XYZ )
   ids = grid%sd31 ; ide = grid%ed31 ;
   jds = grid%sd32 ; jde = grid%ed32 ;
   kds = grid%sd33 ; kde = grid%ed33 ;

   ims = grid%sm31 ; ime = grid%em31 ;
   jms = grid%sm32 ; jme = grid%em32 ;
   kms = grid%sm33 ; kme = grid%em33 ;

   its = grid%sp31 ; ite = grid%ep31 ;   ! note that tile is entire patch
   jts = grid%sp32 ; jte = grid%ep32 ;   ! note that tile is entire patch
   kts = grid%sp33 ; kte = grid%ep33 ;   ! note that tile is entire patch
         CASE ( DATA_ORDER_XZY )
   ids = grid%sd31 ; ide = grid%ed31 ;
   kds = grid%sd32 ; kde = grid%ed32 ;
   jds = grid%sd33 ; jde = grid%ed33 ;

   ims = grid%sm31 ; ime = grid%em31 ;
   kms = grid%sm32 ; kme = grid%em32 ;
   jms = grid%sm33 ; jme = grid%em33 ;

   its = grid%sp31 ; ite = grid%ep31 ;   ! note that tile is entire patch
   kts = grid%sp32 ; kte = grid%ep32 ;   ! note that tile is entire patch
   jts = grid%sp33 ; jte = grid%ep33 ;   ! note that tile is entire patch

   END SELECT

   CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags )

  ideal_constants: SELECT CASE ( model_config_rec%ideal_case )
  CASE ( hill2d_x )
   hm = 100.
   xa = 5.0

   icm = ide/2


   xa1  = 5000./500.
   xal1 = 4000./500.
   pii  = 2.*asin(1.0)
   hm1  = 250.
!   hm1  = 1000.


   stretch_grid = .true.
   delt = 0.
!   z_scale = .40
   z_scale = 8000./config_flags%ztop
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/2
   nyc = (jde-jds)/2

  CASE ( quarter_ss, squall2d_x, squall2d_y )
  ! domain initialization

   stretch_grid = .true. 
   delt = 0.02   ! perturbation temperature
!   z_scale = .40
   z_scale = 8000./config_flags%ztop
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/2
   nyc = (jde-jds)/2

  CASE (grav2d_x)

   stretch_grid = .true.
!    z_scale = 1.675
   z_scale = 8000./config_flags%ztop
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/2
   nyc = (jde-jds)/2

  CASE (convrad)
   stretch_grid = .true.
   delt = 1.
!   z_scale = .32
   z_scale = 8000./config_flags%ztop
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/2
   nyc = jde/2
   icm = ide/2
! lm is the half width of the land in terms of grid points
   lm = 25
   write(6,*) 'lm,icm-lm,icm+lm = ', lm,icm-lm,icm+lm

  CASE (b_wave)

! Note: htbub, radbub, radz, tpbub set above as parameters

   piov2 = 2.*atan(1.0)
   icen = ide/4
   jcen = jde/2

   stretch_grid = .true.
   delt = 0.
!   z_scale = .50
   z_scale = 8000./config_flags%ztop
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/4
   nyc = (jde-jds)/2

  CASE (seabreeze2d_x)

   stretch_grid = .true.
   delt = 6.
!   z_scale = .40
   z_scale = 8000./config_flags%ztop
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/2
   nyc = jde/2
   icm = ide/2
! lm is the half width of the land in terms of grid points
   lm = 25
   write(6,*) 'lm,icm-lm,icm+lm = ', lm,icm-lm,icm+lm

  CASE (les)
  !  FOR LES, set stretch to false
   stretch_grid = .false.
   delt = 3.
!   z_scale = .40
   z_scale = 8000./config_flags%ztop  ! not used if stretch_grid = .false.
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/2
   nyc = (jde-jds)/2


  CASE DEFAULT

      WRITE( wrf_err_message , * ) 'Need to choose valid non-zero ideal_case:  ideal_case = ', model_config_rec%ideal_case
      CALL wrf_error_fatal ( wrf_err_message )

  END SELECT ideal_constants




! here we check to see if the boundary conditions are set properly

   CALL boundary_condition_check( config_flags, bdyzone, error, grid%id )

   moisture_init = .true.

    grid%itimestep=0

#ifdef DM_PARALLEL
   CALL wrf_dm_bcast_bytes( icm , IWORDSIZE )
   CALL wrf_dm_bcast_bytes( jcm , IWORDSIZE )
#endif

  ideal_landmap: SELECT CASE ( model_config_rec%ideal_case )

   CASE(hill2d_x, squall2d_x, squall2d_y, grav2d_x)

    CALL nl_set_mminlu(1,'    ')
    CALL nl_set_iswater(1,0)
    CALL nl_set_cen_lat(1,40.)
    CALL nl_set_cen_lon(1,-105.)
    CALL nl_set_truelat1(1,0.)
    CALL nl_set_truelat2(1,0.)
    CALL nl_set_moad_cen_lat (1,0.)
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_pole_lon (1,0.)
    CALL nl_set_pole_lat (1,90.)
    CALL nl_set_map_proj(1,0)


!  here we initialize data that currently is not initialized
!  in the input data

    DO j = jts, jte
      DO i = its, ite
         grid%msftx(i,j)    = 1.
         grid%msfty(i,j)    = 1.
         grid%msfux(i,j)    = 1.
         grid%msfuy(i,j)    = 1.
         grid%msfvx(i,j)    = 1.
         grid%msfvx_inv(i,j)= 1.
         grid%msfvy(i,j)    = 1.
         grid%sina(i,j)     = 0.
         grid%cosa(i,j)     = 1.
         grid%e(i,j)        = 0.
         grid%f(i,j)        = 0.

      END DO
   END DO

   CASE(b_wave)

    CALL nl_set_mminlu(1,'    ')
    CALL nl_set_iswater(1,0)
    CALL nl_set_cen_lat(1,40.)
    CALL nl_set_cen_lon(1,-105.)
    CALL nl_set_truelat1(1,0.)
    CALL nl_set_truelat2(1,0.)
    CALL nl_set_moad_cen_lat (1,0.)
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_pole_lon (1,0.)
    CALL nl_set_pole_lat (1,90.)
    CALL nl_set_map_proj(1,0)


!  here we initialize data that currently is not initialized
!  in the input data

    DO j = jts, jte
      DO i = its, ite
         grid%msftx(i,j)    = 1.
         grid%msfty(i,j)    = 1.
         grid%msfux(i,j)    = 1.
         grid%msfuy(i,j)    = 1.
         grid%msfvx(i,j)    = 1.
         grid%msfvx_inv(i,j)= 1.
         grid%msfvy(i,j)    = 1.
         grid%sina(i,j)     = 0.
         grid%cosa(i,j)     = 1.
         grid%e(i,j)        = 0.
         grid%f(i,j)        = 1.e-04

      END DO
   END DO

   CASE(convrad, quarter_ss) ! @lkugler: add radiation in quarter_ss
    mminlu2 = ' '
    mminlu2(1:4) = 'USGS'
    CALL nl_set_mminlu(1, mminlu2)
!   CALL nl_set_mminlu(1, 'USGS')
    CALL nl_set_iswater(1,16)
    CALL nl_set_isice(1,3)
    CALL nl_set_cen_lat(1,45.)
    CALL nl_set_cen_lon(1,0.)
    CALL nl_set_truelat1(1,45.)
    CALL nl_set_truelat2(1,45.)
    CALL nl_set_moad_cen_lat (1,45.)
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_pole_lon (1,0.)
    CALL nl_set_pole_lat (1,90.)
    CALL nl_set_map_proj(1,0)
!   CALL model_to_grid_config_rec(1,model_config_rec,config_flags)
    CALL nl_get_iswater(1,grid%iswater)

!  here we initialize data that currently is not initialized
!  in the input data

    DO j = jts, jte
      DO i = its, ite
         grid%msft(i,j)     = 1.
         grid%msfu(i,j)     = 1.
         grid%msfv(i,j)     = 1.
         grid%msftx(i,j)    = 1.
         grid%msfty(i,j)    = 1.
         grid%msfux(i,j)    = 1.
         grid%msfuy(i,j)    = 1.
         grid%msfvx(i,j)    = 1.
         grid%msfvy(i,j)    = 1.
         grid%msfvx_inv(i,j)= 1.
         grid%sina(i,j)     = 0.
         grid%cosa(i,j)     = 1.
         grid%e(i,j)        = 0.
         grid%xlat(i,j)     = 45.
         grid%f(i,j)        = 0.  ! @lkugler
         grid%xlong(i,j)     = 0.
! Hard-wire the ocean-land configuration
! all ocean
!       if (i .ge. (icm-lm) .and. i .lt. (icm+lm)) then
!        grid%xland(i,j)     = 1.
!        grid%lu_index(i,j)  = 18
!else
         ! @lkugler: use LAND, 'Dryland Cropland and Pasture'
         grid%xland(i,j)     = 1.
         grid%lu_index(i,j)  = 2
!end if
      END DO
   END DO

! for Noah LSM, additional variables need to be initialized

   convrad_masked_fields : SELECT CASE ( model_config_rec%sf_surface_physics(grid%id) )

      CASE (SLABSCHEME)

      CASE (LSMSCHEME)

        DO j = jts , MIN(jde-1,jte)
           DO i = its , MIN(ide-1,ite)
              IF (grid%xland(i,j) .lt. 1.5) THEN
                 grid%vegfra(i,j) = 50.
                 grid%canwat(i,j) = 0.
                 grid%ivgtyp(i,j) = 2  ! @lkugler
                 grid%isltyp(i,j) = 6  ! @lkugler LOAM
                 grid%xice(i,j) = 0.
                 grid%snow(i,j) = 0.
              ELSE
                 grid%vegfra(i,j) = 0.
                 grid%canwat(i,j) = 0.
                 grid%ivgtyp(i,j) = 16
                 grid%isltyp(i,j) = 14
                 grid%xice(i,j) = 0.
                 grid%snow(i,j) = 0.
              ENDIF
           END DO
        END DO

      CASE (RUCLSMSCHEME)

   END SELECT convrad_masked_fields
! Process the soil; note that there are some things hard-wired into share/module_soil_pre.F
      CALL process_soil_ideal(grid%xland,grid%xice,grid%vegfra,grid%snow,grid%canwat, &
                     grid%ivgtyp,grid%isltyp,grid%tslb,grid%smois, &
                     grid%tsk,grid%tmn,grid%zs,grid%dzs,model_config_rec%num_soil_layers, &
                     model_config_rec%sf_surface_physics(grid%id), &
                                   ids,ide, jds,jde, kds,kde,&
                                   ims,ime, jms,jme, kms,kme,&
                                   its,ite, jts,jte, kts,kte )
   CASE(seabreeze2d_x)

    mminlu2 = ' '
    mminlu2(1:4) = 'USGS'
    CALL nl_set_mminlu(1, mminlu2)
!   CALL nl_set_mminlu(1, 'USGS')
    CALL nl_set_iswater(1,16)
    CALL nl_set_isice(1,3)
    CALL nl_set_cen_lat(1,20.)
    CALL nl_set_cen_lon(1,-105.)
    CALL nl_set_truelat1(1,0.)
    CALL nl_set_truelat2(1,0.)
    CALL nl_set_moad_cen_lat (1,0.)
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_pole_lon (1,0.)
    CALL nl_set_pole_lat (1,90.)
    CALL nl_set_map_proj(1,0)
!   CALL model_to_grid_config_rec(1,model_config_rec,config_flags)
    CALL nl_get_iswater(1,grid%iswater)

!  here we initialize data that currently is not initialized
!  in the input data

    DO j = jts, jte
      DO i = its, ite
         grid%msft(i,j)     = 1.
         grid%msfu(i,j)     = 1.
         grid%msfv(i,j)     = 1.
         grid%msftx(i,j)    = 1.
         grid%msfty(i,j)    = 1.
         grid%msfux(i,j)    = 1.
         grid%msfuy(i,j)    = 1.
         grid%msfvx(i,j)    = 1.
         grid%msfvy(i,j)    = 1.
         grid%msfvx_inv(i,j)= 1.
         grid%sina(i,j)     = 0.
         grid%cosa(i,j)     = 1.
         grid%e(i,j)        = 0.
         grid%f(i,j)        = 0.
         grid%xlat(i,j)     = 30.
         grid%xlong(i,j)     = 0.
! Hard-wire the ocean-land configuration
        if (i .ge. (icm-lm) .and. i .lt. (icm+lm)) then
         grid%xland(i,j)     = 1.
         grid%lu_index(i,j)  = 18
         grid%tsk(i,j) = 280.0
         grid%tmn(i,j) = 280.0
        else
         grid%xland(i,j)     = 2.
         grid%lu_index(i,j)  = 16
         grid%tsk(i,j) = 287.0
         grid%tmn(i,j) = 280.0
        end if
      END DO
   END DO

! for Noah LSM, additional variables need to be initialized

   seabreeze_masked_fields : SELECT CASE ( model_config_rec%sf_surface_physics(grid%id) )

      CASE (SLABSCHEME)

      CASE (LSMSCHEME)

        DO j = jts , MIN(jde-1,jte)
           DO i = its , MIN(ide-1,ite)
              IF (grid%xland(i,j) .lt. 1.5) THEN
                 grid%vegfra(i,j) = 50.
                 grid%canwat(i,j) = 0.
                 grid%ivgtyp(i,j) = 18
                 grid%isltyp(i,j) = 8
                 grid%xice(i,j) = 0.
                 grid%snow(i,j) = 0.
              ELSE
                 grid%vegfra(i,j) = 0.
                 grid%canwat(i,j) = 0.
                 grid%ivgtyp(i,j) = 16
                 grid%isltyp(i,j) = 14
                 grid%xice(i,j) = 0.
                 grid%snow(i,j) = 0.
              ENDIF
           END DO
        END DO

      CASE (RUCLSMSCHEME)

   END SELECT seabreeze_masked_fields

! Process the soil; note that there are some things hard-wired into share/module_soil_pre.F
      CALL process_soil_ideal(grid%xland,grid%xice,grid%vegfra,grid%snow,grid%canwat, &
                     grid%ivgtyp,grid%isltyp,grid%tslb,grid%smois, &
                     grid%tsk,grid%tmn,grid%zs,grid%dzs,model_config_rec%num_soil_layers, &
                     model_config_rec%sf_surface_physics(grid%id), &
                                   ids,ide, jds,jde, kds,kde,&
                                   ims,ime, jms,jme, kms,kme,&
                                   its,ite, jts,jte, kts,kte )

   CASE (les)

    CALL nl_set_mminlu(1, '    ')
    CALL nl_set_iswater(1,0)
    CALL nl_set_cen_lat(1,40.)
    CALL nl_set_cen_lon(1,-105.)
    CALL nl_set_truelat1(1,0.)
    CALL nl_set_truelat2(1,0.)
    CALL nl_set_moad_cen_lat (1,0.)
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_pole_lon (1,0.)
    CALL nl_set_pole_lat (1,90.)
    CALL nl_set_map_proj(1,0)


!  here we initialize data that currently is not initialized
!  in the input data

    DO j = jts, jte
      DO i = its, ite
         grid%msftx(i,j)    = 1.
         grid%msfty(i,j)    = 1.
         grid%msfux(i,j)    = 1.
         grid%msfuy(i,j)    = 1.
         grid%msfvx(i,j)    = 1.
         grid%msfvx_inv(i,j)= 1.
         grid%msfvy(i,j)    = 1.
         grid%sina(i,j)     = 0.
         grid%cosa(i,j)     = 1.
         grid%e(i,j)        = 0.
!  for LES, include Coriolis force
         grid%f(i,j)        = 1.e-4

      END DO
   END DO

   END SELECT ideal_landmap

    DO j = jts, jte
    DO k = kts, kte
      DO i = its, ite
         grid%ww(i,k,j)     = 0.
      END DO
   END DO
   END DO

   grid%step_number = 0

! set up the grid

   ideal_levels: SELECT CASE ( model_config_rec%ideal_case )
   CASE(hill2d_x, squall2d_x, squall2d_y, grav2d_x, b_wave)
   IF (stretch_grid) THEN ! exponential stretch for eta (nearly constant dz)
     DO k=1, kde
      grid%znw(k) = (exp(-(k-1)/float(kde-1)/z_scale) - exp(-1./z_scale))/ &
                                (1.-exp(-1./z_scale))
     ENDDO
   ELSE
     DO k=1, kde
      grid%znw(k) = 1. - float(k-1)/float(kde-1)
     ENDDO
   ENDIF

   CASE(convrad, seabreeze2d_x, quarter_ss)  ! lkugler add quarter
   IF (stretch_grid) THEN ! read eta_levels
     DO k=1, kde
      grid%znw(k) = model_config_rec%eta_levels(k)
     ENDDO

 ! Check the value of the first and last eta level for our domain,
 ! then check that the vector of eta levels is only decreasing

     IF (model_config_rec%eta_levels(1) .NE. 1.0) THEN
        CALL wrf_error_fatal("--- ERROR: the first specified eta_level is not 1.0")
     ENDIF
     IF (model_config_rec%eta_levels(kde) .NE. 0.0) THEN
        CALL wrf_error_fatal("--- ERROR: the last specified eta_level is not 0.0")
     ENDIF
     DO k=2,kde
       IF (model_config_rec%eta_levels(k) .GT. model_config_rec%eta_levels(k-1)) THEN
          CALL wrf_error_fatal("--- ERROR: specified eta_levels are not uniformly decreasing from 1.0 to 0.0")
       ENDIF
     ENDDO
   ELSE
     DO k=1, kde
      grid%znw(k) = 1. - float(k-1)/float(kde-1)
     ENDDO
   ENDIF
   CASE(les)
       IF (model_config_rec%eta_levels(1) .EQ. -1) THEN !we do not have eta_levels from namelist
   IF (stretch_grid) THEN ! exponential stretch for eta (nearly constant dz)
     DO k=1, kde
      grid%znw(k) = (exp(-(k-1)/float(kde-1)/z_scale) - exp(-1./z_scale))/ &
                                (1.-exp(-1./z_scale))
     ENDDO
   ELSE
     DO k=1, kde
      grid%znw(k) = 1. - float(k-1)/float(kde-1)
     ENDDO
   ENDIF
      ELSE
          CALL wrf_debug(0,"module_initialize_les: vertical nesting is enabled, using eta_levels specified in namelist.input")
          ks = 0
          DO id=1,grid%id
             ks = ks+model_config_rec%e_vert(id)
          ENDDO
          IF (ks .GT. max_eta) THEN
             CALL wrf_error_fatal("too many vertical levels, increase max_eta in frame/module_driver_constants.F")
          ENDIF
          !Now set the eta_levels to what we specified in the namelist. We've
          !packed all the domains' eta_levels into a 'vector' and now we need
          !to pull only the section of the vector associated with our domain
          !of interest, which is between indicies ks and ke.
          IF (grid%id .EQ. 1) THEN
            ks = 1
            ke = model_config_rec%e_vert(1)
          ELSE
            id = 1
            ks = 1
            ke = 0
            DO WHILE (grid%id .GT. id)
              id = id+1
              ks = ks+model_config_rec%e_vert(id-1)
              ke = ks+model_config_rec%e_vert(id)
            ENDDO
          ENDIF
          DO k=1,kde
            grid%znw(k) = model_config_rec%eta_levels(ks+k-1)
          ENDDO
          !Check the value of the first and last eta level for our domain,
          !then check that the vector of eta levels is only decreasing
          IF (grid%znw(1) .NE. 1.0) THEN
            CALL wrf_error_fatal("error with specified eta_levels, first level is not 1.0")
          ENDIF
          IF (grid%znw(kde) .NE. 0.0) THEN
            CALL wrf_error_fatal("error with specified eta_levels, last level is not 0.0")
          ENDIF
          DO k=2,kde
            IF (grid%znw(k) .GT. grid%znw(k-1)) THEN
              CALL wrf_error_fatal("eta_levels are not uniformly decreasing from 1.0 to 0.0")
            ENDIF
          ENDDO
      ENDIF
   END SELECT ideal_levels

   DO k=1, kde-1
    grid%dnw(k) = grid%znw(k+1) - grid%znw(k)
    grid%rdnw(k) = 1./grid%dnw(k)
    grid%znu(k) = 0.5*(grid%znw(k+1)+grid%znw(k))
   ENDDO
   DO k=2, kde-1
    grid%dn(k) = 0.5*(grid%dnw(k)+grid%dnw(k-1))
    grid%rdn(k) = 1./grid%dn(k)
    grid%fnp(k) = .5* grid%dnw(k  )/grid%dn(k)
    grid%fnm(k) = .5* grid%dnw(k-1)/grid%dn(k)
   ENDDO

   cof1 = (2.*grid%dn(2)+grid%dn(3))/(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(2)
   cof2 =     grid%dn(2)        /(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(3)
   grid%cf1  = grid%fnp(2) + cof1
   grid%cf2  = grid%fnm(2) - cof1 - cof2
   grid%cf3  = cof2

   grid%cfn  = (.5*grid%dnw(kde-1)+grid%dn(kde-1))/grid%dn(kde-1)
   grid%cfn1 = -.5*grid%dnw(kde-1)/grid%dn(kde-1)
   grid%rdx = 1./config_flags%dx
   grid%rdy = 1./config_flags%dy

!  get the sounding from the ascii sounding file, first get dry sounding and
!  calculate base state
   ideal_sounding: SELECT CASE ( model_config_rec%ideal_case )

   CASE (b_wave)

  write(6,*) ' reading input jet sounding '
  call read_input_jet( u_jet, rho_jet, th_jet, z_jet, nz_jet, ny_jet )

  write(6,*) ' getting dry sounding for base state '
  write(6,*) ' using middle column in jet sounding, j = ',ny_jet/2

  dry_sounding   = .true.
!  debug = .true.  !  this will produce print of the sounding (last argument)
  CALL get_sounding_b_wave( zk, p_in, pd_in, theta, rho, u, v, qv, dry_sounding, &
                      nl_max, nl_in, u_jet, rho_jet, th_jet, z_jet,      &
                      nz_jet, ny_jet, ny_jet/2, .true.                   )

  write(6,*) ' returned from reading sounding, nl_in is ',nl_in

!  find ptop for the desired ztop (ztop is input from the namelist),
!  and find surface pressure

!  For the jet, using the middle column for the base state means that
!  we will be extrapolating above the highest height data to the south
!  of the centerline.


   CASE DEFAULT

  IF ( wrf_dm_on_monitor() ) THEN
  write(6,*) ' getting dry sounding for base state '
  dry_sounding = .true.
  CALL get_sounding( zk, p_in, pd_in, theta, rho, u, v, qv, dry_sounding, &
                     nl_max, nl_in, theta_surf)
  ENDIF
  CALL wrf_dm_bcast_real( zk , nl_max )
  CALL wrf_dm_bcast_real( p_in , nl_max )
  CALL wrf_dm_bcast_real( pd_in , nl_max )
  CALL wrf_dm_bcast_real( theta , nl_max )
  CALL wrf_dm_bcast_real( rho , nl_max )
  CALL wrf_dm_bcast_real( u , nl_max )
  CALL wrf_dm_bcast_real( v , nl_max )
  CALL wrf_dm_bcast_real( qv , nl_max )
  CALL wrf_dm_bcast_integer ( nl_in , 1 )

  write(6,*) ' returned from reading sounding, nl_in is ',nl_in


!  find ptop for the desired ztop (ztop is input from the namelist),
!  and find surface pressure
   END SELECT ideal_sounding

  grid%p_top = interp_0( p_in, zk, config_flags%ztop, nl_in )

!  Fill in the hybrid coordinate coefficients

   DO k=1, kde
      IF      ( config_flags%hybrid_opt .EQ. 0 ) THEN
         grid%c3f(k) = grid%znw(k)
      ELSE IF ( config_flags%hybrid_opt .EQ. 1 ) THEN
         grid%c3f(k) = grid%znw(k)
      ELSE IF ( config_flags%hybrid_opt .EQ. 2 ) THEN
         B1 = 2. * grid%etac**2 * ( 1. - grid%etac )
         B2 = -grid%etac * ( 4. - 3. * grid%etac - grid%etac**3 )
         B3 = 2. * ( 1. - grid%etac**3 )
         B4 = - ( 1. - grid%etac**2 )
         B5 = (1.-grid%etac)**4
         grid%c3f(k) = ( B1 + B2*grid%znw(k) + B3*grid%znw(k)**2 + B4*grid%znw(k)**3 ) / B5
         IF ( grid%znw(k) .LT. grid%etac ) THEN
            grid%c3f(k) = 0.
         END IF
         IF ( k .EQ. kds ) THEN
            grid%c3f(k) = 1.
         ELSE IF ( k .EQ. kde ) THEN
            grid%c3f(k) = 0.
         END IF
      ELSE IF ( config_flags%hybrid_opt .EQ. 3 ) THEN
         IF ( grid%znw(k) .GE. grid%etac ) THEN
            sin_arg = (1./(1.-grid%etac))*(grid%znw(k)-1.)+1
            grid%c3f(k) = (sin(sin_arg*3.14159265358/2.))**2
         ELSE
            grid%c3f(k) = 0.
         END IF
         IF ( k .EQ. kds ) THEN
            grid%c3f(k) = 1.
         ELSE IF ( k .EQ. kds ) THEN
            grid%c3f(kde) = 0.
         END IF
      ELSE
         CALL wrf_error_fatal ( 'ERROR: --- hybrid_opt=0 ===> Standard WRF Coordinate; hybrid_opt>=1 ===> Hybrid Vertical Coordinate' )
      END IF
   END DO

   DO k=1, kde
      grid%c4f(k) = ( grid%znw(k) - grid%c3f(k) ) * ( p1000mb - grid%p_top )
   ENDDO

   !  Now on half levels, just add up and divide by 2 (for c3h).  Use (eta-c3)*(p00-pt) for c4 on half levels.

   DO k=1, kde-1
      grid%c3h(k) = ( grid%c3f(k+1) + grid%c3f(k) ) * 0.5
      grid%c4h(k) = ( grid%znu(k) - grid%c3h(k) ) * ( p1000mb - grid%p_top )
   ENDDO

   !  c1 = d(B)/d(eta).  We define c1f as c1 on FULL levels.  For a vertical difference,
   !  we need to use B and eta on half levels.  The k-loop ends up referring to the
   !  full levels, neglecting the top and bottom.

   DO k=kds+1, kde-1
      grid%c1f(k) = ( grid%c3h(k) - grid%c3h(k-1) ) / ( grid%znu(k) - grid%znu(k-1) )
   ENDDO

   !  The boundary conditions to get the coefficients:
   !  1) At k=kts: define d(B)/d(eta) = 1.  This gives us the same value of B and d(B)/d(eta)
   !     when doing the sigma-only B=eta.
   !  2) At k=kte: with the new vertical coordinate, define d(B)/d(eta) = 0.  The curve B SMOOTHLY
   !     goes to zero, and at the very top, B continues to SMOOTHLY go to zero.  Note that for
   !     almost all cases of non B=eta, B is ALREADY=ZERO at the top, so this is a reasonable BC to
   !     assume.
   !  3) At k=kte: when trying to mimic the original vertical coordinate, since B = eta, then
   !     d(B)/d(eta) = 1.

   grid%c1f(kds) = 1.
   IF      ( ( config_flags%hybrid_opt .EQ. 0 ) .OR. ( config_flags%hybrid_opt .EQ. 1 ) ) THEN
      grid%c1f(kde) = 1.
   ELSE
      grid%c1f(kde) = 0.
   END IF

   !  c2 = ( 1. - c1(k) ) * (p00 - pt).  There is no vertical differencing, so we can do the
   !  full kds to kde looping.

   DO k=kds, kde
      grid%c2f(k) = ( 1. - grid%c1f(k) ) * ( p1000mb - grid%p_top )
   END DO

   !  Now on half levels for c1 and c2.  The c1h will result from the full level c3 and full
   !  level eta differences.  The c2 value use the half level c1(k).

   DO k=1, kde-1
      grid%c1h(k) = ( grid%c3f(k+1) - grid%c3f(k) ) / ( grid%znw(k+1) - grid%znw(k) )
      grid%c2h(k) = ( 1. - grid%c1h(k) ) * ( p1000mb - grid%p_top )
   END DO

#if 0
   DO k=1, kde
      grid%c3f(k) = grid%znw(k)
      grid%c4f(k) = 0.
      grid%c3h(k) = grid%znu(k)
      grid%c4h(k) = 0.
      grid%c1f(k) = 1.
      grid%c2f(k) = 0.
      grid%c1h(k) = 1.
      grid%c2h(k) = 0.
   END DO
#endif

  ideal_terrain: SELECT CASE ( model_config_rec%ideal_case )
  CASE (hill2d_x)
  DO j=jts,jte
  DO i=its,ite
    grid%ht(i,j) = hm/(1.+(float(i-icm)/xa)**2)
!!    grid%ht(i,j) = hm1*exp(-(( float(i-icm)/xa1)**2))   &
!!               *( (cos(pii*float(i-icm)/xal1))**2 )
    grid%phb(i,1,j) = g*grid%ht(i,j)
    grid%php(i,1,j) = 0.
    grid%ph0(i,1,j) = grid%phb(i,1,j)
  ENDDO
  ENDDO

  CASE (quarter_ss, convrad, squall2d_x, squall2d_y, grav2d_x, b_wave, seabreeze2d_x)
  DO j=jts,jte
  DO i=its,ite
    grid%ht(i,j) = 489.  ! @lkugler surface altitude
    grid%phb(i,1,j) = g*grid%ht(i,j)
    grid%php(i,1,j) = 0.  ! perturbation geopotential unaffected by z0 != 0 
    grid%ph0(i,1,j) = grid%phb(i,1,j)
  ENDDO
  ENDDO

  CASE (les)
  DO j=jts,jte
  DO i=its,ite
    grid%ht(i,j) = 0.
  ENDDO
  ENDDO

  xs=ide/2 -3
  xs=ids   -3
  xe=xs + 6
  ys=jde/2 -3
  ye=ys + 6
  mtn_ht = 500
#ifdef MTN
  DO j=max(ys,jds),min(ye,jde-1)
  DO i=max(xs,ids),min(xe,ide-1)
     grid%ht(i,j) = mtn_ht * 0.25 * &
               ( 1. + COS ( 2*pi/(xe-xs) * ( i-xs ) + pi ) ) * &
               ( 1. + COS ( 2*pi/(ye-ys) * ( j-ys ) + pi ) )
  ENDDO
  ENDDO
#endif
#ifdef EW_RIDGE
  DO j=max(ys,jds),min(ye,jde-1)
  DO i=ids,ide
     grid%ht(i,j) = mtn_ht * 0.50 * &
               ( 1. + COS ( 2*pi/(ye-ys) * ( j-ys ) + pi ) )
  ENDDO
  ENDDO
#endif
#ifdef NS_RIDGE
  DO j=jds,jde
  DO i=max(xs,ids),min(xe,ide-1)
     grid%ht(i,j) = mtn_ht * 0.50 * &
               ( 1. + COS ( 2*pi/(xe-xs) * ( i-xs ) + pi ) )
  ENDDO
  ENDDO
#endif
  DO j=jts,jte
  DO i=its,ite
    grid%phb(i,1,j) = g * grid%ht(i,j)
    grid%ph0(i,1,j) = g * grid%ht(i,j)
  ENDDO
  ENDDO

  END SELECT ideal_terrain

  DO J = jts, jte
  DO I = its, ite
    p_surf = interp_0( p_in, zk, grid%phb(i,1,j)/g, nl_in )
    grid%MUB(i,j) = p_surf-grid%p_top

   !  this is dry hydrostatic sounding (base state), so given grid%p (coordinate),
   !  interp theta (from interp) and compute 1/rho from eqn. of state
   !write(*,*) 'p_surf, grid%p_top', p_surf, grid%p_top
    DO K = 1, kte-1
      p_level = grid%c3h(k)*(p_surf - grid%p_top) + grid%c4h(k) + grid%p_top
      grid%pb(i,k,j) = p_level
      grid%t_init(i,k,j) = interp_0( theta, p_in, p_level, nl_in ) - t0
      grid%alb(i,k,j) = (r_d/p1000mb)*(grid%t_init(i,k,j)+t0)*(grid%pb(i,k,j)/p1000mb)**cvpm
    ENDDO
!  calc hydrostatic balance (alternatively we could interp the geopotential from the
!  sounding, but this assures that the base state is in exact hydrostatic balance with
!  respect to the model eqns.
    DO kk  = 2,kte
      k=kk - 1
      grid%phb(i,kk,j) = grid%phb(i,kk-1,j) - grid%dnw(kk-1)*(grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))*grid%alb(i,kk-1,j)
    ENDDO
  ENDDO
  ENDDO
  IF ( wrf_dm_on_monitor() ) THEN
  write(6,*) ' ptop is ',grid%p_top
  write(6,*) ' base state grid%MUB(1,1), p_surf is ',grid%MUB(1,1),grid%c3f(kts)*grid%MUB(1,1)+grid%c4f(kts)+grid%p_top
  ENDIF
!  calculate full state for each column - this includes moisture.
  write(6,*) ' getting moist sounding for full state '
  IF ( model_config_rec%ideal_case .EQ. b_wave )THEN
    dry_sounding = .true.
    IF (config_flags%mp_physics /= 0)  dry_sounding = .false.
  ELSE
    dry_sounding = .false.
    CALL get_sounding( zk, p_in, pd_in, theta, rho, u, v, qv, dry_sounding, &
                     nl_max, nl_in, theta_surf )
  ENDIF
  DO J = jts, min(jde-1,jte)

  IF ( model_config_rec%ideal_case .EQ. b_wave )THEN
!  get sounding for this point
!   debug = .false.  !  this will turn off print of the sounding
    CALL get_sounding_b_wave( zk, p_in, pd_in, theta, rho, u, v, qv, dry_sounding, &
                      nl_max, nl_in, u_jet, rho_jet, th_jet, z_jet,      &
                      nz_jet, ny_jet, j, .false.                          )
  ENDIF
  DO I = its, min(ide-1,ite)
!  At this point grid%p_top is already set. find the DRY mass in the column
!  by interpolating the DRY pressure.
   pd_surf = interp_0( pd_in, zk, grid%phb(i,1,j)/g, nl_in )
!  compute the perturbation mass and the full mass
    grid%MU_1(i,j) = pd_surf-grid%p_top - grid%MUB(i,j)
    grid%MU_2(i,j) = grid%MU_1(i,j)
    grid%MU0(i,j) = grid%MU_1(i,j) + grid%MUB(i,j)
! given the dry pressure and coordinate system, interp the potential
! temperature and qv
    do k=1,kde-1
      p_level = grid%c3h(k)*(pd_surf - grid%p_top) + grid%c4h(k) + grid%p_top
      moist(i,k,j,P_QV) = interp_0( qv, pd_in, p_level, nl_in )
      grid%t_1(i,k,j)          = interp_0( theta, pd_in, p_level, nl_in ) - t0
      grid%t_2(i,k,j)          = grid%t_1(i,k,j)

    enddo
!  integrate the hydrostatic equation (from the RHS of the bigstep
!  vertical momentum equation) down from the top to get grid%p.
!  first from the top of the model to the top pressure
    kk = kte-1  ! top level
    k=kk+1
    qvf1 = 0.5*(moist(i,kk,j,P_QV)+moist(i,kk,j,P_QV))
    qvf2 = 1./(1.+qvf1)
    qvf1 = qvf1*qvf2
    write(*,*) grid%c1f(k), grid%c2f(k), grid%Mu_1(i,j), grid%Mub(i,j), &
               grid%rdnw(kk), qvf1, qvf2
    grid%p(i,kk,j) =  - 0.5*((grid%c1f(k)*grid%Mu_1(i,j))+qvf1*(grid%c1f(k)*grid%Mub(i,j)+grid%c2f(k)))/grid%rdnw(kk)/qvf2
    qvf = 1. + rvovrd*moist(i,kk,j,P_QV)
    grid%alt(i,kk,j) = (r_d/p1000mb)*(grid%t_1(i,kk,j)+t0)*qvf* &
                (((grid%p(i,kk,j)+grid%pb(i,kk,j))/p1000mb)**cvpm)
    grid%al(i,kk,j) = grid%alt(i,kk,j) - grid%alb(i,kk,j)
!  down the column
    do kk=kte-2,1,-1
      k = kk + 1
      qvf1 = 0.5*(moist(i,kk,j,P_QV)+moist(i,kk+1,j,P_QV))  ! lkugler: full level qvapor top of layer
      qvf2 = 1./(1.+qvf1)   ! lkugler: this is the specific humdity (=wrt. moist air) 
      qvf1 = qvf1*qvf2      ! lkugler: why is this necessary?
      grid%p(i,kk,j) =  grid%p(i,kk+1,j) - ((grid%c1f(k)*grid%Mu_1(i,j)) + qvf1*(grid%c1f(k)*grid%Mub(i,j)+grid%c2f(k)))/qvf2/grid%rdn(kk+1)
      qvf = 1. + rvovrd*moist(i,kk,j,P_QV)
      grid%alt(i,kk,j) = (r_d/p1000mb)*(grid%t_1(i,kk,j)+t0)*qvf* &
                  (((grid%p(i,kk,j)+grid%pb(i,kk,j))/p1000mb)**cvpm)
      grid%al(i,kk,j) = grid%alt(i,kk,j) - grid%alb(i,kk,j)
    enddo
!  this is the hydrostatic equation used in the model after the
!  small timesteps.  In the model, grid%al (inverse density)
!  is computed from the geopotential.
    grid%ph_1(i,1,j) = 0.
    DO kk  = 2,kte
      k = kk-1
      grid%ph_1(i,kk,j) = grid%ph_1(i,kk-1,j) - (grid%dnw(kk-1))*(       &
                   ((grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))+(grid%c1h(k)*grid%mu_1(i,j)))*grid%al(i,kk-1,j)+ &
                    (grid%c1h(k)*grid%mu_1(i,j))*grid%alb(i,kk-1,j)  )

      grid%ph_2(i,kk,j) = grid%ph_1(i,kk,j)
      grid%ph0(i,kk,j) = grid%ph_1(i,kk,j) + grid%phb(i,kk,j)
    ENDDO
    IF ( wrf_dm_on_monitor() ) THEN
    if((i==2) .and. (j==2)) then
     k=1
     write(6,*) ' grid%ph_1 k=1 calc ',grid%ph_1(2,k,2),&
                              (grid%c1h(k)*grid%mu_1(2,2))+(grid%c1h(k)*grid%mub(2,2)+grid%c2h(k)),(grid%c1h(k)*grid%mu_1(2,2)), &
                              grid%alb(2,k,2),grid%rdnw(k)
     k=2
     write(6,*) ' grid%ph_1 k=2 calc ',grid%ph_1(2,k,2),&
                              (grid%c1h(k)*grid%mu_1(2,2))+(grid%c1h(k)*grid%mub(2,2)+grid%c2h(k)),(grid%c1h(k)*grid%mu_1(2,2)), &
                              grid%alb(2,k,2)
    endif
    ENDIF
  IF ( model_config_rec%ideal_case .EQ. b_wave )THEN
    DO K = 1, kte
      p_level = grid%c3h(k)*(p_surf - grid%p_top) + grid%c4h(k) + grid%p_top
      grid%u_1(i,k,j) = interp_0( u, p_in, p_level, nl_in )
      grid%u_2(i,k,j) = grid%u_1(i,k,j)
    ENDDO
  ENDIF

  ENDDO
  ENDDO

!===============
  ideal_pert: SELECT CASE ( model_config_rec%ideal_case )
  CASE (quarter_ss)
! -------------------------------------
! @lkugler: add thermal perturbations

!  random low-level thermal perturbation to kick off convection
! call random_seed
write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
seed = 1

  DO J = jts, min(jde-1,jte)
    DO I = its, min(ide-1,ite)

      mu = 0.
      sigma = 0.02
      a = -10.
      b = 10.
      call truncated_normal_ab_sample( mu, sigma, a, b, seed, rnd ) !RANDN(rnd)
      DO K = 1, 20  ! will and should never reach 20 with enough eta_levels
          ! WRITE(*,*) 'PB', grid%pb(i,k,j)
          IF (grid%pb(i,k,j).LT.(grid%pb(i,1,j)-10000.)) EXIT  ! lowest 100hPa
          rnd2 = rnd*EXP((grid%pb(i,k,j)-grid%pb(i,1,j))/2500.)  ! same sign in one column, smoothes out
          if ((J.eq.jts).and.(I.eq.its)) then
              write(*,*) 'Tpert  k: ', k, ', rnd:', rnd2
          end if

          grid%t_1(i,k,j)=grid%t_1(i,k,j)+rnd2  ! add perturbation
          grid%t_2(i,k,j)=grid%t_1(i,k,j)
          qvf = 1. + rvovrd*moist(i,k,j,P_QV)
          grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
          grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
      ENDDO

      mu = 0.
      sigma = 0.02
      a = -1.
      b = 1.
      call truncated_normal_ab_sample( mu, sigma, a, b, seed, rnd ) !RANDN(rnd)
      DO K = 1, 20
          !WRITE(*,*) 'PB', grid%pb(i,k,j)
          IF (grid%pb(i,k,j).LT.(grid%pb(i,1,j)-10000.)) EXIT

          rnd2 = rnd*EXP(FLOAT(1-K)/2)  ! same sign in one column, smoothes out
          if ((J.eq.jts).and.(I.eq.its)) then
              write(*,*) 'Wpert  k: ', k, ', rnd:', rnd2
          end if

          ! perturb vertical velocity
          grid%w_1(i,k,j) = grid%w_1(i,k,j)+rnd2
          grid%w_2(i,k,j) = grid%w_1(i,k,j)
      ENDDO

!  rebalance hydrostatically

      DO k  = 2,kte
        grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (grid%dnw(k-1))*(       &
                     ((grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))+(grid%c1h(k)*grid%mu_1(i,j)))*grid%al(i,k-1,j)+ &
                      (grid%c1h(k)*grid%mu_1(i,j))*grid%alb(i,k-1,j)  )

        grid%ph_2(i,k,j) = grid%ph_1(i,k,j)
        grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO

    ENDDO
  ENDDO
  ! -------------------------------------
  CASE (squall2d_x)
!  QSS thermal perturbation to kick off convection

  write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
  write(6,*) ' delt for perturbation ',delt

  DO J = jts, min(jde-1,jte)
!   yrad = config_flags%dy*float(j-nyc)/4000.
    yrad = 0.
    DO I = its, min(ide-1,ite)
      xrad = config_flags%dx*float(i-nxc)/4000.
!     xrad = 0.
      DO K = 1, kte-1

!  put in preturbation theta (bubble) and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

        zrad = 0.5*(grid%ph_1(i,k,j)+grid%ph_1(i,k+1,j)  &
                   +grid%phb(i,k,j)+grid%phb(i,k+1,j))/g
        zrad = (zrad-1500.)/1500.
        RAD=SQRT(xrad*xrad+yrad*yrad+zrad*zrad)
        IF(RAD <= 1.) THEN
           grid%t_1(i,k,j)=grid%t_1(i,k,j)+delt*COS(.5*PI*RAD)**2
           grid%t_2(i,k,j)=grid%t_1(i,k,j)
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
        ENDIF
      ENDDO

!  rebalance hydrostatically

      DO k  = 2,kte
        grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (grid%dnw(k-1))*(       &
                     ((grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))+(grid%c1h(k)*grid%mu_1(i,j)))*grid%al(i,k-1,j)+ &
                      (grid%c1h(k)*grid%mu_1(i,j))*grid%alb(i,k-1,j)  )

        grid%ph_2(i,k,j) = grid%ph_1(i,k,j)
        grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO

    ENDDO
  ENDDO
  CASE (squall2d_y)
!  QSS thermal perturbation to kick off convection

  write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
  write(6,*) ' delt for perturbation ',delt

  DO J = jts, min(jde-1,jte)
    yrad = config_flags%dy*float(j-nyc)/4000.
!   yrad = 0.
    DO I = its, min(ide-1,ite)
!     xrad = config_flags%dx*float(i-nxc)/4000.
      xrad = 0.
      DO K = 1, kte-1

!  put in preturbation theta (bubble) and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

        zrad = 0.5*(grid%ph_1(i,k,j)+grid%ph_1(i,k+1,j)  &
                   +grid%phb(i,k,j)+grid%phb(i,k+1,j))/g
        zrad = (zrad-1500.)/1500.
        RAD=SQRT(xrad*xrad+yrad*yrad+zrad*zrad)
        IF(RAD <= 1.) THEN
           grid%t_1(i,k,j)=grid%t_1(i,k,j)+delt*COS(.5*PI*RAD)**2
           grid%t_2(i,k,j)=grid%t_1(i,k,j)
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
        ENDIF
      ENDDO

!  rebalance hydrostatically

      DO k  = 2,kte
        grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (grid%dnw(k-1))*(       &
                     ((grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))+(grid%c1h(k)*grid%mu_1(i,j)))*grid%al(i,k-1,j)+ &
                      (grid%c1h(k)*grid%mu_1(i,j))*grid%alb(i,k-1,j)  )

        grid%ph_2(i,k,j) = grid%ph_1(i,k,j)
        grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO

    ENDDO
  ENDDO
  CASE (convrad)
!  random low-level thermal perturbation to kick off convection
        call random_seed
  write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
  write(6,*) ' delt for perturbation ',delt

  DO J = jts, min(jde-1,jte)
    DO I = its, min(ide-1,ite)
      DO K = 1, 10

        call RANDOM_NUMBER(rnd)
          grid%t_1(i,k,j)=grid%t_1(i,k,j)+delt*(rnd-0.5)
         !  grid%t_1(i,k,j)=grid%t_1(i,k,j)+delt*COS(.5*PI*RAD)**2
           grid%t_2(i,k,j)=grid%t_1(i,k,j)
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
      ENDDO

!  rebalance hydrostatically

      DO k  = 2,kte
        grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (grid%dnw(k-1))*(       &
                     ((grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))+(grid%c1h(k)*grid%mu_1(i,j)))*grid%al(i,k-1,j)+ &
                      (grid%c1h(k)*grid%mu_1(i,j))*grid%alb(i,k-1,j)  )
        grid%ph_2(i,k,j) = grid%ph_1(i,k,j)
        grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO

    ENDDO
  ENDDO
  CASE (grav2d_x)
!  cold bubble input  (from straka et al, IJNMF, vol 17, 1993 pp 1-22)

  t_min = grid%t_1(its,kts,jts)
  t_max = t_min
  u_mean = 00.

  xpos = config_flags%dx*nxc - u_mean*900.
  xposml = xpos - config_flags%dx*(ide-1)
  xpospl = xpos + config_flags%dx*(ide-1)

  DO J = jts, min(jde-1,jte)
    DO I = its, min(ide-1,ite)
!      xrad = config_flags%dx*float(i-nxc)/4000.  !  4000 meter horizontal radius
!                                    !  centered in the domain

       xrad = min( abs(config_flags%dx*float(i)-xpos),   &
                   abs(config_flags%dx*float(i)-xposml), &
                   abs(config_flags%dx*float(i)-xpospl))/4000.

      DO K = 1, kte-1

!  put in preturbation theta (bubble) and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

        zrad = 0.5*(grid%ph_1(i,k,j)+grid%ph_1(i,k+1,j)  &
                   +grid%phb(i,k,j)+grid%phb(i,k+1,j))/g
        zrad = (zrad-3000.)/2000. !  2000 meter vertical radius,
                                  !  centered at z=3000,
        RAD=SQRT(xrad*xrad+zrad*zrad)
        IF(RAD <= 1.) THEN

           !  perturbation temperature is 15 C, convert to potential temperature

           delt = -15.0 / ((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**rcp

           grid%T_1(i,k,j)=grid%T_1(i,k,j)+delt*(COS(PI*RAD)+1.0)/2.
           grid%T_2(i,k,j)=grid%T_1(i,k,j)
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
        ENDIF

        t_min = min(t_min, grid%t_1(i,k,j))
        t_max = max(t_max, grid%t_1(i,k,j))
      ENDDO

!  rebalance hydrostatically

      DO k  = 2,kte
        grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (grid%dnw(k-1))*(       &
                     ((grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))+(grid%c1h(k)*grid%mu_1(i,j)))*grid%al(i,k-1,j)+ &
                      (grid%c1h(k)*grid%mu_1(i,j))*grid%alb(i,k-1,j)  )

        grid%ph_2(i,k,j) = grid%ph_1(i,k,j)
        grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO

    ENDDO
  ENDDO

  write(6,*) ' min and max theta perturbation ',t_min,t_max

! -- end bubble insert
  CASE (b_wave)
!  thermal perturbation to kick off convection

  write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
  write(6,*) ' delt for perturbation ',tpbub

  DO J = jts, min(jde-1,jte)
    yrad = config_flags%dy*float(j-jde/2-1)/radbub
    DO I = its, min(ide-1,ite)
      xrad = float(i-1)/float(ide-ids)

      DO K = 1, kte-1

!  put in preturbation theta (bubble) and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

        zrad = 0.5*(grid%ph_1(i,k,j)+grid%ph_1(i,k+1,j)  &
                   +grid%phb(i,k,j)+grid%phb(i,k+1,j))/g
        zrad = (zrad-htbub)/radz
        RAD=SQRT(yrad*yrad+zrad*zrad)
        IF(RAD <= 1.) THEN
           tp = tpbub*cos(rad*piov2)*cos(rad*piov2)*cos(xrad*2*pi+pi)
           grid%t_1(i,k,j)=grid%t_1(i,k,j)+tp
           grid%t_2(i,k,j)=grid%t_1(i,k,j)
           qvf = 1. + rvovrd*grid%moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
        ENDIF
      ENDDO

!  rebalance hydrostatically

      DO k  = 2,kte
        grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (grid%dnw(k-1))*(       &
                     ((grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))+(grid%c1h(k)*grid%mu_1(i,j)))*grid%al(i,k-1,j)+ &
                      (grid%c1h(k)*grid%mu_1(i,j))*grid%alb(i,k-1,j)  )

        grid%ph_2(i,k,j) = grid%ph_1(i,k,j)
        grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO

    ENDDO
  ENDDO

  CASE (les)
!  thermal perturbation to kick off convection

  write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
  write(6,*) ' delt for perturbation ',delt

! For LES, change the initial random perturbations
! For 2D test, call randx outside I-loop
! For 3D runs, call randx inside both I-J loops

  DO J = jts, min(jde-1,jte)
!   yrad = config_flags%dy*float(j-nyc)/10000.
    yrad = 0.
    DO I = its, min(ide-1,ite)
!     xrad = config_flags%dx*float(i-nxc)/10000.
      xrad = 0.
      call random_number (randx)
      randx = randx - 0.5
!     DO K = 1, kte-1
      DO K = 1, 4

!  No bubbles for LES!
!  put in preturbation theta (bubble) and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

!       zrad = 0.5*(grid%ph_1(i,k,j)+grid%ph_1(i,k+1,j)  &
!                  +grid%phb(i,k,j)+grid%phb(i,k+1,j))/g
!       zrad = (zrad-1500.)/1500.
        zrad = 0.
        RAD=SQRT(xrad*xrad+yrad*yrad+zrad*zrad)
        IF(RAD <= 1.) THEN
!          grid%t_1(i,k,j)=grid%t_1(i,k,j)+delt*COS(.5*PI*RAD)**2
           grid%t_1(i,k,j)=grid%t_1(i,k,j)+ 0.1 *randx
           grid%t_2(i,k,j)=grid%t_1(i,k,j)
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
        ENDIF
      ENDDO

!  rebalance hydrostatically


      DO  kk  = 2,kte
        k = kk - 1
        grid%ph_1(i,kk,j) = grid%ph_1(i,kk-1,j) - (grid%dnw(kk-1))*(       &
                     ((grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))+(grid%c1h(k)*grid%mu_1(i,j)))*grid%al(i,kk-1,j)+ &
                      (grid%c1h(k)*grid%mu_1(i,j))*grid%alb(i,kk-1,j)  )

        grid%ph_2(i,kk,j) = grid%ph_1(i,kk,j)
        grid%ph0(i,kk,j) = grid%ph_1(i,kk,j) + grid%phb(i,kk,j)
      ENDDO

    ENDDO
  ENDDO

  END SELECT ideal_pert

#if 0
!===============

!  Test for resting atmosphere

  DO iter_loop = 1, 100
     DO J = jts, min(jde-1,jte)
        DO I = its, min(ide-1,ite)
            DO K = 2, kte-1
               z1 = (grid%phb(i,k+1,j)+grid%ph_1(i,k+1,j))/g
               z2 = (grid%phb(i,k  ,j)+grid%ph_1(i,k  ,j))/g
               z=(z1+z2)/2.
               IF ( z .LT. 10000 ) THEN
                  Nsq = 1.E-4
               ELSE
                  Nsq = 4.E-4
               END IF
               grid%t_1(i,k,j)=(grid%t_1(i,k-1,j)+t0) * EXP(Nsq / g * (z1-z2)) - t0
               grid%t_2(i,k,j)=grid%t_1(i,k,j)
               qvf = 1. + rvovrd*moist(i,k,j,P_QV)
               grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                                 (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
               grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
            ENDDO

            DO kk  = 2,kte
              k = kk - 1
              grid%ph_1(i,kk,j) = grid%ph_1(i,kk-1,j) - (grid%dnw(kk-1))*(       &
                           ((grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))+(grid%c1h(k)*grid%mu_1(i,j)))*grid%al(i,kk-1,j)+ &
                            (grid%c1h(k)*grid%mu_1(i,j))*grid%alb(i,kk-1,j)  )

              grid%ph_2(i,kk,j) = grid%ph_1(i,kk,j)
              grid%ph0(i,kk,j) = grid%ph_1(i,kk,j) + grid%phb(i,kk,j)
            ENDDO
        ENDDO
     ENDDO
  END DO

!===============
#endif
   IF ( wrf_dm_on_monitor() ) THEN
   k=1
   write(6,*) ' grid%mu_1 from comp ', (grid%c1h(k)*grid%mu_1(1,1))
   write(6,*) ' full state sounding from comp, ph, grid%p, grid%al, grid%t_1, qv '
   do k=1,kde-1
     write(6,'(i3,1x,5(1x,1pe10.3))') k, grid%ph_1(1,k,1)+grid%phb(1,k,1), &
                                      grid%p(1,k,1)+grid%pb(1,k,1), grid%alt(1,k,1), &
                                      grid%t_1(1,k,1)+t0, moist(1,k,1,P_QV)
   enddo

   write(6,*) ' pert state sounding from comp, grid%ph_1, pp, alp, grid%t_1, qv '
   do k=1,kde-1
     write(6,'(i3,1x,5(1x,1pe10.3))') k, grid%ph_1(1,k,1), &
                                      grid%p(1,k,1), grid%al(1,k,1), &
                                      grid%t_1(1,k,1), moist(1,k,1,P_QV)
   enddo
   ENDIF


  IF ( model_config_rec%ideal_case .EQ. b_wave )THEN
! set v (u already interpolated)

  DO J = jts, jte
  DO I = its, min(ide-1,ite)

    DO K = 1, kte
      grid%v_1(i,k,j) = 0.
      grid%v_2(i,k,j) = grid%v_1(i,k,j)
    ENDDO

  ENDDO
  ENDDO

!  fill out last i row for u

  DO J = jts, min(jde-1,jte)
  DO I = ite, ite

    DO K = 1, kte
      grid%u_1(i,k,j) = grid%u_1(its,k,j)
      grid%u_2(i,k,j) = grid%u_2(its,k,j)
    ENDDO

  ENDDO
  ENDDO

  ELSE
! interp v
  DO J = jts, jte
  DO I = its, min(ide-1,ite)

    IF (j == jds) THEN
      z_at_v = grid%phb(i,1,j)/g
    ELSE IF (j == jde) THEN
      z_at_v = grid%phb(i,1,j-1)/g
    ELSE
      z_at_v = 0.5*(grid%phb(i,1,j)+grid%phb(i,1,j-1))/g
    END IF

    p_surf = interp_0( p_in, zk, z_at_v, nl_in )

    DO K = 1, kte
      p_level = grid%c3h(k)*(p_surf - grid%p_top) + grid%c4h(k) + grid%p_top
      grid%v_1(i,k,j) = interp_0( v, p_in, p_level, nl_in )
      grid%v_2(i,k,j) = grid%v_1(i,k,j)
    ENDDO
  ENDDO
  ENDDO
! interp u
  DO J = jts, min(jde-1,jte)
  DO I = its, ite
    IF (i == ids) THEN
      z_at_u = grid%phb(i,1,j)/g
    ELSE IF (i == ide) THEN
      z_at_u = grid%phb(i-1,1,j)/g
    ELSE
      z_at_u = 0.5*(grid%phb(i,1,j)+grid%phb(i-1,1,j))/g
    END IF
    p_surf = interp_0( p_in, zk, z_at_u, nl_in )
    DO K = 1, kte
      p_level = grid%c3h(k)*(p_surf - grid%p_top) + grid%c4h(k) + grid%p_top
      grid%u_1(i,k,j) = interp_0( u, p_in, p_level, nl_in )
      grid%u_2(i,k,j) = grid%u_1(i,k,j)
    ENDDO

  ENDDO
  ENDDO
  ENDIF

!  set w

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte
  DO I = its, min(ide-1,ite)
    grid%w_1(i,k,j) = 0.
    grid%w_2(i,k,j) = 0.
  ENDDO
  ENDDO
  ENDDO

!  set a few more things

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte-1
  DO I = its, min(ide-1,ite)
    grid%h_diabatic(i,k,j) = 0.
  ENDDO
  ENDDO
  ENDDO

  IF ( wrf_dm_on_monitor() ) THEN
  DO k=1,kte-1
    grid%t_base(k) = grid%t_1(1,k,1)
    grid%qv_base(k) = moist(1,k,1,P_QV)
    grid%u_base(k) = grid%u_1(1,k,1)
    grid%v_base(k) = grid%v_1(1,k,1)
    grid%z_base(k) = 0.5*(grid%phb(1,k,1)+grid%phb(1,k+1,1)+grid%ph_1(1,k,1)+grid%ph_1(1,k+1,1))/g
  ENDDO
  ENDIF
  CALL wrf_dm_bcast_real( grid%t_base , kte )
  CALL wrf_dm_bcast_real( grid%qv_base , kte )
  CALL wrf_dm_bcast_real( grid%u_base , kte )
  CALL wrf_dm_bcast_real( grid%v_base , kte )
  CALL wrf_dm_bcast_real( grid%z_base , kte )

  ideal_surfacet: SELECT CASE ( model_config_rec%ideal_case )
  CASE(hill2d_x, squall2d_x, squall2d_y, grav2d_x, b_wave)
  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)
     thtmp   = grid%t_2(i,1,j)+t0
     ptmp    = grid%p(i,1,j)+grid%pb(i,1,j)
     temp(1) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = grid%t_2(i,2,j)+t0
     ptmp    = grid%p(i,2,j)+grid%pb(i,2,j)
     temp(2) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = grid%t_2(i,3,j)+t0
     ptmp    = grid%p(i,3,j)+grid%pb(i,3,j)
     temp(3) = thtmp * (ptmp/p1000mb)**rcp

     grid%tsk(I,J)=grid%cf1*temp(1)+grid%cf2*temp(2)+grid%cf3*temp(3)
     grid%tmn(I,J)=grid%tsk(I,J)-0.5
  ENDDO
  ENDDO

  CASE(seabreeze2d_x)
  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)
     thtmp   = grid%t_2(i,1,j)+t0
     ptmp    = grid%p(i,1,j)+grid%pb(i,1,j)
     temp(1) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = grid%t_2(i,2,j)+t0
     ptmp    = grid%p(i,2,j)+grid%pb(i,2,j)
     temp(2) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = grid%t_2(i,3,j)+t0
     ptmp    = grid%p(i,3,j)+grid%pb(i,3,j)
     temp(3) = thtmp * (ptmp/p1000mb)**rcp

!    grid%tsk(I,J)=grid%cf1*temp(1)+grid%cf2*temp(2)+grid%cf3*temp(3)
     grid%tmn(I,J)=grid%tsk(I,J)-0.5
  ENDDO
  ENDDO

  CASE(convrad, quarter_ss)
  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)
         grid%tsk(i,j) = theta_surf * (p_surf/p1000mb)**rcp
         grid%tmn(i,j) = grid%tsk(i,j)
  ENDDO
  ENDDO

  CASE(les)
  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)
         grid%tsk(i,j) = theta_surf * (p_surf/p1000mb)**rcp
         grid%tmn(i,j) = grid%tsk(i,j)-0.5
  ENDDO
  ENDDO

  END SELECT ideal_surfacet


  trajectories: SELECT CASE ( model_config_rec%ideal_case )
  CASE (quarter_ss)
  !  Template for initializing trajectories.  The i, j, and k starting locations
  !  are specified.  Right now, a small plane in the middle of the domain is
  !  selected.

  grid%traj_i    = -9999
  grid%traj_j    = -9999
  grid%traj_k    = -9999
  grid%traj_lat  = -9999
  grid%traj_long = -9999

  IF (config_flags%num_traj .gt. 0 .and. config_flags%traj_opt .gt. 0) THEN
     icount = 1
     DO j = (jde + jds)/2 - 2, (jde + jds)/2 + 2, 1
        DO i = (ide + ids)/2 - 2, (ide + ids)/2 + 2, 1
           IF ( its .LE. i    .and. ite .GE. i   .and.  jts .LE. j    .and. jte .GE. j ) THEN
              grid%traj_i   (icount) = i
              grid%traj_j   (icount) = j
              grid%traj_k   (icount) = 10
              grid%traj_lat (icount) = grid%xlat(i,j)
              grid%traj_long(icount) = grid%xlong(i,j)
           END IF

#ifdef DM_PARALLEL
           grid%traj_i   (icount) = wrf_dm_max_real ( grid%traj_i   (icount) )
           grid%traj_j   (icount) = wrf_dm_max_real ( grid%traj_j   (icount) )
           grid%traj_k   (icount) = wrf_dm_max_real ( grid%traj_k   (icount) )
           grid%traj_lat (icount) = wrf_dm_max_real ( grid%traj_lat (icount) )
           grid%traj_long(icount) = wrf_dm_max_real ( grid%traj_long(icount) )
#endif

           icount = icount + 1
           IF (icount .GT. config_flags%num_traj) THEN
              EXIT
           END IF
        END DO
     END DO
  END IF
  END SELECT trajectories

  tracers: SELECT CASE ( config_flags%tracer_opt )
  CASE (tracer_test1)

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte-1
  DO I = its, min(ide-1,ite)
    grid%h_diabatic(i,k,j) = 0.
    if(k.eq.kts)tracer(i,k,j,p_tr17_1)=1.
    if(k.eq.kts.and.grid%xland(i,j).lt.1.5)tracer(i,k,j,p_tr17_2)=1.
    if(k.eq.kts.and.grid%xland(i,j).gt.1.5)tracer(i,k,j,p_tr17_3)=1.
    if(k.le.5)tracer(i,k,j,p_tr17_4)=1.
    if(k.le.5.and.grid%xland(i,j).lt.1.5)tracer(i,k,j,p_tr17_5)=1.
    if(k.le.5.and.grid%xland(i,j).gt.1.5)tracer(i,k,j,p_tr17_6)=1.
    if(k.le.10)tracer(i,k,j,p_tr17_7)=1.
    if(k.le.10.and.k.gt.5)tracer(i,k,j,p_tr17_8)=1.
  ENDDO
  ENDDO
  ENDDO

  END SELECT tracers

      !  Save the dry perturbation potential temperature.

      DO j = jts, min(jde-1,jte)
         DO k = kts, kte
            DO i = its, min(ide-1,ite)
               grid%th_phy_m_t0(i,k,j) = grid%t_2(i,k,j)
            END DO
         END DO
      END DO

  !  Turn dry potential temperature into moist potential temperature
      !  at the very end of this routine
      !  This field will be in the model IC and and used to construct the
      !  BC file.

      IF ( ( config_flags%use_theta_m .EQ. 1 ) .AND. (P_Qv .GE. PARAM_FIRST_SCALAR) ) THEN
      DO J  = jts, min(jde-1,jte)
         DO K = kts, kte-1
            DO I = its, min(ide-1,ite)
               grid%t_2(i,k,j) = ( grid%t_2(i,k,j) + T0 ) * (1. + (R_v/R_d) * moist(i,k,j,p_qv)) - T0
            END DO
         END DO
      END DO
      ENDIF

  RETURN

 END SUBROUTINE init_domain_rk

   SUBROUTINE init_module_initialize
   END SUBROUTINE init_module_initialize

!---------------------------------------------------------------------

!  test driver for get_sounding
!
!      implicit none
!      integer n
!      parameter(n = 1000)
!      real zk(n),p(n),theta(n),rho(n),u(n),v(n),qv(n),pd(n)
!      logical dry
!      integer nl,k
!
!      dry = .false.
!      dry = .true.
!      call get_sounding( zk, p, pd, theta, rho, u, v, qv, dry, n, nl )
!      write(6,*) ' input levels ',nl
!      write(6,*) ' sounding '
!      write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
!      do k=1,nl
!        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), pd(k), theta(k), rho(k), u(k), v(k), qv(k)
!      enddo
!      end
!
!---------------------------------------------------------------------------

      subroutine get_sounding( zk, p, p_dry, theta, rho, &
                               u, v, qv, dry, nl_max, nl_in, th_surf )
      implicit none

      integer nl_max, nl_in
      real zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max), &
           u(nl_max), v(nl_max), qv(nl_max), p_dry(nl_max)
      logical dry

      integer n, iz
      parameter(n=1000)
      logical debug
      parameter( debug = .true.)
      character*256 message

! input sounding data

      real p_surf, th_surf, qv_surf
      real pi_surf, pi(n)
      real h_input(n), th_input(n), qv_input(n), u_input(n), v_input(n)

! diagnostics

      real rho_surf, p_input(n), rho_input(n)
      real pm_input(n)  !  this are for full moist sounding

! local data

      real r
      parameter (r = r_d)
      integer k, it, nl
      real qvf, qvf1, dz

!  first, read the sounding

      call read_sounding( p_surf, th_surf, qv_surf, &
                          h_input, th_input, qv_input, u_input, v_input,n, nl, debug )

!        iz = 1
!        do k=2,nl
!          if(h_input(k) .lt. 12000.) iz = k
!        enddo
!        write(6,*) " tropopause ",iz,h_input(iz)
!        if(dry) then
!        write(6,*) ' nl is ',nl
!        do k=1,nl
!          th_input(k) = th_input(k)+10.+10*float(k)/nl
!        enddo
!        write(6,*) ' finished adjusting theta '
!        endif

!        do k=1,nl
!          u_input(k) = 2*u_input(k)
!        enddo
!
!      end if

      if(dry) then
       do k=1,nl
         qv_input(k) = 0.
       enddo
      endif

      if(debug) write(6,*) ' number of input levels = ',nl

        nl_in = nl
        if(nl_in .gt. nl_max ) then
          write(6,*) ' too many levels for input arrays ',nl_in,nl_max
          call wrf_error_fatal ( ' too many levels for input arrays ' )
        end if

!  compute diagnostics,
!  first, convert qv(g/kg) to qv(g/g)

      do k=1,nl
        qv_input(k) = 0.001*qv_input(k)
      enddo

      p_surf = 100.*p_surf  ! convert to pascals
      qvf = 1. + rvovrd*qv_input(1)
      rho_surf = 1./((r/p1000mb)*th_surf*qvf*((p_surf/p1000mb)**cvpm))
      pi_surf = (p_surf/p1000mb)**(r/cp)

      if(debug) then
        write(6,*) ' surface density is ',rho_surf
        write(6,*) ' surface pi is      ',pi_surf
      end if


!  integrate moist sounding hydrostatically, starting from the
!  specified surface pressure
!  -> first, integrate from surface to lowest level

          qvf = 1. + rvovrd*qv_input(1)
          qvf1 = 1. + qv_input(1)
          rho_input(1) = rho_surf
          dz = h_input(1) - 489.   ! @lkugler: first entry is not layer depth but absolute height-> subtract surface altitude
          do it=1,10
            pm_input(1) = p_surf &
                    - 0.5*dz*(rho_surf+rho_input(1))*g*qvf1
            rho_input(1) = 1./((r/p1000mb)*th_input(1)*qvf*((pm_input(1)/p1000mb)**cvpm))
          enddo

! integrate up the column

          do k=2,nl
            rho_input(k) = rho_input(k-1)
            dz = h_input(k)-h_input(k-1)
            qvf1 = 0.5*(2.+(qv_input(k-1)+qv_input(k)))
            qvf = 1. + rvovrd*qv_input(k)   ! qv is in g/kg here

            do it=1,10
              pm_input(k) = pm_input(k-1) &
                      - 0.5*dz*(rho_input(k)+rho_input(k-1))*g*qvf1
              IF(pm_input(k) .LE. 0. )THEN
                CALL wrf_message("Integrated pressure has gone negative - too cold for chosen height")
                WRITE(message,*)'k,pm_input(k),h_input(k),th_input(k) = ',k,pm_input(k),h_input(k),th_input(k)
                CALL wrf_error_fatal ( message )
              ENDIF
              rho_input(k) = 1./((r/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm))
            enddo
          enddo

!  we have the moist sounding

!  next, compute the dry sounding using p at the highest level from the
!  moist sounding and integrating down.

        p_input(nl) = pm_input(nl)

          do k=nl-1,1,-1
            dz = h_input(k+1)-h_input(k)
            p_input(k) = p_input(k+1) + 0.5*dz*(rho_input(k)+rho_input(k+1))*g
          enddo


        do k=1,nl

          zk(k) = h_input(k)
          p(k) = pm_input(k)
          p_dry(k) = p_input(k)
          theta(k) = th_input(k)
          rho(k) = rho_input(k)
          u(k) = u_input(k)
          v(k) = v_input(k)
          qv(k) = qv_input(k)

        enddo

     if(debug) then
      write(6,*) ' sounding '
      write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
      do k=1,nl
        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), p_dry(k), theta(k), rho(k), u(k), v(k), qv(k)
      enddo

     end if

      end subroutine get_sounding

!-------------------------------------------------------

      subroutine read_sounding( ps,ts,qvs,h,th,qv,u,v,n,nl,debug )
      implicit none
      integer n,nl
      real ps,ts,qvs,h(n),th(n),qv(n),u(n),v(n)
      logical end_of_file
      logical debug

      integer k

      open(unit=10,file='input_sounding',form='formatted',status='old')
      rewind(10)
      read(10,*) ps, ts, qvs
      if(debug) then
        write(6,*) ' input sounding surface parameters '
        write(6,*) ' surface pressure (mb) ',ps
        write(6,*) ' surface pot. temp (K) ',ts
        write(6,*) ' surface mixing ratio (g/kg) ',qvs
      end if

      end_of_file = .false.
      k = 0

      do while (.not. end_of_file)

        read(10,*,end=100) h(k+1), th(k+1), qv(k+1), u(k+1), v(k+1)
        k = k+1
        if(debug) write(6,'(1x,i3,5(1x,e10.3))') k, h(k), th(k), qv(k), u(k), v(k)
        go to 110
 100    end_of_file = .true.
 110    continue
      enddo

      nl = k

      close(unit=10,status = 'keep')

      end subroutine read_sounding

!------------------------------------------------------------------

    subroutine get_sounding_b_wave( zk, p, p_dry, theta, rho,       &
                             u, v, qv, dry, nl_max, nl_in,  &
                             u_jet, rho_jet, th_jet, z_jet, &
                             nz_jet, ny_jet, j_point, debug )
    implicit none

    integer nl_max, nl_in
    real zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max), &
         u(nl_max), v(nl_max), qv(nl_max), p_dry(nl_max)
    logical dry

    integer nz_jet, ny_jet, j_point
    real, dimension(nz_jet, ny_jet) :: u_jet, rho_jet, th_jet, z_jet

    integer n
    parameter(n=1000)
    logical debug

! input sounding data

    real p_surf, th_surf, qv_surf
    real pi_surf, pi(n)
    real h_input(n), th_input(n), qv_input(n), u_input(n), v_input(n)

! diagnostics

    real rho_surf, p_input(n), rho_input(n)
    real pm_input(n)  !  this are for full moist sounding

! local data

    real r
    parameter (r = r_d)
    integer k, it, nl
    real qvf, qvf1, dz

!  first, read the sounding

!    call read_sounding( p_surf, th_surf, qv_surf, &
!                          h_input, th_input, qv_input, u_input, v_input,n, nl, debug )

   call calc_jet_sounding( p_surf, th_surf, qv_surf,                             &
                           h_input, th_input, qv_input, u_input, v_input,        &
                           n, nl, debug, u_jet, rho_jet, th_jet, z_jet, j_point, &
                           nz_jet, ny_jet, dry                                  )

   nl = nz_jet

    if(dry) then
     do k=1,nl
       qv_input(k) = 0.
     enddo
    endif

    if(debug) write(6,*) ' number of input levels = ',nl

      nl_in = nl
      if(nl_in .gt. nl_max ) then
        write(6,*) ' too many levels for input arrays ',nl_in,nl_max
        call wrf_error_fatal ( ' too many levels for input arrays ' )
      end if

!  compute diagnostics,
!  first, convert qv(g/kg) to qv(g/g)
!
!      do k=1,nl
!        qv_input(k) = 0.001*qv_input(k)
!      enddo
!      p_surf = 100.*p_surf  ! convert to pascals

    qvf = 1. + rvovrd*qv_input(1)
    rho_surf = 1./((r/p1000mb)*th_surf*qvf*((p_surf/p1000mb)**cvpm))
    pi_surf = (p_surf/p1000mb)**(r/cp)

    if(debug) then
      write(6,*) ' surface density is ',rho_surf
      write(6,*) ' surface pi is    ',pi_surf
    end if


!  integrate moist sounding hydrostatically, starting from the
!  specified surface pressure
!  -> first, integrate from surface to lowest level

        qvf = 1. + rvovrd*qv_input(1)
        qvf1 = 1. + qv_input(1)
        rho_input(1) = rho_surf
        dz = h_input(1)
        do it=1,10
          pm_input(1) = p_surf &
                  - 0.5*dz*(rho_surf+rho_input(1))*g*qvf1
          rho_input(1) = 1./((r/p1000mb)*th_input(1)*qvf*((pm_input(1)/p1000mb)**cvpm))
        enddo

! integrate up the column

        do k=2,nl
          rho_input(k) = rho_input(k-1)
          dz = h_input(k)-h_input(k-1)
          qvf1 = 0.5*(2.+(qv_input(k-1)+qv_input(k)))
          qvf = 1. + rvovrd*qv_input(k)   ! qv is in g/kg here

          do it=1,10
            pm_input(k) = pm_input(k-1) &
                    - 0.5*dz*(rho_input(k)+rho_input(k-1))*g*qvf1
            rho_input(k) = 1./((r/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm))
          enddo
        enddo

!  we have the moist sounding

!  next, compute the dry sounding using p at the highest level from the
!  moist sounding and integrating down.

        p_input(nl) = pm_input(nl)

          do k=nl-1,1,-1
            dz = h_input(k+1)-h_input(k)
            p_input(k) = p_input(k+1) + 0.5*dz*(rho_input(k)+rho_input(k+1))*g
          enddo


        do k=1,nl

          zk(k) = h_input(k)
          p(k) = pm_input(k)
          p_dry(k) = p_input(k)
          theta(k) = th_input(k)
          rho(k) = rho_input(k)
          u(k) = u_input(k)
          v(k) = v_input(k)
          qv(k) = qv_input(k)

        enddo

     if(debug) then
      write(6,*) ' sounding '
      write(6,*) '  k  height(m)  press (Pa)   pd(Pa)   theta (K)  den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
      do k=1,nl
        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), p_dry(k), theta(k), rho(k), u(k), v(k), qv(k)
      enddo

     end if

     end subroutine get_sounding_b_wave

!------------------------------------------------------------------

  subroutine calc_jet_sounding( p_surf, th_surf, qv_surf,      &
                                h, th, qv, u, v, n, nl, debug, &
                                u_jet, rho_jet, th_jet, z_jet, &
                                jp, nz_jet, ny_jet, dry       )
  implicit none
  integer :: n, nl, jp, nz_jet, ny_jet

  real, dimension(nz_jet, ny_jet) :: u_jet, rho_jet, th_jet, z_jet
  real, dimension(n) :: h,th,qv,u,v
  real :: p_surf, th_surf, qv_surf
  logical :: debug, dry

  real, dimension(1:nz_jet) :: rho, rel_hum, p
  integer :: k

!  some local stuff

  real :: tmppi, es, qvs, temperature

!  get sounding from column jp

   do k=1,nz_jet
     h(k)  = z_jet(k,jp)
     th(k) = th_jet(k,jp)
     qv(k) = 0.
     rho(k) = rho_jet(k,jp)
     u(k) = u_jet(k,jp)
     v(k) = 0.
   enddo

   if (.not.dry) then
     DO k=1,nz_jet
       if(h(k) .gt. 8000.) then
         rel_hum(k)=0.1
       else
         rel_hum(k)=(1.-0.90*(h(k)/8000.)**1.25)
       end if
       rel_hum(k) = min(0.7,rel_hum(k))
     ENDDO
   else
     do k=1,nz_jet
       rel_hum(k) = 0.
     enddo
   endif

!  next, compute pressure

   do k=1,nz_jet
     p(k) = p1000mb*(R_d*rho(k)*th(k)/p1000mb)**cpovcv
   enddo

!  here we adjust for fixed moisture profile

     IF (.not.dry)  THEN

!  here we assume the input theta is th_v, so we reset theta accordingly

       DO k=1,nz_jet
         tmppi=(p(k)/p1000mb)**rcp
         temperature = tmppi*th(k)
         if (temperature .gt. svpt0) then
            es  = 1000.*svp1*exp(svp2*(temperature-svpt0)/(temperature-svp3))
            qvs = ep_2*es/(p(k)-es)
         else
            es  = 1000.*svp1*exp( 21.8745584*(temperature-273.16)/(temperature-7.66) )
            qvs = ep_2*es/(p(k)-es)
         endif
         qv(k) = rel_hum(k)*qvs
         th(k) = th(k)/(1.+.61*qv(k))
       ENDDO

     ENDIF

!  finally, set the surface data. We'll just do a simple extrapolation

   p_surf = 1.5*p(1) - 0.5*p(2)
   th_surf = 1.5*th(1) - 0.5*th(2)
   qv_surf = 1.5*qv(1) - 0.5*qv(2)

   end subroutine calc_jet_sounding

!---------------------------------------------------------------------

 SUBROUTINE read_input_jet( u, r, t, zk, nz, ny )
 implicit none

 integer, intent(in) :: nz,ny
 real, dimension(nz,ny), intent(out) :: u,r,t,zk
 integer :: ny_in, nz_in, j,k
 real, dimension(ny,nz) :: field_in
 character*256 message

! this code assumes it is called on processor 0 only

   OPEN(unit=10, file='input_jet', form='unformatted', status='old' )
   REWIND(10)
   read(10) ny_in,nz_in
   if((ny_in /= ny ) .or. (nz_in /= nz)) then
     write(message,*) ' error in input jet dimensions '
     CALL wrf_message (message)
     write(message,*) ' ny, ny_input, nz, nz_input ', ny, ny_in, nz,nz_in
     CALL wrf_message (message)
     write(message,*) ' error exit '
     CALL wrf_message (message)
     call wrf_error_fatal ( ' error in input jet dimensions ' )
   end if
   read(10) field_in
   do j=1,ny
   do k=1,nz
     u(k,j) = field_in(j,k)
   enddo
   enddo
   read(10) field_in
   do j=1,ny
   do k=1,nz
     t(k,j) = field_in(j,k)
   enddo
   enddo

   read(10) field_in
   do j=1,ny
   do k=1,nz
     r(k,j) = field_in(j,k)
   enddo
   enddo

   do j=1,ny
   do k=1,nz
     zk(k,j) = 125. + 250.*float(k-1)
   enddo
   enddo

 end subroutine read_input_jet


END MODULE module_initialize_ideal
