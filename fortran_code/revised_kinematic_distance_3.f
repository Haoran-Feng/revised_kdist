      program revised_kinematic_distance

c     Calculates kinematic distances using the revised prescription given in 
c     Reid et al 2009 (Paper VI):  
c     Version 2 (projects Vlsr to Galactic Plane)
c     Version 3 (fixes bug: used Theta instead of Theta_0 in rootterm)

      implicit real*8 (a-h,o-z)

      character*48  file_name

      character*12  src
      character*1   c1
      equivalence  (src, c1)

c     Logical unit numbers...
      lu_out       = 6
      lu_dat       = 8
      lu_control   = 9

      pi = 4.d0 * atan(1.d0)
      deg_rad = pi/180.d0

c     =============================================================
c     Get Galaxy, Solar Motion, and average source parameters from ascii 
c     text file named "parameter_file.inp".  See example file for details

      file_name = 'parameter_file.inp' 
      call galaxy_parameters ( lu_control, file_name, 
     +         Ro, To, dTdr, Uo, Vo, Wo, Us, Vs, Ws, err_v_lsr )

      write (lu_out,7000)
 7000 format('! Program: revised_kineamtic_distance.f:',
     +       ' Version 3; M. Reid; November 2011')
      write (lu_out,7005) Ro, To, dTdr
 7005 format('! Assumed Galactic paramters:  Ro, To, dTdr',3f7.2)
      write (lu_out,7010) Uo, Vo, Wo
 7010 format('! Assumed Solar Motion parms:  Uo, Vo, Wo  ',3f7.2)
      write (lu_out,7020) Us, Vs, Ws
 7020 format('! Assumed Source Motion parms: Us, Vs, Ws  ',3f7.2)
      write (lu_out,7030) err_v_lsr
 7030 format('! Assumed uncertainties in V_lsr           ',f7.2)
      write (lu_out,7090)
 7090 format('! Source     Gal Long  Gal Lat    V_lsr     V_rev  ',
     +       '  Rev. D_k     +/-  ',/
     +       '!              (deg)    (deg)    (km/s)    (km/s)  ',
     +       '   (kpc)      (kpc)')

c     ============================================================
c     Get input source data from ascii text file "source_file.dat"
c     See example file for details

      file_name = 'source_file.dat'
      call open_ascii_file ( lu_dat, file_name )

      ieof    = 0
      do while ( ieof.ge.0 )

c        Input data for a source ...
         read(lu_dat,*,iostat=ieof) src, 
     +             ra_hhmmss, dec_ddmmss, v_lsr, farnear

c        Check for end-of-file mark or "!" (skip this line)...
         if ( ieof.ge.0 .and. c1.ne.'!' ) then

c           Convert coordinates to radians and then to hours, deg
            call hmsrad ( ra_hhmmss, ra_rad )
            call dmsrad (dec_ddmmss,dec_rad )
            ra = ra_rad * 12.d0  / pi        ! hours
            dec=dec_rad * 180.d0 / pi        ! degrees

c           Get Galactic coordinates of source
            call radec_to_galactic ( ra, dec,  ell, bee )

c           Calculate "adjusted" Vlsr, including effects of (Us,Vs).
c           Since applying (Us,Vs) needs Dk, so do both at same time.
            call calc_Dk ( lu_out,
     +             Ro, To, dTdr, Uo, Vo, Wo, Us, Vs, Ws,
     +             src, ra, dec, farnear, v_lsr,
     +             v_lsr_rev, Dk_best )

c           -------------------------------------------------------
c           Estimation of (asymmetric) uncertainties based only on
c           the uncertainty in Vlsr supplied in "parameter_file.inp".
c           Assumes Galactic & Source parameters have no uncertainties.

c           Check if v_lsr exceeds tangent point speed.  If so,
c           use tangent point speed as center for velocity offsets
c           used for error estimate.

c           Calculate tangent point speed (T_tp)
            sin_l = sin( ell*deg_rad )
            R_tp = Ro * sin_l                   ! kpc
            T_tp = To + dTdR * (R_tp - Ro)      ! km/s

c           Check of v_lsr_rev exceeds tangent point speed
            del_v = 0.d0
            v_tp  = -9999.9d0
            if ( ell .lt. 90.d0  ) then         ! quadrant 1
               v_tp = +T_tp - To*sin_l          ! km/s 
               if ( v_lsr_rev .gt. v_tp ) del_v = v_lsr_rev - v_tp
            endif
            if ( ell .gt. 270.d0 ) then         ! quadrant 4
               v_tp = -T_tp - To*sin_l          ! km/s 
               if ( v_lsr_rev .lt. v_tp ) del_v = v_lsr_rev - v_tp
            endif
c           Adjust v_lsr (if needed) so as not to exceed tangent point
            v_lsr_adj = v_lsr - del_v

c           Get lower velocity distance uncertainty
            v_lsr_low = v_lsr_adj - err_v_lsr 
            call calc_Dk ( lu_out,
     +             Ro, To, dTdr, Uo, Vo, Wo, Us, Vs, Ws,
     +             src, ra, dec, farnear, v_lsr_low,
     +             v_lsr_rev_low, D1 )
c           Flag if close to tangent point velocity
            if ( abs(v_lsr_rev_low - v_tp) .lt. err_v_lsr ) D1=Dk_best

c           Get higher velocity distance uncertainty
            v_lsr_high = v_lsr_adj + err_v_lsr 
            call calc_Dk ( lu_out,
     +             Ro, To, dTdr, Uo, Vo, Wo, Us, Vs, Ws,
     +             src, ra, dec, farnear, v_lsr_high,
     +             v_lsr_rev_high, D2 )
c           Flag if close to tangent point velocity
            if ( abs(v_lsr_rev_high - v_tp) .lt. err_v_lsr ) D2=Dk_best

c           Calculate + (high) and - (low) errors
            if ( D2 .gt. D1 ) then 
                 d_err_low  = D1  - Dk_best
                 d_err_high = D2  - Dk_best
               else
                 d_err_low  = D2  - Dk_best
                 d_err_high = D1  - Dk_best
            endif

c           If flagged distance (error=0), use other error estimate
            if ( abs(d_err_low) .eq. 0.d0 ) d_err_low  = -d_err_high
            if ( d_err_high     .eq. 0.d0 ) d_err_high = -d_err_low

c           Check for pathalogical cases...
            if ( Dk_best .le. 0.0 ) then 
               Dk_best    = 0.0
               d_err_low  = 0.0
               d_err_high = 0.0
            endif
            if ( Dk_best+d_err_low .lt. 0.0 ) d_err_low = -Dk_best

            write (lu_out,9000) src, ell, bee, v_lsr, v_lsr_rev,
     +                       Dk_best, d_err_high, d_err_low
 9000       format(1x,a12,2f8.3,2f10.1,f10.2,2f7.2)

         endif

      enddo

      stop
      end

c======================================================================
	subroutine open_ascii_file ( lu_in, ascii_file )

c	Opens ascii text file 
c       Reads and skips comments (if first character is "!") 
c       Leaves file open, ready to read first data line

	character*48       ascii_file

	character*80       comment
	character*1        c_1
	equivalence       (c_1, comment)

C	Open input ascii file...
	open (unit=lu_in, file=ascii_file, status='old')

C       Read all comment lines (but stop at data)
	c_1 = '!'
	do while ( c_1 .eq. '!' )

	   read  (lu_in,1000) comment
 1000	   format(a80)

	   if ( c_1 .ne. '!' ) backspace (unit=lu_in)
	       
	enddo

	return
	end

c===================================================================
      subroutine galaxy_parameters( lu_control, file_name,
     +           Ro, To, dTdr, Uo, Vo, Wo, Us, Vs, Ws, err_v_lsr )

c     Opens file and reads parameters needed for revised kinematic distance
 
      implicit real*8 (a-h,o-z)

      character*48    file_name

      call open_ascii_file ( lu_control, file_name )

c     Galactic structure parameters
      read (lu_control,*) Ro       ! Ro [kpc]
      read (lu_control,*) To       ! Theta_0 [km/s]
      read (lu_control,*) dTdr     ! dTheta/dr [km/s/kpc]

c     Solar Motion parameters
      read (lu_control,*) Uo       ! Uo toward G.C. [km/s]
      read (lu_control,*) Vo       ! Vo toward Gal. rotation [km/s]
      read (lu_control,*) Wo       ! Wo toward North Gal.Pole [km/s]

c     Source peculiar motion in its local Galactocentric frame
c     (G.C.is as viewed by source, not Sun)
      read (lu_control,*) Us       ! Us toward G.C. [km/s]
      read (lu_control,*) Vs       ! Vs toward Gal. rotation [km/s]
      read (lu_control,*) Ws       ! Ws toward North Gal.Pole [km/s]

c     Uncertainty in V_lsr (used to estimate uncertainty in distance)
      read (lu_control,*) err_v_lsr ! [km/s]
        
      close(unit=lu_control)

      return
      end

c===================================================================
      subroutine calc_Dk ( lu_out,
     +        Ro, To, dTdr, Uo, Vo, Wo, Us, Vs, Ws,  
     +        src, ra, dec, farnear, v_lsr,
     +        v_lsr_rev, Dk )

c     Calculate revised Vlsr by converting standard Vlsr back to
c     heliocentric, apply modern Solar Motion values (Uo,Vo,Wo), and 
c     remove effects of average source non-ciruclar motion (Us,Vs,Ws).  
c     Then calculate kinematic distance using the linear rotation 
c     curve specified by Ro, To, and dTdR

      implicit real*8 (a-h,o-z)

      character*12   src

      pi     = 4.d0 * atan(1.d0)
      deg_rad = pi/180.d0

      n_iter_max = 1000

c     ============================================================
c     Get source galactic coordinates...

      call radec_to_galactic (ra, dec, gal_long, gal_lat)

      gal_long_rad = gal_long * deg_rad       ! radians
      cos_l = cos( gal_long_rad )
      sin_l = sin( gal_long_rad )

      gal_lat_rad  = gal_lat * deg_rad        ! radians
      cos_b = cos( gal_lat_rad )
      sin_b = sin( gal_lat_rad )

c     -------------------------------------------------------------
c     Convert to true Heliocentric frame (v_rad)
c     Add back Sun's peculiar motion to radial velocity
c     Use old Standard Solar Motion (Sun moves at 20 km/s toward
c     18h, +30deg in 1900 coordinates), since this has been used to 
c     define V_lsr:

      Uo_IAU = 10.27d0            ! km/s precessed to J2000
      Vo_IAU = 15.32d0
      Wo_IAU =  7.74d0

      v_helio = v_lsr - (Vo_IAU*sin_l + Uo_IAU*cos_l)*cos_b 
     +                -  Wo_IAU*sin_b

c     -------------------------------------------------------------
c     Make "new" V(LSR) using best Solar Motion
c      eg, Hipparcos (Dehnen & Binney 1997) gives 
c      Uo = 10.00d0                ! km/s
c      Vo =  5.25d0
c      Wo =  7.17d0
  
      v_newlsr = v_helio + (Vo*sin_l + Uo*cos_l)*cos_b 
     +                   +  Wo*sin_b

c     --------------------------------------------------------
c     Remove effects of common peculiar motions specified in
c     "parameter_file.inp"

c     If dTdr.ne.0, need to know distance to get source Galactocentric 
c     radius (Rs) to evalute rotation curve.   So must iterate...
      n_iter =  0
      del_d  = 99.d0
      Dk     =  3.d0

      do while ( del_d.gt.0.01d0 .and. n_iter.lt.n_iter_max )

c        Save old value of kinematic distance
         Dk_old = Dk

c        Calculate "gamma" angle and projected Galactocentric radius 
         d_proj = Dk * cos_b                     ! kpc in Gal Plane
         r_sq   = Ro**2 + d_proj**2 - 2.d0*Ro * d_proj * cos_l
         r_proj = sqrt( r_sq )                   ! kpc in Gal Plane

c        Calculate Galactocentric longitude (beta in paper)...
         sin_beta =   d_proj  * sin_l     / r_proj
         cos_beta = ( Ro - d_proj*cos_l ) / r_proj
         beta     = atan2( sin_beta, cos_beta )  ! radians
         beta_deg = beta / deg_rad               ! deg

c        Calculate Sun-Maser-GC angle...
         gamma = pi - gal_long_rad - beta        ! radians
         cos_gamma = cos( gamma )
         sin_gamma = sin( gamma )

         v_fixed = v_newlsr - (Vs*sin_gamma - Us*cos_gamma)*cos_b
     +                      -  Ws*sin_b          ! km/s

c        -----------------------------------------------------------------
c        Calculate a kinematic distance using best Ro, To and dTdr

         Rs = r_proj
         V_proj = v_fixed * cos_b
         call kinematic_distance ( v_proj, gal_long, Ro, To, 
     +                             r_proj, dTdr,
     +                             D_near, D_far )

         Dk = D_near
         if ( farnear .ne. 0.d0 ) Dk = D_far

c        Ignore "farnear" flag if one of the values is zero
         if ( D_near .le. 0.d0 .and. D_far  .gt. 0.d0 ) Dk = D_far
         if ( D_far  .le. 0.d0 .and. D_near .gt. 0.d0 ) Dk = D_near

         del_d = abs( Dk - Dk_old )
         n_iter = n_iter + 1

      enddo

      v_lsr_rev = v_fixed

      return
      end

c     =======================================================
      subroutine kinematic_distance ( Vlsr, gal_long, Ro, To, 
     +                                Rs, dTdr,
     +                                D_near, D_far )

c     Caluclate kinematic distance given Vlsr (projected in Gal plane)
c     and information required to construct the kinematic model.  
c     Returns both near and far distances (when appropriate)

      implicit real*8 (a-h,o-z)

      pi     = 4.0d0*atan(1.0d0)	
      deg_rad = pi/180.d0	     ! convert degrees to radians

      glongrad = gal_long*deg_rad 

      cos_l = cos(glongrad)
      sin_l = sin(glongrad)

      Rosinl = Ro * sin_l
      Rocosl = Ro * cos_l

c     Non-flat rotation curve allowed...
      Theta = To + dTdr*(Rs - Ro)
      Tosinl= To    * sin_l
      Tsinl = Theta * sin_l

      rootterm = Rocosl**2 + ( Tsinl / ( Tosinl/Ro + Vlsr/Ro ) )**2 
     +         - Ro**2
      if ( rootterm .lt. 0 ) rootterm = 0.d0

      if ( gal_long.ge.0.d0   .and. gal_long.lt.90.d0 ) then
         D_near = Rocosl - sqrt( rootterm )
         D_far  = Rocosl + sqrt( rootterm )
      endif

      if ( gal_long.ge.90.d0  .and. gal_long.le.270.d0 ) then
         D_near = Rocosl + sqrt( rootterm )
         D_far  = D_near
      endif

      if ( gal_long.gt.270.d0 .and. gal_long.lt.360.d0 ) then
         D_near = Rocosl - sqrt( rootterm )
         D_far  = Rocosl + sqrt( rootterm )
      endif

      return
      end

c     ========================================================
      subroutine radec_to_galactic ( ra, dec, ell, bee )

c     Converts (ra,dec) in J2000 coordinates to Galactic (ell,bee)

c     inputs:    ra  (decimal hours)  in J2000 !
c                dec (decimal degrees)

c     outputs:   ell (galactic longitude in degrees)
c                bee (galactic latitude in degrees)

      implicit real*8 (a-h,o-z)

      pi        = 4.0d0*atan(1.0d0)
      hr_rad    = pi/12.d0      ! convert hours to radians
      deg_rad   = pi/180.d0     ! convert degrees to radians
      hr_deg    = 15.d0         ! convert hours to degrees

c     B1950 Galactic plane-pole angles (IAU standard)
c      g_alpha  = 282.25d0      
c      g_dec    = 62.6d0        
c      g_33     = 33.0d0        

c     J2000 Galactic plane-pole angles (modified)
      g_pole_ra = 12.d0 + 51.d0/60.d0 + 26.2817d0/3600.d0     ! hrs
      g_alpha   = (g_pole_ra + 6.d0) * hr_deg                 ! deg
                       ! = 282.859507
      g_pole_dec= 27.d0 + 07.d0/60.d0 + 42.013d0/3600.d0      ! deg
      g_dec     = 90.d0 - g_pole_dec                          ! deg
                       ! = 62.8716631
      g_33      = 33.d0 - 0.068d0 ! Longitude of ascending node of gal plane
                       ! = 32.932
                       ! Determined g_33 by requiring that the
                       ! B1950 defined G.C., when precessed to J2000,
                       ! still gives l=0.0 and b=0.0. (works to <0.001 deg)

c     Useful quantities...
      cos_a     = cos( ra*hr_rad - g_alpha*deg_rad )
      sin_a     = sin( ra*hr_rad - g_alpha*deg_rad )
      cos_d     = cos( dec*deg_rad )
      sin_d     = sin( dec*deg_rad )
      cos_gd    = cos( g_dec*deg_rad )
      sin_gd    = sin( g_dec*deg_rad )

c     Now calculate central galactic coordinates given ra,dec...
      sin_b   = sin_d*cos_gd - cos_d*sin_a*sin_gd
      bee     = asin( sin_b )
      cos_b   = cos( bee )

      cos_ell = cos_d*cos_a / cos_b
      sin_ell = ( cos_d*sin_a*cos_gd + sin_d*sin_gd ) / cos_b
      ell     = atan2( sin_ell,cos_ell )

      ell     = ell/deg_rad + g_33
      bee     = bee/deg_rad

      if ( ell .lt. 0.d0 ) ell = ell + 360.d0

      return
      end

c     ========================================================
      subroutine standard_vlsr_to_helio ( ra_rad, dec_rad, V_proj )

c     Enter (RA,Dec) = (ra_rad,dec_rad)  in radians

c     Gives V_proj   =  V(Heliocentric) - V(LSR)  in km/s
c         ie, add V_proj to V(LSR) to get V(Heliocentric)

c     Uses Standard Solar Motion (old values)

      implicit real*8 (a-h,o-z)

      pi = 4.d0 * atan( 1.d0 )
      hr_rad  = pi / 12.d0
      deg_rad = pi / 180.d0

c     LSR defined by removing peculiar Solar Motion of
c     20.0 km/s toward 18.0 hours, +30.0 degrees (RA,Dec)
c     (Actually defined in 1900.0 system, technically should
c     precess to J2000)

      Vsun     = 20.d0                  ! km/s
c     RA_Vsun  = 18.d0                  ! hours      1900.0
c     Dec_Vsun = 30.d0                  ! degrees
      RA_Vsun  = 18.0640d0              ! hours      J2000.0
      Dec_Vsun = 30.0024d0              ! degrees

      cos_RA = cos( RA_Vsun * hr_rad )
      sin_RA = sin( RA_Vsun * hr_rad )

      cos_Dec= cos( Dec_Vsun* deg_rad )
      sin_Dec= sin( Dec_Vsun* deg_rad )

c     In equatorial Cartesian frame the Solar Motion is
      Xo = Vsun * cos_RA * cos_Dec       ! toward RA=0h in equat plane
      Yo = Vsun * sin_RA * cos_Dec       ! toward RA=6h in equat plane
      Zo = Vsun * sin_Dec                ! toward Dec=90 (North)

c     Projection of Solar Motion toward a source at RA,Dec
      cos_alpha = cos( ra_rad )
      sin_alpha = sin( ra_rad )

      cos_delta = cos( dec_rad )
      sin_delta = sin( dec_rad )

      V_proj = -Xo*cos_alpha*cos_delta - Yo*sin_alpha*cos_delta -
     +          Zo*sin_delta

      return
      end

c     ========================================================
      subroutine hmsrad ( hms, hms_rad )

C     Converts hhmmss.sss format to radians

      implicit real*8 (a-h,o-z)

      pi = 4.d0 * atan(1.d0)

      xhms = abs(hms)

      ih = xhms/10000.d0 + 1.d-10
      im = dmod(xhms,10000.d0)/100.d0 + 1.d-10
      s  = dmod(xhms,100.d0)

      hms_rad = dfloat(ih) + dfloat(im)/60.d0 + s/3600.d0
      hms_rad = hms_rad * pi / 12.d0
      if ( hms .lt. 0.d0 )  hms_rad = -hms_rad

      return
      end

c     ========================================================
      subroutine dmsrad ( dms, dms_rad )

c     Converts dd''""."" format to radians

      implicit real*8 (a-h,o-z)

      pi = 4.d0 * atan(1.d0)
      deg_rad = pi / 180.d0

      xdms = abs(dms)

      id = xdms/10000.d0
      im = dmod(xdms,10000.d0)/100.d0
      s  = dmod(xdms,100.d0)

      dms_rad = dfloat(id) + dfloat(im)/60.d0 + s/3600.d0
      dms_rad = dms_rad * deg_rad
      if ( dms .lt. 0.d0 )  dms_rad = -dms_rad

      return
      end

