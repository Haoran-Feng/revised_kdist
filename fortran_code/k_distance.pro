      pro k_distance,input,coord1,coord2,vlsr,nearfars,c_code,dk_full,output=output,parfile=parfile,silent=silent

;c     Calculates kinematic distances using the revised prescription given in 
;c     Reid et al 2009 (Paper VI):  
;c     Version 2 (projects Vlsr to Galactic Plane)
;c     Version 3 (fixes bug: used Theta instead of Theta_0 in rootterm)
;      tranlated into IDL by Z. Jiang 2019.06.19

      pi = !dpi
      deg_rad = pi/180.d0

;c     =============================================================
;c     Get Galaxy, Solar Motion, and average source parameters from ascii 
;c     text file named "parameter_file.inp".  See example file for details

      if not keyword_set(parfile) then file_name = '/Users/zbjiang/IDL/fortran/parameter_file.inp'  else file_name=parfile
      galaxy_parameters, file_name, $
             Ro, To, dTdr, Uo, Vo, Wo, Us, Vs, Ws, err_v_lsr 
      if keyword_set(output) then openw,lun,output,/get_lun else lun=-1
      if ~keyword_set(silent) then begin
      print, $
        '! Program: revised_kineamtic_distance:', $
           ' Version 3; (M. Reid; November 2011)'
      print,'! Assumed Galactic paramters:  Ro, To, dTdr,  '   
      print, Ro, To, dTdr ,format='(3f7.2)'
        
      print,'! Assumed Solar Motion parms:  Uo, Vo, Wo  ,'
      print, Uo, Vo, Wo, format='(3f7.2)'
      print,'! Assumed Source Motion parms: Us, Vs, Ws  ,'
      print, Us, Vs, Ws, format='(3f7.2)'
      print,'! Assumed uncertainties in V_lsr          '
      print, err_v_lsr, format='(f7.2)'
      
      printf,lun, '! Source     Gal Long  Gal Lat    V_lsr     V_rev  ', $
            '  Rev. D_k     +/-  '
      printf, lun,     '!              (deg)    (deg)    (km/s)    (km/s)  ', $
            '   (kpc)      (kpc))'
      endif
;c     ============================================================
;c     Get input source data from ascii text file "source_file.dat"
;c     See example file for details
    if n_params() lt 1 then $
      file_name = 'source_file.txt' else file_name=input
;c        Input data for a source ...
;         indicating the coordinate code 1:Galactic, 2:RA_DEC in xx:xx:xx.x style, 3: RA_DEC in degrees
    if n_params() le 1 then readcol,file_name,srcs, coord1, coord2, v_lsr_s, nearfars,c_code,$
        format='AAAFII',comment='!',silent=silent $
        else begin
        srcs=input  & v_lsr_s=vlsr 
        endelse
        if n_elements(srcs) lt 2 then srcs=replicate(srcs,n_elements(coord1))
;     if more than one sources, begin a loop

            for i=0,n_elements(srcs)-1 do begin

;c           Convert coordinates to radians and then to hours, deg if input coordinates in RA_DEC
               src=srcs[i]
               farnear=nearfars[i]
               v_lsr=v_lsr_s[i]
             case c_code[i] of 
             1: begin
                ell=coord1[i]
                bee=coord2[i]
                end
             2: begin
                ra=hex2dec(coord1[i], /hour)
                dec=hex2dec(coord2[i])
                end
             3: begin
                ra=coord1[i]
                dec=coord2[i]
                end
             else:             
             endcase   
             glactc,ra,dec,2000,ell,bee,(c_code[i] gt 1)?1:2,/degree        


;               if c_code[i] gt 1 then begin
;                if c_code[i] eq 2 then begin
;                ra=hex2dec(coord1[i], /hour)
;                dec=hex2dec(coord2[i])
;                    print,'c_code',c_code[i],2
;                endif else begin
;                    ra=coord1[i]
;                    dec=coord2[i]
;                    print,'c_code',c_code[i],3
;                    endelse
;                glactc,ra,dec,2000.,ell,bee,1
;                endif else begin
;                    ell=coord1[i] & bee=coord2[i]
;                    glactc,ra,dec,2000.,ell,bee,2
;                    print,'c_code',c_code[i],1
;                    endelse
;                    
;                    
;c           Get Galactic coordinates of source
;            call radec_to_galactic ( ra, dec,  ell, bee )

;c           Calculate "adjusted" Vlsr, including effects of (Us,Vs).
;c           Since applying (Us,Vs) needs Dk, so do both at same time.
             calc_Dk, $
                  Ro, To, dTdr, Uo, Vo, Wo, Us, Vs, Ws, $
                  src, ra, dec, farnear, v_lsr, $
                  v_lsr_rev, Dk_best 

;c           -------------------------------------------------------
;c           Estimation of (asymmetric) uncertainties based only on
;c           the uncertainty in Vlsr supplied in "parameter_file.inp".
;c           Assumes Galactic & Source parameters have no uncertainties.
;
;c           Check if v_lsr exceeds tangent point speed.  If so,
;c           use tangent point speed as center for velocity offsets
;c           used for error estimate.
;
;c           Calculate tangent point speed (T_tp)
            sin_l = sin( ell*deg_rad )
            R_tp = Ro * sin_l                   ;! kpc
            T_tp = To + dTdR * (R_tp - Ro)      ;! km/s

;c           Check of v_lsr_rev exceeds tangent point speed
            del_v = 0.d0
            v_tp  = -9999.9d0
            if ( ell lt 90.d0  ) then begin         ;! quadrant 1
               v_tp = +T_tp - To*sin_l          ;! km/s 
               if ( v_lsr_rev gt v_tp )  then del_v = v_lsr_rev - v_tp
            endif
            if ( ell gt 270.d0 ) then begin        ;! quadrant 4
               v_tp = -T_tp - To*sin_l          ;! km/s 
               if ( v_lsr_rev lt v_tp ) then del_v = v_lsr_rev - v_tp
            endif
;c           Adjust v_lsr (if needed) so as not to exceed tangent point
            v_lsr_adj = v_lsr - del_v

;c           Get lower velocity distance uncertainty
            v_lsr_low = v_lsr_adj - err_v_lsr 
            calc_Dk , $
                  Ro, To, dTdr, Uo, Vo, Wo, Us, Vs, Ws, $
                  src, ra, dec, farnear, v_lsr_low,$
                  v_lsr_rev_low, D1 
;c           Flag if close to tangent point velocity
            if ( abs(v_lsr_rev_low - v_tp) lt err_v_lsr ) then D1=Dk_best

;c           Get higher velocity distance uncertainty
            v_lsr_high = v_lsr_adj + err_v_lsr 
            calc_Dk, $
                  Ro, To, dTdr, Uo, Vo, Wo, Us, Vs, Ws, $
                  src, ra, dec, farnear, v_lsr_high, $
                  v_lsr_rev_high, D2 
;c           Flag if close to tangent point velocity
            if ( abs(v_lsr_rev_high - v_tp) lt err_v_lsr ) then D2=Dk_best

;c           Calculate + (high) and - (low) errors
            if ( D2 gt D1 ) then  begin
                 d_err_low  = D1  - Dk_best
                 d_err_high = D2  - Dk_best
               endif else begin
                 d_err_low  = D2  - Dk_best
                 d_err_high = D1  - Dk_best
            endelse

;c           If flagged distance (error=0), use other error estimate
            if ( abs(d_err_low) eq 0.d0 )  then d_err_low  = -d_err_high
            if ( d_err_high     eq 0.d0 )  then d_err_high = -d_err_low

;c           Check for pathalogical cases...
            if ( Dk_best le 0.0 ) then begin
               Dk_best    = 0.0
               d_err_low  = 0.0
               d_err_high = 0.0
            endif
            if ( Dk_best+d_err_low lt 0.0 )  then d_err_low = -Dk_best

            if ~keyword_set(silent) then printf,lun, src, ell, bee, v_lsr, v_lsr_rev, $
                            Dk_best, d_err_high, d_err_low, $
            format='(1x,a12,2f8.3,2f10.1,f10.2,2f7.2)'

        endfor
        if keyword_set(output) then free_lun,lun
        dk_full=[dk_best,d_err_high,d_err_low]
      end


;c===================================================================
      pro galaxy_parameters, file_name, $
                Ro, To, dTdr, Uo, Vo, Wo, Us, Vs, Ws, err_v_lsr

        readcol,file_name,xx,yy,format='F,F',comment='!',/silent
      Ro=xx[0]
      To=xx[1]
      dTdr=xx[2]
      Uo=xx[3]
      Vo=xx[4]
      Wo=xx[5]
      Us=xx[6]
      Vs=xx[7]
      Ws=xx[8]
      err_v_lsr=xx[9]

      return
      end

;c===================================================================
      pro calc_Dk, $
             Ro, To, dTdr, Uo, Vo, Wo, Us, Vs, Ws,  $
             src, ra, dec, farnear, v_lsr, $
             v_lsr_rev, Dk 

;c     Calculate revised Vlsr by converting standard Vlsr back to
;c     heliocentric, apply modern Solar Motion values (Uo,Vo,Wo), and 
;c     remove effects of average source non-ciruclar motion (Us,Vs,Ws).  
;c     Then calculate kinematic distance using the linear rotation 
;c     curve specified by Ro, To, and dTdR

      pi     = !dpi
      deg_rad = pi/180.d0

      n_iter_max = 1000

;c     ============================================================
;c     Get source galactic coordinates...

      glactc,ra, dec, 2000,gal_long, gal_lat,1,/degree

      gal_long_rad = gal_long * deg_rad       ; radians
      cos_l = cos( gal_long_rad )
      sin_l = sin( gal_long_rad )

      gal_lat_rad  = gal_lat * deg_rad        ;! radians
      cos_b = cos( gal_lat_rad )
      sin_b = sin( gal_lat_rad )

;c     -------------------------------------------------------------
;c     Convert to true Heliocentric frame (v_rad)
;c     Add back Sun's peculiar motion to radial velocity
;c     Use old Standard Solar Motion (Sun moves at 20 km/s toward
;c     18h, +30deg in 1900 coordinates), since this has been used to 
;c     define V_lsr:

      Uo_IAU = 10.27d0            ;! km/s precessed to J2000
      Vo_IAU = 15.32d0
      Wo_IAU =  7.74d0

      v_helio = v_lsr - (Vo_IAU*sin_l + Uo_IAU*cos_l)*cos_b $
                     -  Wo_IAU*sin_b

;c     -------------------------------------------------------------
;c     Make "new" V(LSR) using best Solar Motion
;c      eg, Hipparcos (Dehnen & Binney 1997) gives 
;c      Uo = 10.00d0                ;! km/s
;c      Vo =  5.25d0
;c      Wo =  7.17d0
;  
      v_newlsr = v_helio + (Vo*sin_l + Uo*cos_l)*cos_b $
                        +  Wo*sin_b

;c     --------------------------------------------------------
;c     Remove effects of common peculiar motions specified in
;c     "parameter_file.inp"

;c     If dTdr.ne.0, need to know distance to get source Galactocentric 
;c     radius (Rs) to evalute rotation curve.   So must iterate...
      n_iter =  0
      del_d  = 99.d0
      Dk     =  3.d0

      while ( del_d gt 0.01d0  and n_iter lt n_iter_max ) do begin

;c        Save old value  of kinematic distance
         Dk_old = Dk

;c        Calculate "gamma" angle and projected Galactocentric radius 
         d_proj = Dk * cos_b                   ;  ! kpc in Gal Plane
         r_sq   = Ro^2 + d_proj^2 - 2.d0*Ro * d_proj * cos_l
         r_proj = sqrt( r_sq )                  ; ! kpc in Gal Plane

;c        Calculate Galactocentric longitude (beta in paper)...
         sin_beta =   d_proj  * sin_l     / r_proj
         cos_beta = ( Ro - d_proj*cos_l ) / r_proj
         beta     = atan( sin_beta, cos_beta )  ;! radians
         beta_deg = beta / deg_rad               ;! deg

;c        Calculate Sun-Maser-GC angle...
         gamma = pi - gal_long_rad - beta        ;! radians
         cos_gamma = cos( gamma )
         sin_gamma = sin( gamma )

         v_fixed = v_newlsr - (Vs*sin_gamma - Us*cos_gamma)*cos_b $
                          -  Ws*sin_b          ;! km/s

;c        -----------------------------------------------------------------
;c        Calculate a kinematic distance using best Ro, To and dTdr

         Rs = r_proj
         V_proj = v_fixed * cos_b
         kinematic_distance, v_proj, gal_long, Ro, To, $
                                  r_proj, dTdr, $
                                  D_near, D_far 

         Dk = D_near
         if ( farnear ne 0.d0 )  then Dk = D_far
;print,d_near,d_far
;c        Ignore "farnear" flag if one of the values is zero
         if ( D_near le 0.d0) and (D_far  gt 0.d0 )  then Dk = D_far
         if ( D_far  le 0.d0) and (D_near gt 0.d0 ) then Dk = D_near

         del_d = abs( Dk - Dk_old )
         n_iter = n_iter + 1

      endwhile

      v_lsr_rev = v_fixed

      return
      end

;c     =======================================================
      pro kinematic_distance, Vlsr, gal_long, Ro, To, $
                                     Rs, dTdr, $
                                     D_near, D_far 

;c     Caluclate kinematic distance given Vlsr (projected in Gal plane)
;c     and information required to construct the kinematic model.  
;c     Returns both near and far distances (when appropriate)

;      implicit real*8 (a-h,o-z)

      pi     = 4.0d0*atan(1.0d0)	
      deg_rad = pi/180.d0	     ;! convert degrees to radians

      glongrad = gal_long*deg_rad 

      cos_l = cos(glongrad)
      sin_l = sin(glongrad)

      Rosinl = Ro * sin_l
      Rocosl = Ro * cos_l

;c     Non-flat rotation curve allowed...
      Theta = To + dTdr*(Rs - Ro)
      Tosinl= To    * sin_l
      Tsinl = Theta * sin_l

      rootterm = Rocosl^2 + ( Tsinl / ( Tosinl/Ro + Vlsr/Ro ) )^2  $
              - Ro^2
      if ( rootterm lt 0 ) then rootterm = 0.d0

      if ( gal_long ge 0.d0)   and (gal_long lt 90.d0 ) then begin
         D_near = Rocosl - sqrt( rootterm )
         D_far  = Rocosl + sqrt( rootterm )
      endif

      if ( gal_long ge 90.d0)  and (gal_long le 270.d0 ) then begin
         D_near = Rocosl + sqrt( rootterm )
         D_far  = D_near
      endif

      if ( gal_long gt 270.d0) and (gal_long lt 360.d0 ) then begin
         D_near = Rocosl - sqrt( rootterm )
         D_far  = Rocosl + sqrt( rootterm )
      endif

      return
      end


;c     ========================================================
      pro standard_vlsr_to_helio, ra_rad, dec_rad, V_proj 

;c     Enter (RA,Dec) = (ra_rad,dec_rad)  in radians
;
;c     Gives V_proj   =  V(Heliocentric) - V(LSR)  in km/s
;c         ie, add V_proj to V(LSR) to get V(Heliocentric)
;
;c     Uses Standard Solar Motion (old values)
;
;      implicit real*8 (a-h,o-z)

      pi = 4.d0 * atan( 1.d0 )
      hr_rad  = pi / 12.d0
      deg_rad = pi / 180.d0

;c     LSR defined by removing peculiar Solar Motion of
;c     20.0 km/s toward 18.0 hours, +30.0 degrees (RA,Dec)
;c     (Actually defined in 1900.0 system, technically should
;c     precess to J2000)

      Vsun     = 20.d0                  ;! km/s
;c     RA_Vsun  = 18.d0                  ! hours      1900.0
;c     Dec_Vsun = 30.d0                  ! degrees
      RA_Vsun  = 18.0640d0              ;! hours      J2000.0
      Dec_Vsun = 30.0024d0              ;! degrees

      cos_RA = cos( RA_Vsun * hr_rad )
      sin_RA = sin( RA_Vsun * hr_rad )

      cos_Dec= cos( Dec_Vsun* deg_rad )
      sin_Dec= sin( Dec_Vsun* deg_rad )

;c     In equatorial Cartesian frame the Solar Motion is
      Xo = Vsun * cos_RA * cos_Dec       ;! toward RA=0h in equat plane
      Yo = Vsun * sin_RA * cos_Dec       ;! toward RA=6h in equat plane
      Zo = Vsun * sin_Dec                ;! toward Dec=90 (North)

;c     Projection of Solar Motion toward a source at RA,Dec
      cos_alpha = cos( ra_rad )
      sin_alpha = sin( ra_rad )

      cos_delta = cos( dec_rad )
      sin_delta = sin( dec_rad )

      V_proj = -Xo*cos_alpha*cos_delta - Yo*sin_alpha*cos_delta - $
               Zo*sin_delta

      return
      end


