function [sat_pos_xyz, sat_clk_err_sec, rel_corr_m] = calculate_satellite_position_gps(eph, t_sv_gpst, gps_week_obs, P_obs_for_TGD_correction)


MU_GPS = 3.986005e14;     
OMEGA_EARTH = 7.2921151467e-5;

sat_pos_xyz = [NaN, NaN, NaN];
sat_clk_err_sec = NaN;
rel_corr_m = NaN;

try
    A = eph.sqrt_A^2;
    n0 = sqrt(MU_GPS / A^3); 

    tk = t_sv_gpst - eph.Toe;

    if gps_week_obs > eph.GPSWeek
        if t_sv_gpst < eph.Toe
            tk = tk + 604800; 
        end

    elseif gps_week_obs < eph.GPSWeek
         if t_sv_gpst > eph.Toe 
            tk = tk - 604800; 
         end
    end

    if tk > 302400
        tk = tk - 604800;
    elseif tk < -302400
        tk = tk + 604800;
    end

    n = n0 + eph.Delta_n;
    Mk = eph.M0 + n * tk;

    Ek = Mk; 
    for i = 1:10 
        Ek_old = Ek;
        Ek = Mk + eph.e * sin(Ek_old);
        if abs(Ek - Ek_old) < 1e-12, break; end
    end

    sin_nuk = sqrt(1 - eph.e^2) * sin(Ek) / (1 - eph.e * cos(Ek));
    cos_nuk = (cos(Ek) - eph.e) / (1 - eph.e * cos(Ek));
    nuk = atan2(sin_nuk, cos_nuk);

    phi_k = nuk + eph.omega;

    delta_uk = eph.Cus * sin(2*phi_k) + eph.Cuc * cos(2*phi_k); 
    delta_rk = eph.Crs * sin(2*phi_k) + eph.Crc * cos(2*phi_k);
    delta_ik = eph.Cis * sin(2*phi_k) + eph.Cic * cos(2*phi_k);

    uk = phi_k + delta_uk;
    rk = A * (1 - eph.e * cos(Ek)) + delta_rk; 
    ik = eph.i0 + delta_ik + eph.IDOT * tk; 

    xk_prime = rk * cos(uk);
    yk_prime = rk * sin(uk);

    Omega_k = eph.OMEGA0 + (eph.OMEGA_DOT - OMEGA_EARTH) * tk - OMEGA_EARTH * eph.Toe;

    sat_pos_xyz(1) = xk_prime * cos(Omega_k) - yk_prime * cos(ik) * sin(Omega_k);
    sat_pos_xyz(2) = xk_prime * sin(Omega_k) + yk_prime * cos(ik) * cos(Omega_k);
    sat_pos_xyz(3) = yk_prime * sin(ik);

    dt_clock = t_sv_gpst - eph.toc_second;
    toc_dt_obj = datetime(eph.toc_year, eph.toc_month, eph.toc_day, eph.toc_hour, eph.toc_minute, eph.toc_second);
    [~, ~, toc_sow_eph_week] = date_to_gps_week_sow(toc_dt_obj);
    
    dt_clock = t_sv_gpst - toc_sow_eph_week;
    if eph.GPSWeek > gps_week_obs
        dt_clock = dt_clock + (eph.GPSWeek - gps_week_obs) * 604800;
    elseif eph.GPSWeek < gps_week_obs 
        dt_clock = dt_clock - (gps_week_obs - eph.GPSWeek) * 604800;
    end

    if dt_clock > 302400
        dt_clock = dt_clock - 604800;
    elseif dt_clock < -302400
        dt_clock = dt_clock + 604800;
    end
    
    sat_clk_bias_sec = eph.af0 + eph.af1 * dt_clock + eph.af2 * dt_clock^2;

    F = -4.442807633e-10; 
    rel_corr_sec = F * eph.e * eph.sqrt_A * sin(Ek);

    sat_clk_err_sec = sat_clk_bias_sec + rel_corr_sec - eph.TGD;
    rel_corr_m = rel_corr_sec * 299792458.0; 

catch ME
    fprintf('Error calculating sat pos/clk for PRN %d: %s\n', eph.PRN, ME.message);
    sat_pos_xyz = [NaN, NaN, NaN];
    sat_clk_err_sec = NaN;
    rel_corr_m = NaN;
    return;
end
end

function [gps_week, day_of_week, seconds_of_week] = date_to_gps_week_sow(dt_obj)
    gps_epoch_dt = datetime(1980, 1, 6, 0, 0, 0, 'TimeZone', 'UTC');
     if isempty(dt_obj.TimeZone)
        dt_obj.TimeZone = 'UTC';
    end
    time_diff_seconds = posixtime(dt_obj) - posixtime(gps_epoch_dt);
    gps_week = floor(time_diff_seconds / (7 * 24 * 3600));
    seconds_into_current_week = mod(time_diff_seconds, (7 * 24 * 3600));
    day_of_week = floor(seconds_into_current_week / (24 * 3600));
    seconds_of_week = seconds_into_current_week;
end