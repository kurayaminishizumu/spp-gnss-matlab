function tropo_delay_m = saastamoinen_tropo(h_user_m, el_sat_deg, P_mbar, T_kelvin, Hum_percent)
if isempty(P_mbar) 
    P0 = 1013.25;
    P_mbar = P0 * (1 - 2.2557e-5 * h_user_m)^5.25588;
end
if isempty(T_kelvin) 
    T0 = 288.15; 
    L = 0.0065;
    T_kelvin = T0 - L * h_user_m;
    if T_kelvin < 200, T_kelvin = 200; end 
end
if isempty(Hum_percent)
    Hum_percent = 50; 
end

e_sat_Pa = 6.108 * exp((17.15 * (T_kelvin-273.15)) / (234.7 + (T_kelvin-273.15))); 
e_Pa = (Hum_percent/100) * e_sat_Pa; 
e_mbar = e_Pa / 100; 

el_rad = el_sat_deg * pi/180;
zhd = (0.0022768 * P_mbar) / (1 - 0.00266 * cos(2*0) - 0.00028 * h_user_m/1e3);
zwd = (0.002277 * (1255/T_kelvin + 0.05) * e_mbar);

if el_rad == 0 
    map_func = 100; 
else
    map_func = 1.001 / sqrt(0.002001 + sin(el_rad)^2); 
end

tropo_delay_m = (zhd + zwd) * map_func;
if tropo_delay_m < 0, tropo_delay_m = 0; end
end