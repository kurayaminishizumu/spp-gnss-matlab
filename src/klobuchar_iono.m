function iono_delay_m = klobuchar_iono(gps_sow, lat_user_deg, lon_user_deg, az_sat_deg, el_sat_deg, alpha, beta)
C_LIGHT = 299792458.0;
DEG2RAD = pi/180;
RAD2SC = 1 / (2*pi);
SC2RAD = pi; 

el_rad = el_sat_deg * DEG2RAD;
az_rad = az_sat_deg * DEG2RAD;
lat_user_rad = lat_user_deg * DEG2RAD;
lon_user_rad = lon_user_deg * DEG2RAD;

psi = 0.0137 / (el_rad/SC2RAD + 0.11) - 0.022; 

phi_i = lat_user_rad/SC2RAD + psi * cos(az_rad); 
if phi_i > 0.416, phi_i = 0.416;
elseif phi_i < -0.416, phi_i = -0.416;
end
phi_i_rad = phi_i * SC2RAD; 

lambda_i = lon_user_rad/SC2RAD + (psi * sin(az_rad) / cos(phi_i_rad)); 
lambda_i_rad = lambda_i * SC2RAD; 

phi_m = phi_i + 0.064 * cos(lambda_i_rad - 1.617); 
phi_m_rad = phi_m * SC2RAD; 

t_local = mod(4.32e4 * lambda_i + gps_sow, 86400); 
if t_local >= 86400, t_local = t_local - 86400;
elseif t_local < 0, t_local = t_local + 86400;
end

F = 1.0 + 16.0 * (0.53 - el_rad/SC2RAD)^3; 

PER = beta(1) + beta(2)*phi_m + beta(3)*phi_m^2 + beta(4)*phi_m^3;
if PER < 72000, PER = 72000; end

AMP = alpha(1) + alpha(2)*phi_m + alpha(3)*phi_m^2 + alpha(4)*phi_m^3; 
if AMP < 0, AMP = 0; end

x = 2 * pi * (t_local - 50400) / PER;

if abs(x) < 1.57
    iono_delay_sec = F * (5e-9 + AMP * (1 - x^2/2 + x^4/24));
else
    iono_delay_sec = F * 5e-9; 
end

iono_delay_m = C_LIGHT * iono_delay_sec;
end