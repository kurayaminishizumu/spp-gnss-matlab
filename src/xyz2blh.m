function blh = xyz2blh(xyz, a, f)
X = xyz(1); Y = xyz(2); Z = xyz(3);

e_sq = 2*f - f^2; 
e_prime_sq = e_sq / (1 - e_sq);

p = sqrt(X^2 + Y^2);
theta = atan2(Z*a, p*(a*(1-e_sq)));

lat = atan2(Z + e_prime_sq * a * (1-e_sq)/a * sin(theta)^3, ...
            p - e_sq * a * cos(theta)^3);

lon = atan2(Y, X);

N = a / sqrt(1 - e_sq * sin(lat)^2);
h = p / cos(lat) - N;

if abs(cos(lat)) < 1e-9 
    if Z > 0
        h = Z - N*(1-e_sq);
    else
        h = -Z - N*(1-e_sq);
    end
end

blh = [lat; lon; h];
end