function dops = calculate_dops(Q_matrix)

qxx = Q_matrix(1,1);
qyy = Q_matrix(2,2);
qzz = Q_matrix(3,3);
qtt = Q_matrix(4,4); 

gdop = sqrt(qxx + qyy + qzz + qtt);
pdop = sqrt(qxx + qyy + qzz);
hdop = sqrt(qxx + qyy); 
vdop = sqrt(qzz);
tdop = sqrt(qtt);

dops = [gdop, pdop, hdop, vdop, tdop];
end