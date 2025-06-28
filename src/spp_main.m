clear; clc; close all;

C_L = 299792458.0; 
OMEGA_E = 7.2921151467e-5; 
WGS84_A = 6378137.0; 
WGS84_F = 1/298.257223563;
MU_GPS = 3.986005e14; 

obs_file = fullfile('..', 'data', 'abmf2740.23o');
nav_file = fullfile('..', 'data', 'brdc2740.23p');

approx_pos_xyz = [0, 0, 0];
elevation_mask = 10; 
max_iter = 10; 
iter_threshold = 1e-4;
use_iono_corr = true;  
use_tropo_corr = true; 
 
% --- 1. 读取导航文件 ---
fprintf('正在读取导航文件: %s\n', nav_file);
[ephemeris, iono_params_gps] = read_rinex_nav_gps(nav_file);
if isempty(ephemeris)
    error('无法读取导航文件或未找到GPS星历。');
end
fprintf('导航文件读取完毕。找到 %d 条GPS星历。\n', length(fieldnames(ephemeris)));

% --- 2. 读取观测文件 ---
fprintf('正在读取观测文件: %s\n', obs_file);
[obs_data, approx_pos_rinex_header, ~] = read_rinex_obs_gps(obs_file); 
if isempty(obs_data)
    error('无法读取观测文件或未找到GPS观测数据。');
end
fprintf('观测文件读取完毕。找到 %d 个观测历元。\n', length(obs_data));

if ~all(approx_pos_xyz) && ~isempty(approx_pos_rinex_header) && all(isfinite(approx_pos_rinex_header)) && length(approx_pos_rinex_header)==3
    approx_pos_xyz = approx_pos_rinex_header'; 
    fprintf('使用RINEX文件头中的近似坐标: [%.2f, %.2f, %.2f]\n', approx_pos_xyz(1), approx_pos_xyz(2), approx_pos_xyz(3));
else
    if ~all(approx_pos_xyz) || ~all(isfinite(approx_pos_xyz)) || length(approx_pos_xyz)~=3
       approx_pos_xyz = [WGS84_A, 0, 0];
       fprintf('未提供有效近似坐标或RINEX头中无有效近似坐标，使用默认初值: [%.2f, %.2f, %.2f]\n', approx_pos_xyz(1), approx_pos_xyz(2), approx_pos_xyz(3));
    else
       fprintf('使用用户指定的近似坐标: [%.2f, %.2f, %.2f]\n', approx_pos_xyz(1), approx_pos_xyz(2), approx_pos_xyz(3));
    end
end
if size(approx_pos_xyz, 1) > 1 
    approx_pos_xyz = approx_pos_xyz';
end


% --- 3. 主处理循环 ---
num_epochs = length(obs_data);
results_ecef = zeros(num_epochs, 4);
results_blh = zeros(num_epochs, 3);  
results_dops = zeros(num_epochs, 5); 
processed_epochs = 0;

fprintf('开始处理历元...\n');
for epoch = 1:num_epochs
    current_obs = obs_data{epoch};
    if isempty(current_obs) || ~isfield(current_obs,'time')
        fprintf('历元 %d: 数据为空，跳过。\n', epoch);
        results_ecef(epoch, :) = NaN;
        results_blh(epoch, :) = NaN;
        results_dops(epoch, :) = NaN;
        continue;
    end
    t_obs_gpst = current_obs.time;
    week_obs = current_obs.week;  

    if ~isfield(current_obs, 'sv') || isempty(fieldnames(current_obs.sv))
        fprintf('历元 %d: 无卫星观测数据，跳过。\n', epoch);
        results_ecef(epoch, :) = NaN;
        results_blh(epoch, :) = NaN;
        results_dops(epoch, :) = NaN;
        continue;
    end
    available_sats = fieldnames(current_obs.sv);
    num_available_sats = length(available_sats);

    if num_available_sats < 4
        fprintf('历元 %d: 卫星数不足 (%d)，跳过。\n', epoch, num_available_sats);
        results_ecef(epoch, :) = NaN;
        results_blh(epoch, :) = NaN;
        results_dops(epoch, :) = NaN;
        continue;
    end

    rcv_pos_iter = approx_pos_xyz; 
    rcv_clk_corr_iter = 0;         
    A_final_iter = []; 

    for iter = 1:max_iter
        A = []; 
        L = [];
        valid_sat_count_iter = 0; 

        for i = 1:num_available_sats
            prn_str = available_sats{i};
            if ~startsWith(prn_str, 'G') || length(prn_str) < 2 || isnan(str2double(prn_str(2:end)))
                continue;
            end
            prn = str2double(prn_str(2:end));

            eph_field = sprintf('G%02d', prn);
            if ~isfield(ephemeris, eph_field)
                continue;
            end
            current_eph = ephemeris.(eph_field);

            if isfield(current_obs.sv.(prn_str), 'C1C') && ~isnan(current_obs.sv.(prn_str).C1C)
                P_obs = current_obs.sv.(prn_str).C1C;
            elseif isfield(current_obs.sv.(prn_str), 'P1') && ~isnan(current_obs.sv.(prn_str).P1)
                P_obs = current_obs.sv.(prn_str).P1;
            elseif isfield(current_obs.sv.(prn_str), 'C1P') && ~isnan(current_obs.sv.(prn_str).C1P)
                P_obs = current_obs.sv.(prn_str).C1P;
            else
                continue;
            end

            travel_time_approx = P_obs / C_L;
            t_transmit_gpst = t_obs_gpst - travel_time_approx;

            [sat_pos_xyz_tx, sat_clk_err_sec, ~] = calculate_satellite_position_gps(current_eph, t_transmit_gpst, week_obs, P_obs);
            if any(isnan(sat_pos_xyz_tx))
                continue;
            end

            geom_dist_current_iter = norm(sat_pos_xyz_tx - rcv_pos_iter);
            
            travel_time_precise = geom_dist_current_iter / C_L;
            t_transmit_gpst_precise = t_obs_gpst - travel_time_precise;
            [sat_pos_xyz_tx, sat_clk_err_sec, ~] = calculate_satellite_position_gps(current_eph, t_transmit_gpst_precise, week_obs, P_obs);
            if any(isnan(sat_pos_xyz_tx))
                continue;
            end
            sat_clk_corr_m = sat_clk_err_sec * C_L; 
            sagnac_corr_m = OMEGA_E * (sat_pos_xyz_tx(1)*rcv_pos_iter(2) - sat_pos_xyz_tx(2)*rcv_pos_iter(1)) / C_L; 
            
            geom_dist_final = norm(sat_pos_xyz_tx - rcv_pos_iter);

            [az, el, ~] = topocent(rcv_pos_iter, sat_pos_xyz_tx - rcv_pos_iter, WGS84_A, WGS84_F);
            el_deg = el * 180/pi;
            az_deg = az * 180/pi;

            if el_deg < elevation_mask
                continue;
            end

            iono_delay_m = 0;
            if use_iono_corr && ~isempty(iono_params_gps) && isfield(iono_params_gps, 'alpha') && isfield(iono_params_gps, 'beta')
                rcv_blh_approx_iter = xyz2blh(rcv_pos_iter, WGS84_A, WGS84_F);
                iono_delay_m = klobuchar_iono(t_obs_gpst, rcv_blh_approx_iter(1)*180/pi, rcv_blh_approx_iter(2)*180/pi, ...
                                              az_deg, el_deg, iono_params_gps.alpha, iono_params_gps.beta);
            end

            tropo_delay_m = 0;
            if use_tropo_corr
                rcv_blh_approx_iter = xyz2blh(rcv_pos_iter, WGS84_A, WGS84_F);
                tropo_delay_m = saastamoinen_tropo(rcv_blh_approx_iter(3), el_deg, [], [], []);
            end
            
            l_i = P_obs + sat_clk_corr_m - iono_delay_m - tropo_delay_m - sagnac_corr_m - geom_dist_final;

            unit_vector_x = -(sat_pos_xyz_tx(1) - rcv_pos_iter(1)) / geom_dist_final;
            unit_vector_y = -(sat_pos_xyz_tx(2) - rcv_pos_iter(2)) / geom_dist_final;
            unit_vector_z = -(sat_pos_xyz_tx(3) - rcv_pos_iter(3)) / geom_dist_final;

            A_row = [unit_vector_x, unit_vector_y, unit_vector_z, 1];
            A = [A; A_row];
            L = [L; l_i];
            valid_sat_count_iter = valid_sat_count_iter + 1;
        end 

        if valid_sat_count_iter < 4
            if iter == 1 && epoch <= 10
                fprintf('历元 %d: 迭代 %d, 可用卫星数 %d < 4, 跳过历元。\n', epoch, iter, valid_sat_count_iter);
            end
            A_final_iter = [];
            break; 
        end

        try
            delta_x = (A'*A)\(A'*L); 
        catch ME_ls
            fprintf('历元 %d: 最小二乘解算失败 (可能A''A奇异): %s\n', epoch, ME_ls.message);
            A_final_iter = [];
            break; % 跳出迭代
        end


        rcv_pos_iter = rcv_pos_iter + delta_x(1:3)'; 
        rcv_clk_corr_iter = rcv_clk_corr_iter + delta_x(4); 
        A_final_iter = A; 

        if norm(delta_x(1:3)) < iter_threshold
            break;
        end
        if iter == max_iter
             if epoch <= 10 
                 fprintf('历元 %d: 达到最大迭代次数 %d 未收敛。\n', epoch, max_iter);
             end
             A_final_iter = []; 
        end
    end 
    if ~isempty(A_final_iter) && valid_sat_count_iter >=4 
        results_ecef(epoch, 1:3) = rcv_pos_iter;
        results_ecef(epoch, 4) = rcv_clk_corr_iter;

        blh = xyz2blh(rcv_pos_iter, WGS84_A, WGS84_F);
        results_blh(epoch, :) = [blh(1)*180/pi, blh(2)*180/pi, blh(3)];
        
        try
            Q = inv(A_final_iter'*A_final_iter);
            dops_vals = calculate_dops(Q);
            results_dops(epoch,:) = dops_vals;
        catch ME_inv
            fprintf('历元 %d: 计算DOP时矩阵求逆失败: %s\n', epoch, ME_inv.message);
            results_dops(epoch,:) = NaN;
            dops_vals = [NaN, NaN, NaN, NaN, NaN]; 
        end


        processed_epochs = processed_epochs + 1;
        if mod(processed_epochs, 100) == 0 || processed_epochs == 1
           fprintf('历元 %d/%d: X=%.2f Y=%.2f Z=%.2f c*dt_r=%.2f Nsat=%d GDOP=%.1f\n', ...
                   epoch, num_epochs, rcv_pos_iter(1), rcv_pos_iter(2), rcv_pos_iter(3), rcv_clk_corr_iter, valid_sat_count_iter, dops_vals(1));
        end
        approx_pos_xyz = rcv_pos_iter; 
    else
        results_ecef(epoch, :) = NaN;
        results_blh(epoch, :) = NaN;
        results_dops(epoch, :) = NaN;
    end
end 

fprintf('处理完毕。\n');

valid_indices = ~isnan(results_blh(:,1)) & (results_blh(:,1) ~= 0 | results_blh(:,2) ~= 0 | results_blh(:,3) ~= 0) ; % 过滤掉NaN和全零结果
if ~any(valid_indices)
    disp('没有成功解算的有效历元，无法绘图。');
    return;
end

figure;
subplot(3,1,1);
plot(find(valid_indices), results_blh(valid_indices,1)); ylabel('纬度 (deg)'); title('定位结果 (大地坐标)'); grid on;
subplot(3,1,2);
plot(find(valid_indices), results_blh(valid_indices,2)); ylabel('经度 (deg)'); grid on;
subplot(3,1,3);
plot(find(valid_indices), results_blh(valid_indices,3)); ylabel('高程 (m)'); xlabel('成功解算的历元序号'); grid on;

figure;
plot(find(valid_indices), results_dops(valid_indices,1), 'DisplayName', 'GDOP'); hold on;
plot(find(valid_indices), results_dops(valid_indices,2), 'DisplayName', 'PDOP');
plot(find(valid_indices), results_dops(valid_indices,3), 'DisplayName', 'HDOP');
plot(find(valid_indices), results_dops(valid_indices,4), 'DisplayName', 'VDOP');
plot(find(valid_indices), results_dops(valid_indices,5), 'DisplayName', 'TDOP');
legend show; title('DOP 值'); xlabel('成功解算的历元序号'); ylabel('DOP'); grid on; 

max_gdop_plot = max(results_dops(valid_indices,1), [], 'omitnan');
if isfinite(max_gdop_plot) && max_gdop_plot > 0
    ylim([0, min(max_gdop_plot * 1.1, 50)]); 
else
    ylim([0,50]);
end


% --- 辅助函数: topocent (计算方位角、高度角) ---
function [az, el, d] = topocent(obs_pos_xyz, sat_vec_xyz_rel, const_A_in, const_F_in) 
    
    blh_obs_rad = xyz2blh(obs_pos_xyz, const_A_in, const_F_in);
    lat_rad = blh_obs_rad(1);
    lon_rad = blh_obs_rad(2);

    R_enu = [-sin(lon_rad),             cos(lon_rad),            0;
             -sin(lat_rad)*cos(lon_rad), -sin(lat_rad)*sin(lon_rad), cos(lat_rad);
              cos(lat_rad)*cos(lon_rad),  cos(lat_rad)*sin(lon_rad), sin(lat_rad)];

    if size(sat_vec_xyz_rel,1) == 1 
        sat_vec_xyz_rel_col = sat_vec_xyz_rel';
    else 
        sat_vec_xyz_rel_col = sat_vec_xyz_rel;
    end
    enu_vec = R_enu * sat_vec_xyz_rel_col; 
    
    E = enu_vec(1); 
    N = enu_vec(2); 
    U = enu_vec(3);

    d_sq_horiz = E^2 + N^2; 
    d = sqrt(d_sq_horiz + U^2); 
    
    if d < 1e-6 
        el = pi/2; 
        az = 0;    
    else
        el = atan2(U, sqrt(d_sq_horiz)); 
        az = atan2(E, N);          
    end
    
    if az < 0
        az = az + 2*pi; 
    end
end