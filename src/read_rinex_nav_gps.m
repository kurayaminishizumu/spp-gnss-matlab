function [ephemeris_data, iono_params] = read_rinex_nav_gps(filename)

ephemeris_data = struct();
iono_params = struct('alpha', zeros(1,4), 'beta', zeros(1,4));

fid = fopen(filename, 'r');
if fid == -1
    error('无法打开导航文件: %s', filename);
end

header_done = false;
current_eph = [];
line_counter = 0; 

while ~feof(fid)
    line = fgetl(fid);
    
    if isempty(strtrim(line)), continue; end 
    
    label = '';
    if length(line) >= 61 
        label = strtrim(line(61:end));
    end

    if ~header_done
        % --- 处理电离层校正参数 ---
        if contains(label, 'IONOSPHERIC CORR') 
             sys_char_iono = line(1); 
             parsed_vals = zeros(1,4);
             try
                 parsed_vals(1) = str2double_sci_robust(line, 6, 22);
                 parsed_vals(2) = str2double_sci_robust(line, 23, 39);
                 parsed_vals(3) = str2double_sci_robust(line, 40, 57);
                 parsed_vals(4) = str2double_sci_robust(line, 58, 75);
             catch ME_iono_parse
                 warning('Error parsing IONOSPHERIC CORR line: %s. Detail: %s. Line: %s', ME_iono_parse.message, ME_iono_parse.identifier, line);
                 parsed_vals = NaN(1,4);
             end
            
             if sys_char_iono == 'G' 
                 if contains(line, 'GPSA') 
                     iono_params.alpha = parsed_vals;
                 elseif contains(line, 'GPSB') 
                     iono_params.beta = parsed_vals;
                 end
             end
             
        % --- 处理头文件结束标志 ---
        elseif contains(label, 'END OF HEADER')
            header_done = true;
            continue; 
        end
    else 
        if startsWith(line, 'G') 
            prn_str = line(2:3);
            prn_val = str2double(prn_str);
            if isnan(prn_val) 
                continue;
            end
            if length(line) < 23 
                continue;
            end
            
            try
                year = str2double(line(5:8));
                month = str2double(line(10:11));
                day = str2double(line(13:14));
                hour = str2double(line(16:17));
                minute = str2double(line(19:20));
                second = str2double(line(22:23)); 
                
                if isnan(year) || isnan(month) || isnan(day) || isnan(hour) || isnan(minute) || isnan(second) || ...
                   month < 1 || month > 12 || day < 1 || day > 31 || hour < 0 || hour > 23 || minute < 0 || minute > 59 || second < 0 || second >= 61 % Seconds can be 60 for leap second
                    continue;
                end
            catch ME_toc_parse
                continue;
            end
            if length(line) >= 80
                line_counter = 1; 
                current_eph = struct();
                current_eph.PRN = prn_val;
                
                current_eph.toc_year = year;
                current_eph.toc_month = month;
                current_eph.toc_day = day;
                current_eph.toc_hour = hour;
                current_eph.toc_minute = minute;
                current_eph.toc_second = second;
                
                current_eph.af0 = str2double_sci_robust(line, 24, 42);
                current_eph.af1 = str2double_sci_robust(line, 43, 61);
                current_eph.af2 = str2double_sci_robust(line, 62, 80);
                
            else
                continue; 
            end
        elseif line_counter > 0 && line_counter < 8 && isstruct(current_eph) 
            line_counter = line_counter + 1;
            
            try
                vals(1) = str2double_sci_robust(line, 5, 23);  
                vals(2) = str2double_sci_robust(line, 24, 42);
                vals(3) = str2double_sci_robust(line, 43, 61); 
                vals(4) = str2double_sci_robust(line, 62, 80); 
            catch ME_eph_line_parse
                warning('Error parsing ephemeris data line %d for PRN G%02d: %s. Detail: %s. Line: %s', ...
                        line_counter, current_eph.PRN, ME_eph_line_parse.message, ME_eph_line_parse.identifier, line);
                vals = NaN(1,4); 
            end

            switch line_counter
                case 2 
                    current_eph.IODE = vals(1); current_eph.Crs = vals(2);
                    current_eph.Delta_n = vals(3); current_eph.M0 = vals(4);
                case 3 
                    current_eph.Cuc = vals(1); current_eph.e = vals(2);
                    current_eph.Cus = vals(3); current_eph.sqrt_A = vals(4);
                case 4 
                    current_eph.Toe = vals(1); current_eph.Cic = vals(2);
                    current_eph.OMEGA0 = vals(3); current_eph.Cis = vals(4);
                case 5 
                    current_eph.i0 = vals(1); current_eph.Crc = vals(2);
                    current_eph.omega = vals(3); current_eph.OMEGA_DOT = vals(4);
                case 6
                    current_eph.IDOT = vals(1);
                    current_eph.CodesL2 = vals(2); 
                    current_eph.GPSWeek = vals(3);
                    current_eph.L2Pflag = vals(4); 
                case 7 
                    current_eph.SV_accuracy = vals(1); 
                    current_eph.SV_health = vals(2);
                    current_eph.TGD = vals(3); 
                    current_eph.IODC = vals(4); 
                    
                    eph_field = sprintf('G%02d', current_eph.PRN);
                    ephemeris_data.(eph_field) = current_eph; 
                    line_counter = 0; 
                    current_eph = []; 
            end
        end
    end
end
fclose(fid);
end

% --- 辅助函数 ---
function val = str2double_sci_robust(line_str, start_col, end_col)
    val = NaN; 
    if end_col > length(line_str)
        return; 
    end
    sub_str = strtrim(line_str(start_col:end_col));

    if isempty(sub_str)
        return; 
    end

    sub_str = strrep(sub_str, 'D', 'e'); 
    
    parsed_val = str2double(sub_str);
    if ~isnan(parsed_val)
        val = parsed_val;
    end
end