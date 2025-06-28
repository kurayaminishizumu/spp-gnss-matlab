function [obs_data, approx_pos_xyz, time_system] = read_rinex_obs_gps(filename)
obs_data = {};
approx_pos_xyz = [];
time_system = 'GPS'; 

fid = fopen(filename, 'r');
if fid == -1
    error('无法打开观测文件: %s', filename);
end

current_epoch_data = [];
header_done = false;
obs_types_gps = {}; 

line_buffer_obs_types = {}; 

while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line) || isempty(strtrim(line))
        if ~isempty(line_buffer_obs_types) && ~header_done
            process_obs_types_buffer_internal();
        end
        if ~ischar(line) 
            break; 
        end
        continue; 
    end

    label = '';
    if length(line) >= 61 
        label = strtrim(line(61:end));
    end

    if ~header_done
        if contains(label, 'APPROX POSITION XYZ')
            process_obs_types_buffer_internal();
            approx_pos_xyz = sscanf(line(1:min(60, length(line))), '%f');
            if length(approx_pos_xyz) ~= 3
                warning('APPROX POSITION XYZ 解析不完整或格式不正确，请检查文件或代码。');
                approx_pos_xyz = []; 
            end
        elseif contains(label, 'SYS / # / OBS TYPES')
            current_sys_char = line(1);
            if ~isempty(line_buffer_obs_types)
                first_buffered_sys_char = line_buffer_obs_types{1}(1);
                if current_sys_char ~= first_buffered_sys_char || ... 
                   (length(line_buffer_obs_types{end}) >= 61 && ... 
                    ~strcmp(strtrim(line_buffer_obs_types{end}(61:end)), 'SYS / # / OBS TYPES'))
                    
                    process_obs_types_buffer_internal(); 
                end
            end
            line_buffer_obs_types{end+1} = line; 
        elseif contains(label, 'END OF HEADER')
            process_obs_types_buffer_internal(); 
            header_done = true;
            if isempty(obs_types_gps)
                warning('RINEX头中未找到GPS的 SYS / # / OBS TYPES。');
            end
            line_buffer_obs_types = {}; 
            continue; 
        else
            process_obs_types_buffer_internal();
        end
    else 
        if line(1) == '>' 
            if ~isempty(current_epoch_data) && isfield(current_epoch_data, 'sv') && ~isempty(fieldnames(current_epoch_data.sv))
                obs_data{end+1} = current_epoch_data; 
            end
            current_epoch_data = struct();
            parts = strsplit(strtrim(line(2:29))); 
            year = str2double(parts{1});
            month = str2double(parts{2});
            day = str2double(parts{3});
            hour = str2double(parts{4});
            minute = str2double(parts{5});
            second = str2double(parts{6});

            dt_obj = datetime(year, month, day, hour, minute, second);
            [gps_week, ~, gps_sow] = date_to_gps_week_sow(dt_obj);
            current_epoch_data.time = gps_sow;
            current_epoch_data.week = gps_week;
            current_epoch_data.sv = struct();
        elseif line(1) >= 'A' && line(1) <= 'Z'
            if isempty(current_epoch_data) || ~isfield(current_epoch_data, 'time')
                continue;
            end
            sys_char = line(1);
            if sys_char ~= 'G' 
                continue;
            end
            prn_str = strtrim(line(1:3)); 

            if isempty(obs_types_gps)
                if sys_char == 'G'
                    warning('在处理GPS数据行时，obs_types_gps 为空，无法解析观测值。PRN: %s', prn_str);
                end
                continue; 
            end
            
            sv_obs = struct();
            current_char_pos = 4; 
            for k_obs_type = 1:length(obs_types_gps)
               
                if current_char_pos + 13 <= length(line)
                    val_str = strtrim(line(current_char_pos : current_char_pos+13));
                    if ~isempty(val_str)
                        val = str2double(val_str);
                        obs_type_name = obs_types_gps{k_obs_type};
                        if strcmp(obs_type_name, 'C1C') || strcmp(obs_type_name, 'P1') || strcmp(obs_type_name, 'C1P') || strcmp(obs_type_name, 'C2P') || strcmp(obs_type_name, 'C2W') || strcmp(obs_type_name, 'C5Q')
                           sv_obs.(obs_type_name) = val;
                        end
                    end
                else
                    break; 
                end
                current_char_pos = current_char_pos + 16;
            end
            if ~isempty(fieldnames(sv_obs))
                 current_epoch_data.sv.(prn_str) = sv_obs;
            end
        end
    end
end
if ~isempty(line_buffer_obs_types)
    process_obs_types_buffer_internal();
end

if ~isempty(current_epoch_data) && isfield(current_epoch_data, 'sv') && ~isempty(fieldnames(current_epoch_data.sv))
    obs_data{end+1} = current_epoch_data; 
end
fclose(fid);

% --- 嵌套函数，用于处理观测类型缓冲 ---
    function process_obs_types_buffer_internal()
        if isempty(line_buffer_obs_types), return; end

        first_line_in_buffer = line_buffer_obs_types{1};
        sys_char_of_buffer = first_line_in_buffer(1);

        if sys_char_of_buffer == 'G' 
            num_obs_types_declared_for_G = 0;
            if length(first_line_in_buffer) >=6 
                num_obs_types_declared_for_G = str2double(strtrim(first_line_in_buffer(4:6)));
            else
                warning('SYS / # / OBS TYPES 行过短，无法读取观测类型数量: %s', first_line_in_buffer);
                line_buffer_obs_types = {}; 
                return;
            end
            
            full_obs_types_string_for_G = '';
            for i_line_buf = 1:length(line_buffer_obs_types)
                current_buffered_line = line_buffer_obs_types{i_line_buf};
                if current_buffered_line(1) ~= sys_char_of_buffer || ...
                   (length(current_buffered_line) >= 61 && ~strcmp(strtrim(current_buffered_line(61:end)), 'SYS / # / OBS TYPES'))
                    break;
                end

                start_col_types = 8;
                end_col_types = min(60, length(current_buffered_line)); 
                
                if end_col_types >= start_col_types
                    types_on_this_line = strtrim(current_buffered_line(start_col_types:end_col_types));
                    if ~isempty(types_on_this_line)
                        if isempty(full_obs_types_string_for_G)
                            full_obs_types_string_for_G = types_on_this_line;
                        else
                            full_obs_types_string_for_G = [full_obs_types_string_for_G, ' ', types_on_this_line];
                        end
                    end
                end
            end
            parsed_types_list_for_G_temp = strsplit(full_obs_types_string_for_G, ' ');
            parsed_types_list_for_G = {};
            for pt_idx = 1:length(parsed_types_list_for_G_temp)
                cleaned_type = strtrim(parsed_types_list_for_G_temp{pt_idx});
                if ~isempty(cleaned_type)
                    parsed_types_list_for_G{end+1} = cleaned_type; 
                end
            end
            
            if numel(parsed_types_list_for_G) >= num_obs_types_declared_for_G
                obs_types_gps = parsed_types_list_for_G(1:num_obs_types_declared_for_G);
            else
                warning('RINEX Header: GPS SYS/#/OBS TYPES - Declared %d types, but only parsed %d from all lines. Using parsed count.', ...
                        num_obs_types_declared_for_G, numel(parsed_types_list_for_G));
                obs_types_gps = parsed_types_list_for_G; 
            end
        end
        line_buffer_obs_types = {};
    end
% --- 结束嵌套函数 ---
end

% 辅助函数: datetime to GPS week and SOW
function [gps_week, day_of_week, seconds_of_week] = date_to_gps_week_sow(dt_obj)
    gps_epoch_dt = datetime(1980, 1, 6, 0, 0, 0, 'TimeZone', 'UTC');
    if isempty(dt_obj.TimeZone) || ~strcmp(dt_obj.TimeZone, 'UTC')
        dt_obj.TimeZone = 'UTC'; 
    end
    
    time_diff_seconds = posixtime(dt_obj) - posixtime(gps_epoch_dt);
    
    gps_week = floor(time_diff_seconds / (7 * 24 * 3600));
    seconds_into_current_week = mod(time_diff_seconds, (7 * 24 * 3600));
    
    day_of_week = floor(seconds_into_current_week / (24 * 3600)); 
    seconds_of_week = seconds_into_current_week;
end