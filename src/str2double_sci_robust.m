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