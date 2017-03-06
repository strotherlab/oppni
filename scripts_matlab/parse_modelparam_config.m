function struc_base = parse_modelparam_config( in, struc_base )

if( ischar(in) && ~isempty(strfind(in,':')) && ~isempty(strfind(in,'"')) && isempty(strfind(in,'.json')) )
    disp('modelparam seems to be already in string format.');
    C_str = in;
    
elseif( ~isempty(strfind(in,'.json')) )
    
    disp('modelparam is in a json file; loading.');
    fid = fopen(in);
    C_text=textscan(fid,'%s','delimiter','\n');
    fclose(fid);

    C_str=[];
    for(i=1:length(C_text{1})) C_str=[C_str C_text{1}{i}]; end
    
elseif( ~isempty(in) )
    error('invalid modelparam format');
end

if( isempty(in) )
    
    disp('modelparam is empty. nothing added');
    
elseif( ~isempty(strfind(C_str,'{')) && ~isempty(strfind(C_str,'}')) && ~isempty(strfind(C_str,'"')) && ~isempty(strfind(C_str,':')) )
    
    strucj = p_json(C_str);    
    nm = fields(strucj);
    for(i=1:numel(nm))
        if(isfield(struc_base,nm{i})) 
            error('modelparam trying to overwrite pre-existing field'); 
        end
        struc_base.(nm{i}) = strucj.(nm{i});
    end
else
    error('something wrong with modelparam -- does not fit json format!');
end

