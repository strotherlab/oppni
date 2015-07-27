function [split_info_struct] = Parse_Split_Info(split_info_name)
%
% script to read "Split-Info" file, in either .mat or .txt format
% interprets and reformats for compatibility with PRONTO code
%

[path,prefix,ext] = fileparts(split_info_name);

if(strcmp(ext,'.mat')) %% if in original .mat file format...
    %
    try
        load(split_info_name);
    catch
        warning('loading the split file info: %s The file might be corrupted!',split_info_name);
         split_info_struct = [];
         return;
    end    
    % copy over the output structure
    split_info_struct = split_info;
    
else %% if in textfile format...
    
    split_info_struct=[];
    
    fid = fopen(split_info_name);
    if fid==-1
        warning('loading the split file info: %s The file might be corrupted!',split_info_name);
         return;
        split_info_struct = [];
    end
    % read in first line
    tline = fgetl(fid);
    if ~ischar(tline)
        warning('loading the split file info: %s file appears to be empty.',split_info_name);
        split_info_struct = [];
        return;
    end

    %
    n_name     = 0;
    n_onsets   = 0;
    n_duration = 0;

    %% read input file
    while ischar(tline) 

        % [] bracketing fields
        istart = strfind(tline,'[')+1;
        iend   = strfind(tline,']')-1;   
        
        if( ~isempty(istart) && ~isempty(iend) )
               
            % check for general options
            if(~isempty(strfind(upper(tline),'UNIT')))    split_info_struct.unit = tline(istart:iend);
            end
            if(~isempty(strfind(upper(tline),'TR_MSEC'))) split_info_struct.TR_MSEC = str2num( tline(istart:iend) );            
            end        
            if(~isempty(strfind(upper(tline),'TYPE')))    split_info_struct.type = lower(tline(istart:iend));            
            end     
            % seed for functional connectivity
            if(~isempty(strfind(upper(tline),'SEED')))    split_info_struct.seed_name = tline(istart:iend);
            end
        
            % check for condition-specific options
            if(~isempty(strfind(upper(tline),'NAME')))
                n_name = n_name+1;
                if( ~isempty(strfind(tline(istart:iend),'-')) || ~isempty(strfind(tline(istart:iend),'+')) )
                    sge_exit(100,sprintf('split file info: %s condition names cannot have + or - symbols, as these define contrasts.',split_info_name));
                end
                split_info_struct.cond(n_name).name = tline(istart:iend);
            end
            % 
            if(~isempty(strfind(upper(tline),'ONSETS')))
                n_onsets = n_onsets+1;
                e = regexp( tline(istart:iend), ',','split' );
                onsetsArray{n_onsets} = [];
                for e_counter_temp = 1:length(e)
                    onsetsArray{n_onsets} = [onsetsArray{n_onsets} str2num(e{e_counter_temp})];
                end
            end      
            % 
            if(~isempty(strfind(upper(tline),'DURATION')))
                n_duration = n_duration+1;
                e = regexp( tline(istart:iend), ',','split' );                
                durationArray{n_duration} = [];
                for e_counter_temp = 1:length(e)
                     durationArray{n_duration}= [ durationArray{n_duration} str2num(e{e_counter_temp})];
                end
            end                       
        end

        tline = fgetl(fid);
        if isempty(tline)
            tline = fgetl(fid);
        end
    end
    fclose(fid); 
    
    if( (n_name == n_onsets) && (n_name == n_duration) && (n_duration == n_onsets) )
        
        for(n=1:n_name)
            
            if( length(onsetsArray{n}) == length(durationArray{n}) )
                %
                split_info_struct.cond(n).onsetlist = onsetsArray{n};                
                split_info_struct.cond(n).blklength = durationArray{n};
            else
                error('number of ONSETS and DURATION values must be the same for every condition');
            end
        end
    else
        error('every task condition in split-info. requires a NAME, list of ONSETS and DURATIONS');
    end
end


    
