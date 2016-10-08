function bids_parsejobs( bids_dir, output_dir, level, participant, run_name, contrast, task_design, analysis_model, ndrop, atlasfile )
%
% BIDS_PARSEJOBS: script takes in bids input arguments, then constructs a
% list of files, in preparation of setting up oppni analyses
%
%  Syntax:
%
%     bids_parsejobs( bids_dir, output_dir, level, participant, run_name, contrast, task_design, analysis_model, ndrop, atlasfile )
%

% start out empty
tsvlev=0;

%get the list of all subjects in input file -- allows for missing files
kq=0; sublist=[]; e=dir(bids_dir);
for(i=1:length(e))
    if( ~isempty(strfind(e(i).name,'sub-')) )
        kq=kq+1;
        sublist{kq}=e(i).name;
    end
    % check existence of .tsvs at top level --> should be here!
    if( ~isempty(strfind(e(i).name,run_name)) && ~isempty(strfind(e(i).name,'_bold.json'))  )
        jsonfile = [bids_dir,'/',e(i).name];
    end  
    % check existence of .tsvs at top level --> might not be here
    if( ~isempty(strfind(e(i).name,run_name)) && ~isempty(strfind(e(i).name,'_events.tsv'))  )
        tsvlev  = 1;
    end    
end
% check for sessions -- uses 1st subject as reference
kq=0; seslist=[]; e=dir([bids_dir,'/',sublist{1}]);
for(i=1:length(e))
    if( ~isempty(strfind(e(i).name,'ses-')) )
        kq=kq+1;
        seslist{kq}=e(i).name;
    end  
end
% now get list of files matching run name
if(isempty(seslist))
    kq=0; runlist=[]; 
    e=dir([bids_dir,'/',sublist{1},'/func']);
    for(i=1:length(e))
        if( ~isempty(strfind(e(i).name,run_name)) && ~isempty(strfind(e(i).name,'_bold')) && ~isempty(strfind(e(i).name,'run-')) )
            kq=kq+1;
            runlist{kq}=[e(i).name( strfind(e(i).name,'run-'): strfind(e(i).name,'_bold')-1)];
        end
        % check existence of .tsvs at bottom level --> might not be here
        if( tsvlev ==0 && ~isempty(strfind(e(i).name,run_name)) && ~isempty(strfind(e(i).name,'_events.tsv'))  )
            tsvlev  = 2;
        end            
    end
else
    kq=0; runlist=[];
    for(l=1:length(seslist))
        e=dir([bids_dir,'/',sublist{1},'/',seslist{l},'/func']);
        for(i=1:length(e))
            if( ~isempty(strfind(e(i).name, run_name)) && ~isempty(strfind(e(i).name, '_bold'))  && ~isempty(strfind(e(i).name,'run-'))  )
                kq=kq+1;
                runlist{kq}=[e(i).name( strfind(e(i).name,'run-'): strfind(e(i).name,'_bold')-1)];
            end
            % check existence of .tsvs at bottom level --> might not be here
            if( tsvlev ==0 && ~isempty(strfind(e(i).name,run_name)) && ~isempty(strfind(e(i).name,'_events.tsv'))  )
                tsvlev  = 2;
            end              
        end    
    end
end

if(isempty(runlist)) runlist{1}=''; end

% populate list of fmri IN,OUT, anatomical STRUCT, task files
for(i=1:length(sublist))
    kq=0;
    if(isempty(seslist))
        for(j=1:length(runlist))
            kq=kq+1;
            fmri_in_list{i}{kq}  = [bids_dir,  '/',sublist{i},'/func/', sublist{i},'_task-',run_name,runlist{j},'_bold.nii.gz']; %IN
            fmri_out_list{i}{kq} = [output_dir,'/', sublist{i},'_task-',run_name,runlist{j} ]; %OUT
            e = dir([bids_dir,  '/',sublist{i},'/anat/',sublist{i},'*']); %STRUCT
            struct_list{i}{kq}   = [bids_dir,  '/',sublist{i},'/anat/',e(1).name];
            if(size(struct_list{i}{kq},1)>1) error(['multiple anat files found for ',sublist{i},'.']); end
            % EVENT TSV
            if    (tsvlev==1) tsvfile_list{i}{kq} = [bids_dir, '/task-',run_name,'_events.tsv'];
            elseif(tsvlev==2) tsvfile_list{i}{kq} = [bids_dir,  '/',sublist{i},'/func/', sublist{i},'_task-',run_name,runlist{j},'_events.tsv'];                
            end
        end
    else
        for(j=1:length(seslist))
        for(k=1:length(runlist))
            kq=kq+1;            
            fmri_in_list{i}{kq}  = [bids_dir,  '/',sublist{i},'/',seslist{j},'/func/', sublist{i},'_', seslist{j},'_task-',run_name,runlist{k},'_bold.nii.gz']; %IN
            fmri_out_list{i}{kq} = [output_dir,'/', sublist{i},'_', seslist{j},'_task-',run_name,runlist{k} ]; %OUT
            e = dir([bids_dir,  '/',sublist{i},'/',seslist{j},'/anat/',sublist{i},'_',seslist{j},'*']); %STRUCT
            struct_list{i}{kq}   = [bids_dir,  '/',sublist{i},'/',seslist{j},'/anat/',e(1).name];
            if(size(struct_list{i}{kq},1)>1) error(['multiple anat files found for ',sublist{i},'.']); end
            % EVENT TSV
            if    (tsvlev==1)  tsvfile_list{i}{kq} = [bids_dir, '/task-',run_name,'_events.tsv'];
            elseif(tsvlev==2)  tsvfile_list{i}{kq} = [bids_dir,  '/',sublist{i},'/',seslist{j},'/func/', sublist{i},'_', seslist{j},'_task-',run_name,runlist{k},'_events.tsv'];                
            end
        end
        end
    end
end


if(strcmpi(level,'participant'))
%if individual:
    if(~isempty(participant))
        if(ischar(participant)) participant=str2num(participant); end
        %if ID# provided, keep only this subject [may include multi-runs / sessions]
        fmri_in_list  =  fmri_in_list( participant );
        fmri_out_list = fmri_out_list( participant );
        tsvfile_list  =  tsvfile_list( participant );
        struct_list   =   struct_list( participant );
    end
    %call bids_setupjobs to generate processed data (PART1)
    bids_setupjobs( 'PART1', output_dir, fmri_in_list, fmri_out_list, struct_list, [], ndrop,0, jsonfile, tsvfile_list, [atlasfile], contrast, task_design, analysis_model );
    
elseif(strcmpi(level,'group'))
%if group:
    %call bids_setupjobs (PART2,spnorm,groupmask)
    bids_setupjobs( 'PART2', output_dir, fmri_in_list, fmri_out_list, struct_list, [], ndrop,0, jsonfile, tsvfile_list, [atlasfile], contrast, task_design, analysis_model );
    if( ~isempty(atlasfile) )
        bids_setupjobs( 'SPNORM',output_dir, fmri_in_list, fmri_out_list, struct_list, [], ndrop,0, jsonfile, tsvfile_list, [atlasfile], contrast, task_design, analysis_model );
        bids_setupjobs( 'GMASK', output_dir, fmri_in_list, fmri_out_list, struct_list, [], ndrop,0, jsonfile, tsvfile_list, [atlasfile], contrast, task_design, analysis_model );
    end
end
