function InputStruct = interpret_contrast_list_str(InputStruct,modelparam,analysis_model,contrast_list_str)

% split_info.cond


if ~isempty(modelparam)
    
    if(isempty(strfind(modelparam,' ')))
         paramlist = cellstr(modelparam);
    else paramlist = regexp(modelparam,' ','split');
    end
else
    paramlist=[];
end


for ksub = 1:numel(InputStruct)
    for krun = 1:numel(InputStruct(ksub).run)
        clear split_info
        hdr = load_nii_hdr([InputStruct(ksub).run(krun).Output_nifti_file_path '/' InputStruct(ksub).run(krun).Output_nifti_file_prefix '_baseproc.nii']);
        InputStruct(ksub).run(krun).Nt = hdr.dime.dim(5);
        
        filename = InputStruct(ksub).run(krun).split_info_file;
        try
            load(filename);
        catch
            sge_exit(100,sprintf('loading the split file info: %s The file might be corrupted!',filename));
        end
        
        if ~isfield(split_info,'TR_MSEC')
            sge_exit('Please specify TR in the split info file: %s',filename);
        end
        
        
        
        % read command line switches and to the split_info structure
        for paramcount = 1:2:(fix(length(paramlist)/2)*2)
            temp = str2double(paramlist{paramcount+1});
            if isnan(temp)
                split_info.(paramlist{paramcount}) = paramlist{paramcount+1};
            else
                split_info.(paramlist{paramcount}) = temp;
            end
        end
        
        
        DROP_first = InputStruct(ksub).run(krun).DROP_first;
        
        if strcmpi(contrast_list_str,'NONE') % if contrast is not provided use default
            if     strcmpi(split_info.type,'event')
                contrast_list_str = '10';
            elseif strcmpi(split_info.type,'block')
                contrast_list_str = '12';
            elseif strcmpi(split_info.type,'nocontrast')
                contrast_list_str = '00';
            end
        end
        
        % convert contrast string to contrast matrix
        split_info = extract_contrast_list(split_info,contrast_list_str);
        Contrast_List = split_info.Contrast_List;
        
        
        
        if strcmpi(split_info.type,'event') % if contrast is provided but it is not versus baseline for event-related data exit with error
            for i = 1:size(split_info.Contrast_List,1)
                if min(split_info.Contrast_List(i,:))~=0
                    sge_error(100,'Event-related optimization is only allowed for the task vs baseline contrast,e.g. 1-0, 2-0 ');
                end
            end
        end
        
        
        if isfield(split_info,'cond')   % Check whether it is the new version of split info
            
            if all(split_info.cond(1).blklength==0) % check whether it
                split_info.type = 'event';
                for ktemp = 1:length(split_info.cond)
                    split_info.cond(ktemp).blklength = zeros(size(split_info.cond(ktemp).onsetlist));
                end
            else
                split_info.type = 'block';
            end       
            
            % default do not convolve user provided design_mat
            if ~isfield(split_info,'convolve')
                split_info.convolve = 0;
            end
            % check if the custom design matrix is provided
            % if not build design matrix
            if ~isfield(split_info,'design_mat')
                split_info.design_mat = [];
            end
            if split_info.convolve==1
                split_info.design_mat = design_to_hrf( split_info.design_mat, split_info.TR_MSEC/1000, [5.0 15.0] );
            end
            % check whether if in the design matrix drop_first is considered or
            % not, if design_matrix is larger than
            if size(split_info.design_mat,1)==(InputStruct(ksub).run(krun).Nt+DROP_first)
                split_info.design_mat = split_info.design_mat(DROP_first+1:end,:);
            end
            
            
            
            % If name is not provided use numbers for the conditions
            % If baseline condition is not defined those scans which do not
            % belong to a condition will be considered as baseline.
            num_condition = length(split_info.cond);
            split_info.baseline_index = 0;
            for i = 1:num_condition
                if ~isfield(split_info.cond(i),'name') % Name is optional
                    split_info.cond(i).name = num2str(i);
                else
                    split_info.cond(i).name = lower(split_info.cond(i).name);
                end
                if strcmp(split_info.cond(i).name,'baseline')
                    split_info.baseline_index = i;
                end
            end
            
            
            
            % if units is not defined consider split_info.unit = TR
            if ~isfield(split_info,'unit')
                split_info.unit = 'TR';
            end
            % convert sec to msec
            if strcmpi(split_info.unit,'sec')
                for i = 1:length(split_info.cond)
                    split_info.cond(i).onsetlist    = split_info.cond(i).onsetlist*1000;
                    split_info.cond(i).blklength = split_info.cond(i).blklength*1000;
                end
            end
            % convert TR to msec
            if strcmpi(split_info.unit,'TR')
                for i = 1:length(split_info.cond)
                    split_info.cond(i).onsetlist    = split_info.cond(i).onsetlist*split_info.TR_MSEC;
                    split_info.cond(i).blklength = split_info.cond(i).blklength*split_info.TR_MSEC;
                end
            end
            
            
            % compensate the offset in the split_info file
            % DROP_first
            for i = 1:length(split_info.cond)
                split_info.cond(i).onsetlist = split_info.cond(i).onsetlist - DROP_first*split_info.TR_MSEC;
                % remove those onsets that in the first DROP_first scans
                ind_temp = split_info.cond(i).onsetlist<0;
                split_info.cond(i).onsetlist(ind_temp) = [];
                split_info.cond(i).blklength(ind_temp) = [];
            end
            
            % Generate design_mat
            
            Nsubs  = max([1 round(split_info.TR_MSEC/100)]); % (#samples per TR) catch in case TR<100ms
            design = zeros( InputStruct(ksub).run(krun).Nt*Nsubs, length(split_info.cond));     % initialize design matrix
            % index allocates onsets to appropriate design-points
            
            for cond_counter = 1:length(split_info.cond)
                for onset_counter = 1:length(split_info.cond(cond_counter).onsetlist)
                    st = round(split_info.cond(cond_counter).onsetlist(onset_counter)./(split_info.TR_MSEC/Nsubs));
                    ed = st + round(split_info.cond(cond_counter).blklength(onset_counter)./(split_info.TR_MSEC/Nsubs));
                    st = st + 1;
                    ed = ed + 1;
                    if ed>size(design,1)
                        ed = size(design,1);
                    end
                    design(st:ed,cond_counter) = 1;
                end
            end
            % convolve with HRF (must convert into seconds!)
            design_mat = design_to_hrf( design, (split_info.TR_MSEC/Nsubs)/1000, [5.0 15.0] );
            if split_info.baseline_index~=0
                design_mat(:,split_info.baseline_index) = [];
            end
            design_mat = design_mat(1:InputStruct(ksub).run(krun).Nt*Nsubs,:);
            design_mat = bsxfun(@minus,design_mat,mean(design_mat));
            % now, subsample back to get HRF at actual fMRI sampling rate
            split_info.design_mat = [design_mat( round(Nsubs/2): Nsubs : end, : ) split_info.design_mat ];
            
            
            
            
            % change the units to scan for block design data
            % it remains milisecond for event-related
            if strcmpi(split_info.type,'block')
                for i = 1:length(split_info.cond)
                    
                    split_info.cond(i).onsetlist = (split_info.cond(i).onsetlist/split_info.TR_MSEC)+1; %ADD 1 CHANGE ORIGIN of ONSET to 1
                    split_info.cond(i).blklength = (split_info.cond(i).blklength/split_info.TR_MSEC); %
                end
            end
            
            
            
            %             %%%%%%%%%%%%%% FOR TEST REMOVE PLEASE
            %             if strcmpi(split_info.type,'block')
            %                 split_info.design_mat = split_info.design_mat(:,2)-split_info.design_mat(:,1);
            %             end
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            if strcmpi(split_info.type,'block')   % if type=block automatically split the data in half/
                split_info = split_half(split_info,InputStruct(ksub).run(krun).Nt); % create single and group split halves
            end
            
        end
        
        % convert old version to new version (event-related)
        if isfield(split_info,'onsetlist') % This part is for the compatibility for older version version event_related
            split_info.cond(1).onsetlist = split_info.onsetlist;
            split_info.cond(1).blklength = ones(size(split_info.onsetlist));
        end
        
        
        % add signle and group field to be used for signle subject or group
        % analysis, for the group analysis we do not need splitting, so it
        % merges the splitting content.
        
        if isfield(split_info,'idx_cond1_sp1') % This part is for the compatibility for older version version
            
            for i = 1:size(split_info.Contrast_List,1)
                split_info.single.idx_cond(i,1).sp1 = split_info.(sprintf('idx_cond%d_sp1',split_info.Contrast_List(i,1)));
                split_info.single.idx_cond(i,1).sp2 = split_info.(sprintf('idx_cond%d_sp2',split_info.Contrast_List(i,1)));
                split_info.single.idx_cond(i,2).sp1 = split_info.(sprintf('idx_cond%d_sp1',split_info.Contrast_List(i,2)));
                split_info.single.idx_cond(i,2).sp2 = split_info.(sprintf('idx_cond%d_sp2',split_info.Contrast_List(i,2)));
                split_info.group.idx_cond(i,1).sp   = [split_info.single.idx_cond(i,1).sp1 split_info.single.idx_cond(i,1).sp2];
                split_info.group.idx_cond(i,2).sp   = [split_info.single.idx_cond(i,2).sp1 split_info.single.idx_cond(i,2).sp2];
                split_info.group.unbalanced_idx_cond(i,1).sp = split_info.group.idx_cond(i,1).sp;
                split_info.group.unbalanced_idx_cond(i,2).sp = split_info.group.idx_cond(i,2).sp;
            end
        end
        
        
        
        % Now check whether the final split_info is matched with
        % analysis_model
        InputStruct(ksub).run(krun).split_info =  check_split_info(split_info,analysis_model,filename);
        
        
        
        
    end
end







function split_info = extract_contrast_list(split_info,contrast_list_str) % create single and group split halves

   
if strfind(contrast_list_str,'-')

    ind=find(contrast_list_str==',');
    contrast_list_str(ind) = ' ';
    if(isempty(strfind(contrast_list_str,' ')))
         Contrast_temp = cellstr(contrast_list_str);
    else Contrast_temp = regexp(contrast_list_str,' ','split');
    end
    
    num_contrast = length(Contrast_temp);
    for k = 1:num_contrast
        
        if(isempty(strfind(Contrast_temp{k},'-')))
             Contrast_name = cellstr(Contrast_temp{k});
        else Contrast_name = regexp(Contrast_temp{k},'-','split');
        end
        
        for n = 1:2
            Contrast_name{n} = strrep(Contrast_name{n},' ',''); % remove extra space
        end
        
        Contrast_list{k,1} = Contrast_name{1};
        Contrast_list{k,2} = Contrast_name{2};
               
        split_info = interpret_combined_contrast(split_info,Contrast_list{k,1});
        split_info = interpret_combined_contrast(split_info,Contrast_list{k,2});
        
        for i = 1:length(split_info.cond)
            Condition_name{i} = split_info.cond(i).name;
        end      
        
        ind1 = find(ismember(lower(Condition_name),lower(Contrast_list{k,1})));
        ind2 = find(ismember(lower(Condition_name),lower(Contrast_list{k,2})));
        if isempty(ind1)
            ind1 = str2num(Contrast_list{k,1});
        end
        if isempty(ind2)
            ind2 = str2num(Contrast_list{k,2});
        end
        if strcmpi(Contrast_list{k,1},'baseline')
            ind1 = 0;
        end
        if strcmpi(Contrast_list{k,2},'baseline')
            ind2 = 0;
        end
        Contrast(k,:) =[ind1 ind2];
    end
    split_info.Contrast_List = Contrast;
else
    ind=find(contrast_list_str==',');
    contrast_list_str(ind) = ' ';
    
    if(isempty(strfind(contrast_list_str,' ')))
         Contrast_name = cellstr(contrast_list_str);
    else Contrast_name = regexp(contrast_list_str,' ','split');
    end
    
    for k = 1:length(Contrast_name)
        Contrast(k,:) = [str2num(Contrast_name{k}(1)) str2num(Contrast_name{k}(2))];
    end
    split_info.Contrast_List = Contrast;
end

function split_out = interpret_combined_contrast(split_info,Cond_in)
split_out = split_info;
if strfind(Cond_in,'+')
    Cond_temp = Cond_in;
    Cond_temp(Cond_temp=='+') = '%'; % For Octave Compatibility, Octave is sensitive to +
    
    if(isempty(strfind(Cond_temp,'%')))
         Cond_name = cellstr(Cond_temp);
    else Cond_name = regexp(Cond_temp,'%','split');
    end
    
    for i = 1:length(split_info.cond)
        Condition_name{i} = split_info.cond(i).name;
    end
    onsetlist = [];
    blklength = [];
    for i = 1:length(Cond_name)
        ind       = find(ismember(lower(Condition_name),lower(Cond_name{i})));
        if isempty(ind)
            ind = str2num(Cond_name{i});
        end
        onsetlist = [onsetlist;split_info.cond(ind).onsetlist(:)];
        blklength = [blklength;split_info.cond(ind).blklength(:)];
    end
    n = length(split_out.cond);
    split_out.cond(n+1).onsetlist = onsetlist';
    split_out.cond(n+1).blklength = blklength';
    split_out.cond(n+1).name      = lower(Cond_in);
end

function split_info = split_half(split_info,Nt)

Contrast_List = split_info.Contrast_List;

if split_info.baseline_index==0
    all_stim = [];
    for i = 1:length(split_info.cond)
        for k = 1:length(split_info.cond(i).onsetlist)
            all_stim = [all_stim split_info.cond(i).onsetlist(k):split_info.cond(i).onsetlist(k)+split_info.cond(i).blklength(k)-1];
        end
    end
    split_info.baseline = sort(setdiff(1:Nt,all_stim));   % whatever left in time-series is considered as baseline
else
    all_stim = [];
    for i = split_info.baseline_index
        for k = 1:length(split_info.cond(i).onsetlist)
            all_stim = [all_stim split_info.cond(i).onsetlist(k):split_info.cond(i).onsetlist(k)+split_info.cond(i).blklength(k)-1];
        end
    end
    split_info.baseline = all_stim;   % whatever left in time-series is considered as baseline
end

for i = 1:size(Contrast_List,1)
    for condition = 1:2
        split_info.group.idx_cond(i,condition).sp = [];
        if Contrast_List(i,condition)~=0
            for k = 1:length(split_info.cond(Contrast_List(i,condition)).onsetlist)
                split_info.group.idx_cond(i,condition).sp = [split_info.group.idx_cond(i,condition).sp split_info.cond(Contrast_List(i,condition)).onsetlist(k):split_info.cond(Contrast_List(i,condition)).onsetlist(k)+split_info.cond(Contrast_List(i,condition)).blklength(k)-1];
            end
        else
            split_info.group.idx_cond(i,condition).sp = [split_info.group.idx_cond(i,condition).sp split_info.baseline];
        end
        split_info.group.idx_cond(i,condition).sp = sort(split_info.group.idx_cond(i,condition).sp);
    end
end

for i = 1:size(Contrast_List,1)
    
    max_sp = ceil(Nt/2);
    
    % First split
    class1 = split_info.group.idx_cond(i,1).sp(split_info.group.idx_cond(i,1).sp<=max_sp);
    class2 = split_info.group.idx_cond(i,2).sp(split_info.group.idx_cond(i,2).sp<=max_sp);
    [class1,class2] = make_two_class_similar_length(class1,class2);
    split_info.single.idx_cond(i,1).sp1 = class1;
    split_info.single.idx_cond(i,2).sp1 = class2;
    
    class1 = split_info.group.idx_cond(i,1).sp(split_info.group.idx_cond(i,1).sp>max_sp);
    class2 = split_info.group.idx_cond(i,2).sp(split_info.group.idx_cond(i,2).sp>max_sp);
    [class1,class2] = make_two_class_similar_length(class1,class2);
    split_info.single.idx_cond(i,1).sp2 = class1;
    split_info.single.idx_cond(i,2).sp2 = class2;
    
end
for i = 1:size(Contrast_List,1)
    split_info.group.unbalanced_idx_cond(i,1).sp = split_info.group.idx_cond(i,1).sp;
    split_info.group.unbalanced_idx_cond(i,2).sp = split_info.group.idx_cond(i,2).sp;   
    class1 = split_info.group.idx_cond(i,1).sp;
    class2 = split_info.group.idx_cond(i,2).sp;
    [class1,class2] = make_two_class_similar_length(class1,class2);
    split_info.group.idx_cond(i,1).sp = class1;
    split_info.group.idx_cond(i,2).sp = class2;   
end

function [class_out1,class_out2] = make_two_class_similar_length(class1,class2)

class_out1 = class1;
class_out2 = class2;

dist = abs(bsxfun(@minus,class1,class2'));
[r,c] = size(dist);
if c>r
    [tmp,ind_delete] = sort(min(dist),'descend');
    class_out1 = class1(ind_delete(c-r+1:end));
    class_out2 = class2;
end
if r>c
    [tmp,ind_delete] = sort(min(dist'),'descend');
    class_out1 = class1;
    class_out2 = class2(ind_delete(r-c+1:end));
end
class_out1 = sort(class_out1);
class_out2 = sort(class_out2);

    
    
function split_info_out = check_split_info(split_info,analysis_model,filename)

analysis_model = upper(analysis_model);
split_info_out = split_info;


switch analysis_model
    case 'LDA'
        if ~isfield(split_info,'drf')
            error('Please specify drf in the split info file: %s',filename);
        end
        if ~isfield(split_info,'type')
            warning('Please specify type in the split info file: %s, automatically type=block is considered for the LDA model',filename);
            split_info_out.type = 'block';
        end
        
    case 'GNB'
        if ~isfield(split_info,'decision_model')
            error('Please specify decision_model in the split info file: %s',filename);
        end
        if ~isfield(split_info,'type')
            warning('Please specify type in the split info file: %s, automatically type=block is considered for the GNB model',filename);
            split_info_out.type = 'block';
        end
    case 'GLM'
        if ~isfield(split_info,'design_mat')
            error('Please specify design_mat in the split info file: %s',filename);
        end
        if ~isfield(split_info,'convolve')
            error('Please specify convolve in the split info file: %s',filename);
        end
    case 'ERCVA'
        if ~isfield(split_info,'Nsplit')
            error('Please specify Nsplit in the split info file: %s',filename);
        end
        if ~isfield(split_info,'WIND')
            error('Please specify WIND in the split info file: %s',filename);
        end
        if ~isfield(split_info,'drf')
            error('Please specify drf in the split info file: %s',filename);
        end
        if ~isfield(split_info,'subspace')
            error('Please specify subspace in the split info file: %s',filename);
        end
        if ~isfield(split_info,'type')
            warning('Please specify type in the split info file: %s, automatically type=event is considered for the erCVA model',filename);
            split_info_out.type = 'event';
        end
    case 'ERGNB'
        if ~isfield(split_info,'Nsplit')
            error('Please specify Nsplit in the split info file: %s',filename);
        end
        if ~isfield(split_info,'WIND')
            error('Please specify WIND in the split info file: %s',filename);
        end
        if ~isfield(split_info,'type')
            warning('Please specify type in the split info file: %s, automatically type=event is considered for the erGNB model',filename);
            split_info_out.type = 'event';
        end
    case 'ERGLM'
        if ~isfield(split_info,'type')
            warning('Please specify type in the split info file: %s, automatically type=event is considered for the erGLM model',filename);
            split_info_out.type = 'event';
        end
    case 'SCONN'
        if ~isfield(split_info,'seed_name')
            error('Please specify seed_name in the split info file: %s',filename);
        end
        if ~isfield(split_info,'spm')
            error('Please specify spm in the split info file: %s',filename);
        end
end




