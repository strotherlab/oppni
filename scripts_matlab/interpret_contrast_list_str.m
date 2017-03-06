function InputStruct = interpret_contrast_list_str(InputStruct,modelparam,analysis_model,contrast_list_str)

% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%

for ksub = 1:numel(InputStruct)
    for krun = 1:numel(InputStruct(ksub).run)
        if ~exist(InputStruct(ksub).run(krun).split_info_file,'file')
            InputStruct(ksub).run(krun).split_info = [];
            continue;
        end
        
        clear split_info
        [p,f,e] = fileparts([InputStruct(ksub).run(krun).Input_nifti_filename]);
        if(isempty(strfind(e,'.gz'))) %if not a zip file, read the header direct
            hdr = load_nii_hdr([InputStruct(ksub).run(krun).Input_nifti_filename]);
        else %otherwise need to inflate and load .nii
            v = load_untouch_nii([InputStruct(ksub).run(krun).Input_nifti_filename]);
            hdr=v.hdr; clear v;
        end

        InputStruct(ksub).run(krun).Nt = hdr.dime.dim(5) - InputStruct(ksub).run(krun).DROP_first - InputStruct(ksub).run(krun).DROP_last;
        
        % loading in the split-info filename
        filename = InputStruct(ksub).run(krun).split_info_file;
        [split_info] = Parse_Split_Info(  InputStruct(ksub).run(krun).split_info_file  );
       
        %------------ modification for consistency of condition lists -------------%
        if( ~isempty(contrast_list_str) && ~strcmpi(contrast_list_str,'NONE') )
        
            % conditions being analyzed in pipelines
            conditions_of_interest = contrast_list_str;
            conditions_of_interest = strrep(conditions_of_interest,'-','+');
            conditions_of_interest = strrep(conditions_of_interest,',','+');        
            conditions_of_interest = unique(regexp(conditions_of_interest,'+','split'));

            clear no_interest;
            % (1) drop conditions not present in contrast(s) of interest
            keepfields={};
            for(i=1:length(split_info.cond))
                if( sum(strcmpi(split_info.cond(i).name,conditions_of_interest))==0 )
                    no_interest(i) = true; 
                else
                    no_interest(i) = false;
                    keepfields = [keepfields, {split_info.cond(i).name}];
                end
            end        
            split_info.cond( no_interest ) = []; %% delete fields of no-interest

            % (2) check if all required conditions are present
            for(i=1:length(conditions_of_interest))
                if( sum(strcmpi(conditions_of_interest{i},keepfields))==0 && ~strcmpi(conditions_of_interest{i},'baseline'))
                    sge_exit(100,sprintf('Condition %s in contrast not defined for file: %s',conditions_of_interest{i},InputStruct(ksub).run(krun).split_info_file));
                end
            end
            % and if they all have onsets / etc.
            for(i=1:length(split_info.cond))
                if( isempty( split_info.cond(i).onsetlist ) || isempty( split_info.cond(i).blklength ) )
                    sge_exit(100,sprintf('Condition %s in contrast has no onsets and/or durations in file: %s',split_info.cond(i).name,InputStruct(ksub).run(krun).split_info_file));                
                end
            end
            % (3) reorder all conditions alphabetaically
            [vx ix] = sort( keepfields ); %% reorder alphabetically (since they are unique, ordering is consistent!)
            split_info_temp = split_info;
            for(i=1:length(split_info.cond))
                split_info.cond(i) = split_info_temp.cond(ix(i));
            end
        end
        
        %% READ IN PARAMETER LIST
        if(~isempty(modelparam))
            % appending all parameters
            split_info = parse_modelparam_config( modelparam, split_info );
        end
        
        % TR units must be specified
        if ~isfield(split_info,'TR_MSEC')
            sge_exit(100,sprintf('Please specify split_info.TR_MSEC in the split info file: %s',filename));
        end

        if( ~strcmpi(analysis_model.design_type,'nocontrast') && ~isfield(split_info,'cond') )
            % case with task contrast, no conditions
            sge_exit(100,'you specified contrast-based analysis without any conditions!');

        elseif isfield(split_info,'cond') 
            % case with conditions specified --> currently, we admit for nocontrast/no analysis model
            
            % number of "real" conditions; need to keep track because combined
            % conditions may be specified below
            num_real_condition = length(split_info.cond);

            % convert contrast string to contrast matrix --> numeric indices + 
            % merged conditions for doing combined contrasts
            split_info = extract_contrast_list(split_info,contrast_list_str);

            % Check for event-related: contrast is provided but not versus baseline --> exit with error
            if strcmpi(analysis_model.design_type,'event')
                for i = 1:size(split_info.Contrast_List,1)
                    if min(split_info.Contrast_List(i,:))~=0
                        sge_exit(100,'Event-related optimization only allowed for task v baseline contrast, e.g. A-baseline, A+B-baseline etc. ');
                    end
                end
            end
            
            % If baseline condition is not defined, those scans which do not
            % belong to a condition will be considered as baseline.
            split_info.baseline_index = 0;
            for i = 1:length(split_info.cond)
                split_info.cond(i).name = lower(split_info.cond(i).name);
                if strcmp(split_info.cond(i).name,'baseline')
                    split_info.baseline_index = i;
                end
            end
            
            % checking if onset units specified
            if ~isfield(split_info,'unit')
                error('need to specify units for task onsets!');
            elseif( ~strcmpi(split_info.unit,'sec') && ~strcmpi(split_info.unit,'msec') && ~strcmpi(split_info.unit,'TR') )
                error('split_info.unit must be either "sec" (seconds) "msec" (milliseconds) or "TR" (repetition time)');
            end
            % convert sec to msec
            if strcmpi(split_info.unit,'sec')
                for i = 1:length(split_info.cond)
                    split_info.cond(i).onsetlist = split_info.cond(i).onsetlist*1000;
                    split_info.cond(i).blklength = split_info.cond(i).blklength*1000;
                end
            end
            % convert TR to msec
            if strcmpi(split_info.unit,'TR')
                for i = 1:length(split_info.cond)
                    split_info.cond(i).onsetlist = split_info.cond(i).onsetlist*split_info.TR_MSEC;
                    split_info.cond(i).blklength = split_info.cond(i).blklength*split_info.TR_MSEC;
                end
            end
            
            % compensate the offset in the split_info file, using DROP_first
            for i = 1:length(split_info.cond)
                split_info.cond(i).onsetlist = split_info.cond(i).onsetlist - InputStruct(ksub).run(krun).DROP_first*split_info.TR_MSEC;
                % remove those onsets that occur in the first DROP_first scans
                ind_temp = split_info.cond(i).onsetlist<0;
                split_info.cond(i).onsetlist(ind_temp) = [];
                split_info.cond(i).blklength(ind_temp) = [];
                % check if any blocks fully overrun the fmri run (+4s)
                onsmax=max(split_info.cond(i).onsetlist) + 4000;
                if(onsmax > (InputStruct(ksub).run(krun).Nt*split_info.TR_MSEC)) 
                    error(['condition "',split_info.cond(i).name,'" has blocks that start near (or after) end of fMRI run! check your onsets.']); 
                end
            end
            
            %--------------- start creating design mat ----------------%
            
            % .we create design matrix, with all ORIGINAL conditions included in 
            % .split_info, EXCLUDING combined contrasts. start off with sub-sampled
            % .matrix and interpolate to correct TR intervals
            % .NB: if a "baseline" is specified, this is excluded to ensure non-singular design
            
            Nsubs  = max([1 round(split_info.TR_MSEC/100)]); % subsampled design matrix (#samples per TR) catch in case TR<100ms
            design = zeros( InputStruct(ksub).run(krun).Nt*Nsubs, num_real_condition);     % initialize design matrix
            
            for cond_counter = 1:num_real_condition %% for each real (non-combined) condition
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
            % convolve with HRF - design_to_hrf function requires that we convert TR to sec.
            design_mat = design_to_hrf( design, (split_info.TR_MSEC/Nsubs)/1000, [5.0 15.0] );
            % ensures "baseline" condition is excluded from full regression model
            if split_info.baseline_index~=0
                design_mat(:,split_info.baseline_index) = [];
            end
            design_mat = design_mat(1:InputStruct(ksub).run(krun).Nt*Nsubs,:);
            design_mat = bsxfun(@minus,design_mat,mean(design_mat));
            % now, subsample back to get HRF at "real" fMRI sampling rate
            split_info.design_mat = [design_mat( round(Nsubs/2): Nsubs : end, : )];
            split_info.design_mat=licols(split_info.design_mat);
            
            %--------------- done creating design mat ----------------%
            
            %% for block design, reformat onsets into TRs
            %
            if strcmpi(analysis_model.design_type,'block')
                % change the units to scan for block design data
                % *it remains milisecond for event-related
                for i = 1:length(split_info.cond)
                    % onset is allocated to TR x, if:   onset>(x-1) & onset<=x
                    split_info.cond(i).onsetlist =  ceil(split_info.cond(i).onsetlist/split_info.TR_MSEC);
                    % round blocklengths
                    split_info.cond(i).blklength = round(split_info.cond(i).blklength/split_info.TR_MSEC); %
                end
                
                % quick check to ensure blocks are >=1 TR in length
                for i = 1:size(split_info.Contrast_List,1) % each contrast 
                    for condition = 1:2 % two conditions
                        if split_info.Contrast_List(i,condition)~=0 % only check non-baseline for now
                            if    ( sum(split_info.cond(split_info.Contrast_List(i,condition)).blklength)==0 )
                                error(['cannot run block-design: condition "',split_info.cond(i).name,'" has no stimulus durations >1TR']);
                            elseif( sum(split_info.cond(split_info.Contrast_List(i,condition)).blklength==0)>0 )
                                warning(['condition "',split_info.cond(i).name,'" has some stimulus blocks of duration <1TR (will not be used!']);
                            end
                        end
                    end
                end                      
                
                % if type=block automatically split the data in half/           
                nrun       = numel( InputStruct(ksub).run ); %% number of runs for this subject
                split_info = split_half(split_info,InputStruct(ksub).run(krun).Nt,filename, nrun); % create single and group split halves                
            
            %% for event design, zero-out durations --> these are ignored anyways during windowed analysis
            %
            elseif strcmpi(analysis_model.design_type,'event')
                % for all further analyses, zero out the durations
                for ktemp = 1:length(split_info.cond)
                    split_info.cond(ktemp).blklength = zeros(size(split_info.cond(ktemp).onsetlist));
                end            
            end
        else
            % case with no contrasts specified
            split_info.design_mat = [];    %% empty design matrix placeholder
            split_info.Contrast_List = 0;  %% zero-field if no contrasts specified
        end
        % add signle and group field to be used for signle subject or group
        % analysis, for the group analysis we do not need splitting, so it
        % merges the splitting content.
        
        % now storing parsed split_info into input structure
        InputStruct(ksub).run(krun).split_info = split_info;
        
        % store contrast lists  -- these should be identical across runs
        % since (a) we only keep conditions in contrasts, and (b) they are alphabetized
        Contrast_set_temp{krun} = split_info.Contrast_List;
    end
    
    % sanity check on contrasts -- all pairwise comparisons
    for krun1 = 1:numel(InputStruct(ksub).run)-1
    for krun2 = krun1+1:numel(InputStruct(ksub).run)
        
        if( numel( Contrast_set_temp{krun1} ) ~= numel( Contrast_set_temp{krun2} ) )
            error(['number of contrasts not consistent across runs? For ',InputStruct(ksub).run(krun).split_info_file]);
        elseif( sum( abs( Contrast_set_temp{krun1}(:)-Contrast_set_temp{krun2}(:) ) )>0 )
            error(['specific contrasts not consistent across runs? For ',InputStruct(ksub).run(krun).split_info_file]);
        end        
    end
    end
    clear Contrast_set_temp;
end

%% extract contrast list
function split_info = extract_contrast_list(split_info,contrast_list_str) % create single and group split halves

% if at least one explicit contrast (difference) is defined...
if ~isempty(contrast_list_str) && ~strcmpi(contrast_list_str,'NONE') && isfield(split_info,'cond')
    
    % split up list of contrasts
    ind=find(contrast_list_str==',');
    contrast_list_str(ind) = ' ';
    if(isempty(strfind(contrast_list_str,' ')))
        Contrast_temp = cellstr(contrast_list_str);
    else Contrast_temp = regexp(contrast_list_str,' ','split');
    end
    num_contrast = length(Contrast_temp);
    
    % for each contrast, format for analysis
    for k = 1:num_contrast
        
        if(isempty(strfind(Contrast_temp{k},'-')))
             Contrast_name = {Contrast_temp{k},'baseline'}; %% if not a difference, assume cf. relative to baseline
        else Contrast_name = regexp(Contrast_temp{k},'-','split'); %% if a difference, split into 2 cells
        end
        
        for n = 1:2
            Contrast_name{n} = strrep(Contrast_name{n},' ',''); % remove extra space
        end
        
        % now convert to kth contrast pair
        Contrast_list{k,1} = Contrast_name{1};
        Contrast_list{k,2} = Contrast_name{2};
        
        % integrate into split_info file --> creates new "combined"
        % conditions (with onsets & durations) if specifed by user
        split_info = interpret_combined_contrast(split_info,Contrast_list{k,1});
        split_info = interpret_combined_contrast(split_info,Contrast_list{k,2});
        
        for i = 1:length(split_info.cond)
            Condition_name{i} = split_info.cond(i).name;
        end
        
        % check which named condition appears in the contrast and the associated index                
        % if it matches "baseline" reset index to zero!
        if strcmpi(Contrast_list{k,1},'baseline')
              ind1 = 0;
        else  ind1 = find(ismember(lower(Condition_name),lower(Contrast_list{k,1})));
        end
        if strcmpi(Contrast_list{k,2},'baseline')
              ind2 = 0;
        else  ind2 = find(ismember(lower(Condition_name),lower(Contrast_list{k,2})));            
        end
        
        % check if any conditions can't be matched
        if( isempty(ind1) || isempty(ind2) )
        display(sprintf('ERROR: cannot find element of contrast (%s) among conditions:\n',Contrast_temp{k}));
        display(Condition_name);    
        sge_exit(100);    
        end
        
        % now numerically index the contrast
        Contrast(k,:) =[ind1 ind2];
    end
    
    split_info.Contrast_List = Contrast;
    
else
    split_info.Contrast_List = 0; %% zero-field if no contrasts are defined
end 

%% create "combined" contrasts
function split_out = interpret_combined_contrast(split_info,Cond_in)

split_out = split_info;
if strfind(Cond_in,'+') %% this indicates >1 conditions are to be combined
    
    duplic=0; %% increment if a duplicate is found
    for(i=1:length(split_out.cond))
        duplic = duplic + double( strcmpi(Cond_in,split_out.cond(i).name) );
    end
    
    if( duplic ==0 )
    
        Cond_temp = Cond_in;
        Cond_temp(Cond_temp=='+') = '%'; % For Octave Compatibility, Octave is sensitive to +

        if(isempty(strfind(Cond_temp,'%')))
            Cond_name = cellstr(Cond_temp);
        else Cond_name = regexp(Cond_temp,'%','split');
        end
        % collect ordered array of condition names
        for i = 1:length(split_info.cond)
            Condition_name{i} = split_info.cond(i).name;
        end
        onsetlist = [];
        blklength = [];
        % go through split_info list --> find associated onsets / concatenate them
        for i = 1:length(Cond_name)
            ind       = find(ismember(lower(Condition_name),lower(Cond_name{i})));
            if isempty(ind)
                ind = str2num(Cond_name{i});
            end
            onsetlist = [onsetlist;split_info.cond(ind).onsetlist(:)];
            blklength = [blklength;split_info.cond(ind).blklength(:)];
        end
        % create "new" condition that is combination of individual ones
        n = length(split_out.cond);
        split_out.cond(n+1).onsetlist = onsetlist';
        split_out.cond(n+1).blklength = blklength';
        split_out.cond(n+1).name      = lower(Cond_in);
    
    end
end

%% automatic splitting of blocks
function split_info = split_half(split_info,Nt,filename, nrun)

%collect list of contrasts
Contrast_List = split_info.Contrast_List;

if split_info.baseline_index==0 %% no baseline specified
    all_stim = [];
    for i = 1:length(split_info.cond)
        for k = 1:length(split_info.cond(i).onsetlist)
            all_stim = [all_stim split_info.cond(i).onsetlist(k):split_info.cond(i).onsetlist(k)+split_info.cond(i).blklength(k)-1];
        end
    end
    split_info.baseline = sort(setdiff(1:Nt,all_stim));   % whatever left in time-series is considered as baseline
    
else %% if explicit baseline defined
    all_stim = [];
    for i = split_info.baseline_index
        for k = 1:length(split_info.cond(i).onsetlist)
            all_stim = [all_stim split_info.cond(i).onsetlist(k):split_info.cond(i).onsetlist(k)+split_info.cond(i).blklength(k)-1];
        end
    end
    split_info.baseline = all_stim;   % whatever left in time-series is considered as baseline
end

for i = 1:size(Contrast_List,1) %% iterate through list of contrasts to analyze
    for condition = 1:2
        split_info.group.idx_cond(i,condition).sp = [];
        if Contrast_List(i,condition)~=0 %% if non-baseline
            for k = 1:length(split_info.cond(Contrast_List(i,condition)).onsetlist)
                split_info.group.idx_cond(i,condition).sp = [split_info.group.idx_cond(i,condition).sp split_info.cond(Contrast_List(i,condition)).onsetlist(k):split_info.cond(Contrast_List(i,condition)).onsetlist(k)+split_info.cond(Contrast_List(i,condition)).blklength(k)-1];
            end
        else %% if baseline
            split_info.group.idx_cond(i,condition).sp = [split_info.group.idx_cond(i,condition).sp split_info.baseline];
        end
        % make sure condition time points are ordered numerically
        split_info.group.idx_cond(i,condition).sp = sort(split_info.group.idx_cond(i,condition).sp);
    end
end

% now splitting the list of timepoints into 2 split-halves
for i = 1:size(Contrast_List,1)
    
    max_sp = ceil(Nt/2); %% allocated relative to middle timepoint
    
    % First split
    class1 = split_info.group.idx_cond(i,1).sp(split_info.group.idx_cond(i,1).sp<=max_sp);
    class2 = split_info.group.idx_cond(i,2).sp(split_info.group.idx_cond(i,2).sp<=max_sp);
    % catch bad splitting
    if(numel(class1)<2)
        if(nrun==1)
             sge_exit(100, sprintf(         'file %s, cond. %s: %u scans in 1st half of run -- cannot do single-run analysis.',filename,split_info.cond(Contrast_List(i,1)).name,numel(class1)));
        else          disp(sprintf('WARNING: file %s, cond. %s: %u scans in 1st half of run -- only multi-run will work.',filename,split_info.cond(Contrast_List(i,1)).name,numel(class1)));
        end
    end
    if(numel(class2)<2)
        if(nrun==1)
             sge_exit(100, sprintf(         'file %s, cond. %s: %u scans in 1st half of run -- cannot do single-run analysis.',filename,split_info.cond(Contrast_List(i,2)).name,numel(class2)));
        else          disp(sprintf('WARNING: file %s, cond. %s: %u scans in 1st half of run -- only multi-run will work.',filename,split_info.cond(Contrast_List(i,2)).name,numel(class2)));
        end        
    end
    %
    [class1,class2] = make_two_class_similar_length(class1,class2);
    split_info.single.idx_cond(i,1).sp1 = class1;
    split_info.single.idx_cond(i,2).sp1 = class2;
    
    % Second split
    class1 = split_info.group.idx_cond(i,1).sp(split_info.group.idx_cond(i,1).sp>max_sp);
    class2 = split_info.group.idx_cond(i,2).sp(split_info.group.idx_cond(i,2).sp>max_sp);
    % catch bad splitting
    if(numel(class1)<2)
        if(nrun==1)
             sge_exit(100, sprintf(         'file %s, cond. %s: %u scans in 2nd half of run -- cannot do single-run analysis.',filename,split_info.cond(Contrast_List(i,1)).name,numel(class1)));
        else          disp(sprintf('WARNING: file %s, cond. %s: %u scans in 2nd half of run -- only multi-run will work.',filename,split_info.cond(Contrast_List(i,1)).name,numel(class1)));
        end
    end
    if(numel(class2)<2)
        if(nrun==1)
             sge_exit(100, sprintf(         'file %s, cond. %s: %u scans in 2nd half of run -- cannot do single-run analysis.',filename,split_info.cond(Contrast_List(i,2)).name,numel(class2)));
        else          disp(sprintf('WARNING: file %s, cond. %s: %u scans in 2nd half of run -- only multi-run will work.',filename,split_info.cond(Contrast_List(i,2)).name,numel(class2)));
        end        
    end
    
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

%%
function [class_out1,class_out2] = make_two_class_similar_length(class1,class2)

class_out1 = class1;
class_out2 = class2;
% temporal distance between all paired points (compared between splits)
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

function [Xsub,idx]=licols(X,tol)
%Extract a linearly independent set of columns of a given matrix X
%
%    [Xsub,idx]=licols(X)
%
%in:
%
%  X: The given input matrix
%  tol: A rank estimation tolerance. Default=1e-10
%
%out:
%
% Xsub: The extracted columns of X
% idx:  The indices (into X) of the extracted columns
if ~nnz(X) %X has no non-zeros and hence no independent columns
    Xsub=[]; idx=[];
    return
end
if nargin<2, tol=1e-8; end
[Q, R, E] = qr(X,0);
if ~isvector(R)
    diagr = abs(diag(R));
else
    diagr = R(1);
end
%Rank estimation
r = find(diagr >= tol*diagr(1), 1, 'last');
idx=sort(E(1:r));
Xsub=X(:,idx);