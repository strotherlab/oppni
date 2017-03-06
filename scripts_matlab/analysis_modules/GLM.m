function output = GLM( datamat, split_info, Resampling_Index )
%
% =========================================================================
% MODULE_GLM: module that performs General Linear Model analyses
% =========================================================================
%
%   Syntax:
%           output = module_erGLM( datamat, split_info )
%
%
% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name  = 'GLM';
output.attributes.design_type = 'regression';
output.attributes.model_type  = 'univariate';
output.attributes.num_comp    = 'multi_component';

if(nargin==0)
    disp('no inputs - returning attributes');
    return;
end
%----------------------- Default Parameter Checks ------------------------%
%-------------------------------------------------------------------------%

% ==== start: check for combined contrasts and index ====
%
% collect ordered array of condition names
for i = 1:length(split_info{1}.cond)
    Condition_name{i} = split_info{1}.cond(i).name;
end
n_combin=0;
n_noncomb=size(split_info{1}.design_mat,2); %number of "uncombined" conditions
% go through list of condition names
for( i=1:numel( Condition_name ) )
    %
    % check if "combined"
    if( ~isempty(strfind(Condition_name{i},'+')) )
        %
        if( i<= n_noncomb ) error('a "+" found in a non-combined condition! Check your syntax!'); end
        % split into separate conditions
        Cond_split = regexp(Condition_name{i},'+','split');
        % now assign indices
        for j = 1:length(Cond_split)
            % cell array - condition x (sub-condition)
            index_combin{i}(j)       = find(ismember(lower(Condition_name),lower(Cond_split{j})));
        end
        n_combin = n_combin+1; %% increment #combined contrasts
    else
        index_combin{i}=[];
    end
end
%
% ==== end -- check for combined contrasts and index ====

disp('now running glm...');

%% INDIVIDUAL SUBJECT ANALYSIS
if( ~iscell(datamat) || numel(datamat)==1 )
    
    % in cases where it is a single cell, revert to single-subject anaylsis
    if( iscell(datamat) && length(datamat)==1 ) datamat = datamat{1}; end
    split_info = split_info{1}; %% take entry from single cell

    % matrix dimensions
    [Nvox Ntime] = size( datamat );

    X_design = split_info.design_mat;

    % split the data matrix
    datasplit_1 = datamat(:,1:ceil(Ntime/2));
    datasplit_2 = datamat(:,ceil(Ntime/2)+1:end);
    % split the hrf-convolved design vector
    X_design_1 = X_design(1:ceil(Ntime/2),:);
    X_design_2 = X_design(ceil(Ntime/2)+1:end,:);
    % run split glm analyses
    out1 = GLM_model_fmri( datasplit_1, 0, [], X_design_1, [] );
    out2 = GLM_model_fmri( datasplit_2, 0, [], X_design_2, [] );

    % ==== now create combined-contrast (+) betas, if included
    if( n_combin>0 )
       for( i=1:numel( Condition_name ) )

           if( ~isempty(index_combin{i}) )
               % sum across indices
               out1.Beta_signl(:,i) = mean( out1.Beta_signl(:,index_combin{i}),2 );
               out2.Beta_signl(:,i) = mean( out2.Beta_signl(:,index_combin{i}),2 );                  
           end
       end
    end
    % ==== now create combined-contrast (+) betas, if included
    
    % now store specific contrasts of interest
    for(k=1:size(split_info.Contrast_List,1))
        contr = split_info.Contrast_List(k,:);
        if( min(contr)==0 ) % relative to baseline
        beta_contr1(:,k) = out1.Beta_signl(:,max(contr));
        beta_contr2(:,k) = out2.Beta_signl(:,max(contr));
        else % otherwise
        beta_contr1(:,k) = out1.Beta_signl(:,contr(1)) - out1.Beta_signl(:,contr(2));
        beta_contr2(:,k) = out2.Beta_signl(:,contr(1)) - out2.Beta_signl(:,contr(2));
        end
    end
    
    % reweight the betas
    beta_contr1 = bsxfun(@times,beta_contr1,split_info.spat_weight);
    beta_contr2 = bsxfun(@times,beta_contr2,split_info.spat_weight);

    % reproducibility and SPM:
    [ RR, rSPMZ ] = get_rSPM( beta_contr1, beta_contr2, 1 );

    % reproducibility
    output.metrics.R = RR;
    % optimal eigenimage
    output.images  = rSPMZ;
    % CV score timeseries, on unit-normed eigenimage
    output.temp    = datamat'  * bsxfun(@rdivide,rSPMZ,sqrt(sum(rSPMZ.^2)));

%% GROUP LEVEL ANALYSIS
else

    N_subject = numel(datamat);

    for k = 1:N_subject
        X_design{k} = split_info{k}.design_mat;
    end
    N_resample = size(Resampling_Index,1);
    
    rSPMZ=0; %% initialize rSPMz
    
    for i = 1:N_resample
        
        % split the data matrix
        set1 = Resampling_Index(i,:);
        set2 = setdiff(1:N_subject,set1);

        datasplit1 = [];
        X_design1  = [];
        for k = 1:length(set1)
            datasplit1 = [datasplit1 datamat{set1(k)}];
            X_design1  = [X_design1;X_design{set1(k)}];
        end

        datasplit2 = [];
        X_design2  = [];
        for k = 1:length(set2)
            datasplit2 = [datasplit2 datamat{set2(k)}];
            X_design2  = [X_design2;X_design{set2(k)}];
        end

        % run split glm analyses
        out1 = GLM_model_fmri( datasplit1, 0, [], X_design1, [] );
        out2 = GLM_model_fmri( datasplit2, 0, [], X_design2, [] );

        % ==== now create combined-contrast (+) betas, if included
        if( n_combin>0 )
           for( i=1:numel( Condition_name ) )
              
               if( ~isempty(index_combin{i}) )
                   % sum across indices
                   out1.Beta_signl(:,i) = mean( out1.Beta_signl(:,index_combin{i}),2 );
                   out2.Beta_signl(:,i) = mean( out2.Beta_signl(:,index_combin{i}),2 );                  
               end
           end
        end
        % ==== now create combined-contrast (+) betas, if included
        
        % now store specific contrasts of interest
        for(k=1:size(split_info{1}.Contrast_List,1))
            contr = split_info{1}.Contrast_List(k,:);
            if( min(contr)==0 ) % relative to baseline
            beta_contr1(:,k) = out1.Beta_signl(:,max(contr));
            beta_contr2(:,k) = out2.Beta_signl(:,max(contr));
            else % otherwise
            beta_contr1(:,k) = out1.Beta_signl(:,contr(1)) - out1.Beta_signl(:,contr(2));
            beta_contr2(:,k) = out2.Beta_signl(:,contr(1)) - out2.Beta_signl(:,contr(2));
            end
        end

        % reweight the betas
        beta_contr1 = bsxfun(@times,beta_contr1,split_info{1}.spat_weight);
        beta_contr2 = bsxfun(@times,beta_contr2,split_info{1}.spat_weight);

        % reproducibility and SPM:
        [ RR_temp, rSPMZ_temp ] = get_rSPM( beta_contr1, beta_contr2, 1 );

        rSPMZ = rSPMZ + rSPMZ_temp/N_resample;
        RR_stack(i,:) = RR_temp;
    end

    % reproducibility
    output.metrics.R = median(RR_stack,1);   
    % optimal eigenimage
    output.images    = rSPMZ;
    
    for(k=1:N_subject)
    % CV score timeseries, on unit-normed eigenimage
    output.temp{k} = datamat{k}'  * bsxfun(@rdivide, rSPMZ , sqrt(sum(rSPMZ.^2)) );
    end
end
