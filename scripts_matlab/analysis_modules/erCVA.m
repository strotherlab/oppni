function output = erCVA( datamat, split_info, Resampling_Index )
%
% =========================================================================
% MODULE_ERCVA: module that performs hybrid Canonical Variates Analysis
% with averaging for event-related data
% =========================================================================
%
%   Syntax:
%           output = module_erCVA( datamat, split_info )
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
output.attributes.model_name  = 'erCVA';
output.attributes.design_type = 'event';
output.attributes.model_type  = 'multivariate';
output.attributes.num_comp    = 'one_component';

if(nargin==0)
    disp('no inputs - returning attributes');
    return;
end
%----------------------- Default Parameter Checks ------------------------%
if( ~isfield(split_info{1},'WIND')   || isempty(split_info{1}.WIND)   ) 
    disp('erCVA uses default window-size WIND=10');
    split_info{1}.WIND =10; 
end
if( ~isfield(split_info{1},'Nblock') || isempty(split_info{1}.Nblock) ) 
    disp('erCVA uses default number of blocks Nblock=4');
    split_info{1}.Nblock=4; 
end
if( ~isfield(split_info{1},'norm')   || isempty(split_info{1}.norm)   ) 
    disp('erCVA uses default normalization turned on (norm=1)');
    split_info{1}.norm  =1; 
end
if( ~isfield(split_info{1},'drf') || isempty(split_info{1}.drf) )
    disp('erCVA uses default data reduction drf=0.5');
    split_info{1}.drf = 0.5;
end
%-------------------------------------------------------------------------%
 
% initialized parameters
params.TR     = split_info{1}.TR_MSEC;
params.delay  = split_info{1}.TR_MSEC./2;
params.WIND   = split_info{1}.WIND;
params.norm   = 1;

%% INDIVIDUAL SUBJECT ANALYSIS
if( ~iscell(datamat) || length(datamat)==1 )
    
    % in cases where it is a single cell, revert to single-subject anaylsis
    if( iscell(datamat) && length(datamat)==1 ) datamat = datamat{1}; end
    split_info = split_info{1}; %% take entry from single cell

    % % % count #leftovers, if we do even multiples of splits
    % % leftover = rem( size(datamat,2) , params.Nblock);
    % % % trim "overhang" from the data matrix
    % % datamat_trim  = datamat(:,1:end-leftover);

    %  Time-windowed averaging of data:
    % *NB: simple averaging does nearest-neighbour interpolation to nearest TR
    %      interval. You can replace with "interp_averaging_for_ER" to do
    %      linear time interpolation; the tends to increase reproducibility,
    %      but at the apparent expense of decreased prediction!
    %
    
    % blocks pre-specified
    params.Nblock = split_info.Nblock;

    [windowAvg] = simple_averaging_for_ER( datamat, split_info.onsetlist, params );

    % --- CVAanalysis --- %
    fullPC = floor(split_info.drf*split_info.Nblock*split_info.WIND);
    spltPC = round( fullPC / 2 );

    [outer] = ercva_optimization( windowAvg, [fullPC spltPC],'spat');

    % --- Subspace selection --- %

    % Get HRF reference if requested
    %
    if( ~isfield( split_info, 'lag_match' ) ) 
        %
        split_info.lag_match = 0; 
        %
    elseif( split_info.lag_match > 0 )
        %
        inputv    = zeros( params.WIND, 1 );
        inputv(1) = 1;
        hrf_ref   = design_to_hrf( inputv, params.TR/1000, [(split_info.lag_match/1000) 15.0] );
    end

    % OPTION 1: single component analysis; estimates the primary HRF.
    %           based on spatial reproducibility and single-component prediction
    %        * If user does not specify optimization approach, single-component is default
    if   ( ~isfield( split_info, 'subspace' ) || strcmp(split_info.subspace, 'onecomp') )

            % iterate through pc dimensions
            for(q=1:spltPC)
                %
                if( split_info.lag_match == 0 )
                    %
                    tmpR(q,1)   = outer.R{q}(1);        
                    tmpP(q,1)   = outer.P_global(q); %% outer.P_per_cv{q}(1);
                    tmpSPM(:,q) = outer.SPMs{q}(:,1);
                    tmpTset(:,q)= outer.Tset{q}(:,1); 
                else
                    %
                    cc = corr( outer.Tset{q}, hrf_ref );
                    [vc ic] = max( abs(cc) ); 
                    %
                    tmpR(q,1)   = outer.R{q}(ic);        
                    tmpP(q,1)   = outer.P_global(q); %%outer.P_per_cv{q}(ic);
                    tmpSPM(:,q) = outer.SPMs{q}(:,ic) .* sign( cc(ic) );
                    tmpTset(:,q)= outer.Tset{q}(:,ic) .* sign( cc(ic) );                 
                end
            end

            tmpD      = sqrt( (1-tmpR).^2 + (1-tmpP).^2 );
            [vd id]   = min(tmpD);
            %
            output.images  = tmpSPM(:,id);
            output.temp    = tmpTset(:,id);
            %
            output.metrics.R    =  tmpR(id);
            output.metrics.P    =  tmpP(id);
            output.metrics.dPR  = -vd;

    % OPTION 2: multi component subspace analysis; estimated based on 
    %           spatial reproducibility of stable SPM bases, 
    %           and global prediction accuracy
    elseif( strcmp(split_info.subspace, 'multicomp') )
            %
            % first-pass:
            % (a) discard "unstable", non reproducible components (no significant Z-scores at FDR=.10)
            for(q=1:spltPC) 
                [p thr] = fdr_ex( outer.SPMs{q}, 0.05, 0 );
                thrsum  = sum(thr);

                if( sum(thrsum>0) > 0 )
                %
                outer.SPMs{q} = outer.SPMs{q}(:,thrsum>0);
                outer.Tset{q} = outer.Tset{q}(:,thrsum>0);
                outer.R{q}    = outer.R{q}(thrsum>0);  
                %
                % NOTA BENE: two ways of quantifying overall reproducibility
                % one of these needs to be commented out -- your choice!
                %
                % (1) average reproducibility of the stable components
                %%tmpR(q,1) = mean( outer.R{q} );
                % (2) fraction of W-1 dimensional subspace variance (R^2) that is reproducible
                tmpR(q,1) =  sum( outer.R{q}.^2 ) ./ (split_info.WIND-1);            
                else
                tmpR(q,1) = -1;
                end
            end

            tmpP    = outer.P_global;
            tmpD    = sqrt( (1-tmpR).^2 + (1-tmpP).^2 );
            [vd id] = min(tmpD);
            %
            output.images  = outer.SPMs{id};
            output.temp    = outer.Tset{id};
            %
            output.metrics.R    =  outer.R{id};
            output.metrics.P    =  tmpP(id);
            output.metrics.dPR  = -vd;
    else
        error('invalid component selection method');
    end
    
%% GROUP LEVEL ANALYSIS
else

    % number of subjects
    N_subject = length(datamat);
    Nvox      = size(datamat{1},1);

    if( N_subject >=4 ) %% general case -- runs = resampling units

        % Enforce 1 split per subject for group analysis
        params.Nblock = 1; 

        % initialize data matrix
        windowAvg = zeros( Nvox, split_info{1}.WIND, N_subject );
        % load from averaged blocks from cell arrays
        for( n=1:N_subject )
            % perform stimulus-locked averaging on each subject
            [windowAvg(:,:,n)] = simple_averaging_for_ER( datamat{n}, split_info{n}.onsetlist, params );
        end

        % --- CVA analysis --- %
        %
        fullPC = floor(split_info{1}.drf*N_subject*split_info{1}.WIND);
        spltPC = round( fullPC / 2 );
        %
        [outer] = ercva_optimization( windowAvg, [fullPC spltPC],'temp');

    else %% cases where only 2-3 runs ... need enough for variance estimation within each split (2 samples per class)
        
        % Enforce 2 splits per subject for group analysis
        params.Nblock = 2; 

        % initialize data matrix
        windowAvg = zeros( Nvox, split_info{1}.WIND, 2*N_subject );
        % load from averaged blocks from cell arrays
        for( n=1:N_subject )
            % perform stimulus-locked averaging on each subject
            [windowAvg(:,:,2*n-1:2*n)] = simple_averaging_for_ER( datamat{n}, split_info{n}.onsetlist, params );
        end

        % --- CVA analysis --- %
        %
        fullPC = floor(split_info{1}.drf*2*N_subject*split_info{1}.WIND);
        spltPC = round( fullPC / 2 );
        %
        [outer] = ercva_optimization( windowAvg, [fullPC spltPC],'temp');
    end
    
    
    % --- Subspace selection --- %

    % Get HRF reference if requested
    %
    if( ~isfield( split_info{1}, 'lag_match' ) ) 
        %
        split_info{1}.lag_match = 0; 
        %
    elseif( split_info{1}.lag_match > 0 )
        %
        inputv    = zeros( params.WIND, 1 );
        inputv(1) = 1;
        hrf_ref   = design_to_hrf( inputv, params.TR/1000, [(split_info{1}.lag_match/1000) 15.0] );
    end

    % OPTION 1: single component analysis; estimates the primary HRF.
    %           based on spatial reproducibility and single-component prediction
    %        * If user does not specify optimization approach, single-component is default
    if   ( ~isfield( split_info{1}, 'subspace' ) || strcmp(split_info{1}.subspace, 'onecomp') ) 

            % iterate through pc dimensions
            for(q=1:spltPC)
                %
                if( split_info{1}.lag_match == 0 )
                    %
                    tmpR(q,1)   = outer.R{q}(1);        
                    tmpP(q,1)   = outer.P_global(q); %%outer.P_per_cv{q}(1);
                    tmpSPM(:,q) = outer.SPMs{q}(:,1);
                    tmpTset(:,q)= outer.Tset{q}(:,1); 
                else
                    %
                    cc = corr( outer.Tset{q}, hrf_ref );
                    [vc ic] = max( abs(cc) ); 
                    %
                    tmpR(q,1)   = outer.R{q}(ic);        
                    tmpP(q,1)   = outer.P_global(q); %%outer.P_per_cv{q}(ic);
                    tmpSPM(:,q) = outer.SPMs{q}(:,ic) .* sign( cc(ic) );
                    tmpTset(:,q)= outer.Tset{q}(:,ic) .* sign( cc(ic) );                 
                end
            end

            tmpD      = sqrt( (1-tmpR).^2 + (1-tmpP).^2 );
            [vd id]   = min(tmpD);
            %
            output.images  = tmpSPM(:,id);
            output.temp    = tmpTset(:,id);
            %
            output.metrics.R    =  tmpR(id);
            output.metrics.P    =  tmpP(id);
            output.metrics.dPR  = -vd;

    % OPTION 2: multi component subspace analysis; estimated based on 
    %           spatial reproducibility of stable SPM bases, 
    %           and global prediction accuracy
    elseif( strcmp(split_info{1}.subspace, 'multicomp') )
            %
            % first-pass:
            % (a) discard "unstable", non reproducible components (no significant Z-scores at FDR=.10)
            for(q=1:spltPC) 
                [p thr] = fdr_ex( outer.SPMs{q}, 0.05, 0 );
                thrsum  = sum(thr);

                if( sum(thrsum>0) > 0 )
                %
                outer.SPMs{q} = outer.SPMs{q}(:,thrsum>0);
                outer.Tset{q} = outer.Tset{q}(:,thrsum>0);
                outer.R{q}    = outer.R{q}(thrsum>0);  
                %
                % NOTA BENE: two ways of quantifying overall reproducibility
                % one of these needs to be commented out -- your choice!
                %
                % (1) average reproducibility of the stable components
                %%tmpR(q,1) = mean( outer.R{q} );
                % (2) fraction of W-1 dimensional subspace variance (R^2) that is reproducible
                tmpR(q,1) =  sum( outer.R{q}.^2 ) ./ (split_info{1}.WIND-1);            
                else
                tmpR(q,1) = -1;
                end
            end

            tmpP    = outer.P_global;
            tmpD    = sqrt( (1-tmpR).^2 + (1-tmpP).^2 );
            [vd id] = min(tmpD);
            %
            output.images  = outer.SPMs{id};
            output.temp    = outer.Tset{id};
            %
            output.metrics.R    =  outer.R{id};
            output.metrics.P    =  tmpP(id);
            output.metrics.dPR  = -vd;
    else
        error('invalid component selection method');
    end    
end

%%
function [pcritSet threshMat] = fdr_ex( dataMat, qval, cv )
% 
% trimmed-down FDR code
%

[Ntest Nk] = size(dataMat);
probMat    = 1-normcdf( abs(dataMat) );    
threshMat  = zeros( Ntest, Nk );

for( K=1:Nk )

    % (1) order pvalues smallest->largest
    pvect = sort( probMat(isfinite( probMat(:,K) ), K), 'ascend' );
    Ntest2= length(pvect);
    % (2) find highest index meeting limit criterion
    if(cv == 0) c_V = 1;
    else        c_V = log(Ntest2) + 0.5772 + 1/Ntest2; % approx soln.
    end
    % index vector
    indxd = (1:Ntest2)';
    % get highest index under adaptive threshold
    r = sum( pvect./indxd <= qval/(Ntest2*c_V) );

    if( r > 0 )
        %
        % limiting p-value
        pcrit = pvect(r);
        % threshold matrix values
        threshMat(:,K) = double( probMat(:,K) <= pcrit );
        pcritSet(K,1)  = pcrit;
    else
        threshMat(:,K) = zeros(Ntest,1);
        pcritSet(K,1)  = NaN;
    end
end

