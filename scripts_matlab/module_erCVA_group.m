function output = module_erCVA_group( datamat, split_info )
%
% =========================================================================
% MODULE_ERCVA_GROUP: module that performs group-level hybrid Canonical 
% Variates Analysis with averaging for event-related data
% =========================================================================
%
%   Syntax:
%           output = module_erCVA_group( datamat, split_info )
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
    split_info.drf = 0.5;
end

% initialized parameters
params.TR     = split_info{1}.TR_MSEC;
params.delay  = split_info{1}.TR_MSEC ./ 2;
params.WIND   = split_info{1}.WIND;
params.Nblock = 1; %% now, only 1 split per subject!
params.norm   = 1;
% number of subjects
N_subject = length(datamat);
Nvox      = size(datamat{1},1);

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
