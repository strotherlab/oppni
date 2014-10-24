function output = module_erCVA( datamat, split_info )
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

% initialized parameters
params.TR     = split_info.TR_MSEC;
params.delay  = split_info.TR_MSEC./2;
params.WIND   = split_info.WIND;
params.Nsplit = split_info.Nsplit;
params.norm   = 1;

% % % count #leftovers, if we do even multiples of splits
% % leftover = rem( size(datamat,2) , params.Nsplit);
% % % trim "overhang" from the data matrix
% % datamat_trim  = datamat(:,1:end-leftover);

%  Time-windowed averaging of data:
% *NB: simple averaging does nearest-neighbour interpolation to nearest TR
%      interval. You can replace with "interp_averaging_for_ER" to do
%      linear time interpolation; the tends to increase reproducibility,
%      but at the apparent expense of decreased prediction!
%
[windowAvg] = simple_averaging_for_ER( datamat, split_info.onsetlist, params );

% --- CVAanalysis --- %
fullPC = floor(split_info.drf*split_info.Nsplit*split_info.WIND);
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
        output.metrics.Dneg = -vd;
        
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
        output.metrics.Dneg = -vd;
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

