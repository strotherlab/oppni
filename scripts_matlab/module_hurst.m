function output = module_hurst( datamat, split_info )
%
% =========================================================================
% MODULE_HURST: computes hurst exponent on fmri timeseries
% =========================================================================
%
%   Syntax:
%           output = module_hurst( datamat, split_info )
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

% dimensions
[Nvox Ntime] =  size(datamat); 
% split into two equal halves
datamat_sp1 = datamat(:,1:ceil(Ntime/2));
datamat_sp2 = datamat(:,ceil(Ntime/2)+1:end);
% mean-center splits
datamat_sp1 = bsxfun(@minus,datamat_sp1,mean(datamat_sp1,2));
datamat_sp2 = bsxfun(@minus,datamat_sp2,mean(datamat_sp2,2));

if( (isfield(split_info,'hurst_units') && ~isempty(split_info.hurst_units)) && (isfield(split_info,'hurst_range') && ~isempty(split_info.hurst_range)) )
    %
    % if both units and range correctly specified
    [H_vect1] = DFA_parallel( datamat_sp1, split_info.hurst_units, split_info.hurst_range );
    [H_vect2] = DFA_parallel( datamat_sp2, split_info.hurst_units, split_info.hurst_range );
else
    % otherwise default
    [H_vect1] = DFA_parallel( datamat_sp1 );
    [H_vect2] = DFA_parallel( datamat_sp2 );
end
% catch for bad-valued voxels
H_vect1(~isfinite(H_vect1))=0.5+eps;
H_vect2(~isfinite(H_vect2))=0.5+eps;

% reproducibility estimates on "most independent" splits
[output.metrics.R rSPMZ] = get_rSPM( H_vect1, H_vect2, 1 );

if    ( strcmp( split_info.spm, 'zscore' ) ) output.images  = rSPMZ;
elseif( strcmp( split_info.spm, 'raw'   ) )  output.images  = (H_vect1+H_vect2)./2;
else  error( 'invalid output type for hurst (should be zscore or raw).' );
end

output.temp    = datamat'  * (rSPMZ ./ sqrt(sum(rSPMZ.^2)));

%%
function [H_vect log_FVAR scale_list] = DFA_parallel( DataMat, varargin )
%
% -------------------------------------------------------------------------
% DETRENDED FLUCTUATIONS ANALYSIS: very simple, (relatively) computationally
% efficient code for computing the Hurst exponent simultaneously for a
% matrix of time-series vectors. This measures the long-range temporal
% dependency of your signal estimates, and is robust against linear signal
% drift effects.
% -------------------------------------------------------------------------
%
% Syntax:
%            [H_vect log_FVAR scale_list] = DFA_parallel( DataMat, time_units, scal_range )
%
% Input:
%          DataMat   : a 2D matrix arranged (sample x time), e.g. a row matrix of timeseries.
%          time_units: (optional) integeter specifying the sampling rate
%          scal_range: (optional) 2D vector denoting the [minimum maximum]
%                      timescale interval on which to run DFA. Need to
%                      specify time_units (if time_units=[]), scal_range is
%                      assumed to be in "sampling units"
%
% Output:    
%          H_vect    : (sample x 1) vector of Hurst coefficient estimates,
%                      quantifything fractal scaling for the set of samples.
%          log_FVAR  : (sample x timescale) matrix of fluctuation log-variance 
%                      as a function of time-scale
%          scale_list: (timescale x 1) vector of the corresponding
%                      time-scales (i.e. window sizes). Either in sampling units,
%                      or "time_units" increments, if specified in input
%
%

if(isempty(varargin))
    time_units = [];
    scal_range = [];
elseif( numel(varargin)==1 )
    time_units = varargin{1};
    scal_range = [];
elseif( numel(varargin)==2 )
    time_units = varargin{1};
    scal_range = varargin{2};
    
    if( length(scal_range)~=2 )
        error('scal_range must be a 2D vector');
    end
else
    error('number of arguments exceeds possible inputs'); 
end
    

% matrix dmensions
[Nvox Nsamp] = size( DataMat );

% STEP 0: linear+constant detrending on full data matrix, to minimize
%         any obvious artifact due to scanner drift, etc.
% const+linear regressors
D01 = [ ones(Nsamp,1), linspace(-1,1,Nsamp)' ];
% compute beta (Regression coefficients)
BetaFits = DataMat * D01 * inv( D01' * D01 ); 
% subtract estimated const+linear
DetrMat  = DataMat - (BetaFits * D01');    

% STEP 1: prepare for multiscale analysis
%
% get the cumulative histogram at each voxel
YPER = cumsum( DetrMat, 2 );
% list: #pieces timecourse can be subdivided into;
% need at least 3 timepoints per segment to estimate linear trend
klist   = 1:floor(size(YPER,2)./3);
% get the list of associated segment lengths
seglist = floor(size(YPER,2)./klist);

% identify cases where segments are same length, due to rounding
kep     = [1 double((seglist(2:end) - seglist(1:end-1)) ~= 0)]; % flag when not the same size as prev.
klist   = klist(kep>0);    % drop duplicates from list
seglist = seglist(kep>0);  % drop duplicates from list

% unit specifications?
if    ( isempty(time_units) )   scale_list = seglist;
elseif( isnumeric(time_units) ) scale_list = seglist .* time_units;
else   error('time_units field must be numeric! (or empty)');
end
% adjusting range?
if( ~isempty(scal_range) )

    if( scal_range(1) < min(scale_list) )
        disp(['User specified scaling range below minimum limit of T=',num2str(min(scale_list)),'. Readjusting...']);
        scal_range(1) = min(scale_list);
    end
    if( scal_range(2) > max(scale_list) )
        disp(['User specified scaling range above maximum limit of T=',num2str(max(scale_list)),'. Readjusting...']);
        scal_range(2) = max(scale_list);
    end
    
    disp(['...running DFA on reduced scaling range: [',num2str(scal_range(1)),', ',num2str(scal_range(2)),']']);
    % binary vector specifying timescales to retain
    keepin = (scale_list>=scal_range(1)) & (scale_list<=scal_range(2));
    
    if( sum(keepin>0) < 3 )
        error('scaling range is too constrained, you must have >2 points to fit slope. Try expanding it!');
    end
    
    % now discard points in list
    klist = klist( keepin ); %discard outside of range
    seglist = seglist( keepin ); %discard outside of range
    scale_list = scale_list( keepin ); %discard outside of range
else
    disp(['...running DFA on full scaling range: [',num2str(min(scale_list)),', ',num2str(max(scale_list)),']']);
end

disp(['...computed over K=',num2str(length(klist)),' timescale points']);

% initialize variance data matrix
FVAR    = zeros( Nvox, length(klist));

for(k=1:length(klist)) % k indep. splits

    Nelem = seglist(k); % number of time-points per segment
    Nsplt = klist(k);   % number of segments
    %
    % trim off non-multiple timepoints
    YPER_trim = YPER(:, 1:(Nsplt*Nelem));
    % reshape the matrix: vox x elems in segment x #splits
    YPER_trim = reshape( YPER_trim, Nvox,Nelem,Nsplt );
    % initialize the detrended matrix
    YPER_detr = zeros( Nvox,Nelem,Nsplt );

    % for each segment (split), detrend the matrix
    for(q=1:Nsplt)
        % linear-constant regressors
        D01 = [ ones(Nelem,1), linspace(-1,1,Nelem)' ];
        % regression coefficients
        BetaFits         = YPER_trim(:,:,q) * D01 * inv( D01' * D01 );
        % keep detrended data
        YPER_detr(:,:,q) = YPER_trim(:,:,q) - (BetaFits * D01');    
    end
    
    % compute RSS variation in the detrended splits, then take the average
    FVAR(:,k) = sqrt(sum(sum( YPER_detr.^2, 2), 3))./( Nsplt*Nelem );
end
 
disp('...done!');

% log-transform the RSS variation
log_FVAR = log(FVAR);
% linear fit with log( segment size )
REG = [ ones(length(klist),1), log(seglist(:)) ];
% Get the Beta -- linearity fit with log( segment size )
BetaFits  = log_FVAR * REG * pinv( REG' * REG ); 
% final results take only the linear fit (component 2)
H_vect = BetaFits(:,2);
