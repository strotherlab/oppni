function output = WM_weight( dataMat, dataInfo )
%
% Work In Progress.
%

% Get basic data info
% ====================
% get data dimensions
Nruns = length( dataMat );
for(r=1:Nruns) 
    % dimensions of each split
    [Nvox Ntime(r,1)] = size( dataMat{r} ); 
    % subtract temporal means
    dataMat{r} = bsxfun(@minus,dataMat{r},mean(dataMat{r},2));
end

% parameters for threshold estimation
NumList = (1:Nvox)'; 
Nqart   = round(Nvox/4);
% initialize matrix values
WTset     = zeros(Nvox,Nruns);
RSDset    = zeros(Nvox,Nruns);
percentGM = zeros(Nruns,1);
percentWM = zeros(Nruns,1);

for( r=1:Nruns ) %% estimate on each split

    % Inverse variance
    IvrMat  = 1./var(dataMat{r},0,2); 
    % Get Ordered voxel variance values
    IvrSort = sort(IvrMat);
    IvrIndx = sortrows([IvrMat NumList],1); 
    IvrIndx = IvrIndx(:,2);
    % Get linear-quad fit and estimated deviation
    P            = polyfit(NumList(Nqart:end-Nqart),IvrSort(Nqart:end-Nqart),1);
    IvrLQ        = polyval(P,NumList);
    RSD          = zeros( length( IvrMat ),1 );
    RSD_sorted   = IvrSort-IvrLQ;
    RSD(IvrIndx) = RSD_sorted;
    %
    stdlin = std( RSD_sorted(Nqart:end-Nqart) );
    NCut   = norminv(0.99,0,stdlin);
    %
    MinThr        = 100 * sum( RSD < NCut ) ./ length( RSD );
    % threshold RSD value
    minRSD =  prctile( RSD, MinThr );    
    
% =========================================================================
%   identifying the "non-neuronal tissue" threshold
%
%   --- method of thresholding depends on user choice ---

    if    ( strcmp( dataInfo.thresh_method, 'nothreshold' ) )
        %
        % no voxels declared non-neuronal (most conservative)
        MaxThr = 100;
        %
    elseif( strcmp( dataInfo.thresh_method, 'noprior'     ) )
        %
        % top 5% of voxels declared non-neuronal
        MaxThr = 95.0;
        %
    elseif( strcmp( dataInfo.thresh_method, 'prior' )       )
        %
        if( ~isfield( dataInfo.priorMask ) )
            disp( 'Need to specify binary spatial prior to make this choice' );
            return;
        end
        %
        % find threshold that maximizes overlap with the binary "prior" mask
        k=0;
        % list of percentile thresholds to test
        PCTLIST = MinThr(1):0.20:99.50;
        % initialize list of overlap values
        Overlap = zeros( length(PCTLIST), 1 );
        % test different thresholds
        for( qc=PCTLIST )
            %
            k=k+1;
            % map of supra-threshold voxels (white matter regions)
            Z=double(RSD>prctile(RSD,qc));        
            % Jaccard overlap of (thresholded RSD, binary prior mask)
            Overlap(k,1) = ( sum( Z.*NN_msk ) )./ sum( double((Z+NN_msk)>0) ); 
        end

        % convolve with simple running average smoother (width of 1%), to avoid local minima in curve
        Overlap = conv( Overlap, ones(1,5)./5 ); 
        Overlap = Overlap(3:end-2);
        % point of maximized overlap
        [omax imax]   = max( Overlap ); 
        % percentile of maximized overlap
        MaxThr = PCTLIST(imax);
        % catch for extreme low values
        if( MaxThr < 90.0 ) 
            disp( 'Warning: non-neuronal threshold is overly conservative at_',num2str(MaxThr),'%...readjusting to_90%');
            MaxThr = 90.0; 
        end
    else
        disp( 'Error: need to specify thresholding method.');
        return;
    end

    % threshold RSD value
    maxRSD =  prctile( RSD, MaxThr );
    
% =========================================================================
%   converting into map of spatial weights, of values [0, 1]

    % get non-neuronal tissue mask (general purpose),
    % for regions with any appreciable non-neuronal content
    % * not really necessary in current algorithm - maybe useful though *
    WMmask(:,r) = double(RSD > minRSD);
    
    % convert RSD values into [0,1] scale
    % where 0= greater than maxRSD (white matter threshold)
    %       1= less than minRSD (grey matter threshold)
    HI_dev =   RSD  - minRSD; 
    HI_dev(HI_dev<0) = 0; 
    HI_wt = HI_dev ./ maxRSD; 
    HI_wt(HI_wt>1)=1; 
    HI_wt = 1-HI_wt;

    % recording split information:
    WTset(:,r)   = HI_wt;
    RSDset(:,r)  = RSD;        
    percentGM(r) = MinThr;
    percentWM(r) = 100 - MaxThr;
end

% ====== Output: ====== %

output.WM_weight = mean( WTset, 2 );
output.WM_mask   = prod( WMmask, 2 );
% reproducibility across splits
CC = triu( corr(WTset), 1 );
output.WM_rep  = mean( CC(CC~=0) );
