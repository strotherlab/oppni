function output = WM_weight( dataMat, dataInfo )
%
%     output = WM_weight( dataMat, dataInfo )
%

% Get basic data info
% ====================
% get data dimensions

% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';


Nruns = length( dataMat );
for(r=1:Nruns) 
    % dimensions of each split
    [Nvox Ntime(r,1)] = size( dataMat{r} ); 
    % subtract temporal means
    dataMat{r} = bsxfun(@minus,dataMat{r},mean(dataMat{r},2));
end

% initialize matrix values
WTset     = zeros(Nvox,Nruns);
percentGM = zeros(Nruns,1);
percentWM = zeros(Nruns,1);

for( r=1:Nruns ) %% estimate on each split

    % Inverse variance
    IvrMat  = 1./var(dataMat{r},0,2); 
    % index for non-ill-posed voxels
    ikeep = find( isfinite(IvrMat) );
    % if any ill-posed voxels, drop from the index
    if(length(ikeep)<Nvox)
        disp(['you have ',num2str(Nvox-length(ikeep)),' "empty" voxel timeseries; correcting, but check your masks!']);
        IvrMat = IvrMat(ikeep,:);
    end
    % parameters for threshold estimation
    NumList = (1:length(IvrMat))'; 
    Nqart   = round(Nvox/4);
    
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
        else
            WM_msk = double( dataInfo.priorMask > 0.5 );
            % drop out "ill-posed" voxels from mask
            if(length(ikeep)<Nvox)
               WM_msk = WM_msk(ikeep,1);
            end
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
            Overlap(k,1) = ( sum( Z.*WM_msk ) )./ sum( double((Z+WM_msk)>0) ); 
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

    % convert RSD values into [0,1] scale
    % where 0= greater than maxRSD (white matter threshold)
    %       1= less than minRSD (grey matter threshold)
    HI_dev =   RSD  - minRSD; 
    HI_dev(HI_dev<0) = 0; 
    HI_wt = HI_dev ./ maxRSD; 
    HI_wt(HI_wt>1)=1; 
    HI_wt = 1-HI_wt;

    % if missing voxels, adjust:
    if(length(ikeep)<Nvox)
        tmp = 0.5*ones(Nvox,1); % default setting of 0.5=equal likelihood of either type, all vox
        tmp(ikeep) = HI_wt; % for retained voxels, replace with actual values
        HI_wt = tmp; % swap in for original HI_wt
    end
    
    % recording split information:
    WTset(:,r)   = HI_wt;
    percentGM(r) = MinThr;
    percentWM(r) = 100 - MaxThr;
end

% ====== Output-1: physiological mask parameters ====== %

    output.WM_weight = mean( WTset, 2 ); % average non-neuronal weighting map
    output.WM_mask   = prod( double(WTset>=0.5), 2 ); % intersection of mask estimates
    % reproducibility across splits
    CC = triu( corr(WTset), 1 );
    output.WM_rep  = mean( CC(CC~=0) );
