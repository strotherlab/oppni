function diagnostic_fmri_pca( volname, maskname, mpename, outprefix )
%
% =========================================================================
% DIAGNOSTIC_FMRI_PCA: This script reads in fMRI 4D dataset, brain mask and 
% Motion Parameter Estimates (MPEs), and tests for potential outliers & 
% artifact. This provides a quick, computationally efficient test for major 
% motion spikes and QC issues.
% =========================================================================
%
% Syntax:
%           diagnostic_fmri_pca( volname, maskname, mpename, outprefix )
% Input:
%           volname  = string, giving path/name of fMRI data (NIFTI or ANALYZE format)
%           maskname = string, giving path/name of binary brain mask (excluding non-brain tissue)
%           mpename  = string, giving path/name of MPE file
%                     (e.g. volname = 'mypath/subject_data_1.nii')
%
%          outprefix = optional string specifying output path for 
%                      (1) diagnostic figures, (2) QC output. Default uses 'volname',
%                      just leave this entry empty (eg outprefix = [])
% Output:
%
%  - produces "output" structure, saved to "<outprefix>_QC_output.mat" with subfields:
%
%    [Censor vectors] binary vector (time x 1), where [1=non-outliers] and [0=outliers]. 
%        These are used as input to remove outliers in "interpolate_fmri":
%
%        output.censor_mot   : significant outlier in MPEs
%        output.censor_vol   : significant outlier in fMRI data
%        output.censor_volmot: significant outlier in BOTH fMRI data and MPEs (*)
%
%    [Censor matrices] Alternative outlier testing. Produces binary matrix (time x brain slices)
%        where (1=non-outliers) and (0=outliers). This gives outliers for individual
%        axial brain slices, can be used as input to remove outliers in "interpolate_fmri":
%
%        output.censor_slc   : significant outlier in fMRI slice data
%        output.censor_slcmot: significant outlier in BOTH fMRI slice data and MPEs (*)
%
%    (*) for (censor_volmot) and (censor_slcmot), we discard fMRI outliers
%        if they occur at same time OR 1-TR later than a motion spike. This
%        is to allow for slower fMRI-related signal changes (eg spin-history effects)
%
%    [Other possibly interesting outputs] 
%
%       output.eigfract_fmri : vector (K x 1) of fraction of variance explained by each PC of fMRI data 
%                              (K=PCs explaining 95% of variance total)
%       output.eigimages_fmri: matrix (voxel x K) PCA eigenimages for fMRI data 
%       output.eigvect_fmri  : matrix (time  x 10) PCA eigen-timeseries for fMRI data 
%
%       output.eigfract_mot  : vector (6 x 1) of fraction of variance explained by each PC of MPEs
%       output.eigweights_mot: matrix (6    x 6) PCA weights on motion parameters for MPEs
%       output.eigvect_mot   : matrix (time x 6) PCA eigen-timeseries for MPEs
%
%       output.global_signal : vector (time x 1) of global signal timecourse
%
%  - Also saves QC output figures. Includes:
%       "<outprefix>_diagnostic_plot0.png" : summary results for temporal variance in data
%       "<outprefix>_diagnostic_plot1.png" : summary results from PCA decomposition of data
%       "<outprefix>_diagnostic_plot2.png" : results for estimated motion spikes
%
% ------------------------------------------------------------------------------------------- 
% Description of Testing:
%
%  - Takes 4D fMRI data, converted into (voxels x time) matrix, and performs PCA
%    decomposition. Similarly, MPE matrix (time x 6 parameters) is decomposed 
%    with PCA. This gives an efficient, orthonormal representation of the data.
%  - Plots PCA information about data (eg eigenspectra) to characterize better
%
%  - For outlier detection: 
%    We want to find spikes (brief, large signal changes) in multidimensional data.
%    We estimate outliers in the MPE and fMRI data:
%
%    1. For each datapoint X(t) (1 <= t <= N_time):
%        (a) find median PC-space coordinate Xmed(t), in a 15 TR timewindow centered at t
%        (b) compute displacement by Euclidean distance D(t) = L2( X(t) - Xmed(t) )
%    2. For displacement timeseries D, compute Chi-square distribution on D^2
%        using Maximum Likelihood fit to infer degrees of freedom.
%    3. Identify outliers, significant at p<0.05
%
%  - We measure displacement relative to a 15-TR time window:
%    This has the advantage over (1) a simple derivative, e.g. we can discriminate 
%    displacement AWAY from the main datapoint cluster vs. discplement back TOWARDS
%    this cluster. The windowing also (2) minimizes the impact of slow changes
%    in amplitude over time, which might inflate our displacement estimates.
%
%  - If we do slice-specific outlier testing, this whole process is repeated independently 
%    for each axial brain slice
%  - A conservative apporach is to remove points that are outliers in both MPE and fMRI data
%
% ------------------------------------------------------------------------%
% Author: Nathan Churchill, University of Toronto
%  email: nathan.churchill@rotman.baycrest.on.ca
% ------------------------------------------------------------------------%
% version history: 2013/09/18
% ------------------------------------------------------------------------%
% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%

if( exist('OCTAVE_VERSION','builtin') )
    % load stats packages
    pkg load statistics;
    pkg load optim;
end
 
%% (0) PREPARATION

% formatting the output data name: make sure an output name exists
if( isempty(outprefix) ) 
    % if not, make it a trimmed version of output name
    idx     = strfind(volname,'/');
    % drop the suffix / input directory bits
    if( isempty( idx ) )
            outprefix = volname( 1           :(end-4) ); 
    else    outprefix = volname( (idx(end)+1):(end-4) ); 
    end
end

% (1) load the motion parameters (matrix of [time x 6] elements)
 MPE    = load(mpename);
% mean-center (subtract the mean from each timecourse)
 MPE    = bsxfun(@minus,MPE,mean(MPE));

% (2) load fMRI NIFTI-format data into MatLab
VV     = load_untouch_nii(volname); % load 4D fMRI dataset
MM     = load_untouch_nii(maskname);% load the binary brain mask
% convert mask into 3D matrix volume
mask   = double(MM.img);
% convert 4D fMRI data into (voxels x time) matrix
rawepi = convert_nii_to_mat( VV,MM ); 
% mean-center (subtract mean from each voxel
epimat = bsxfun(@minus,rawepi,mean(rawepi,2));
% get fMRI data matrix dimensions
[Nvox Ntime] = size(rawepi);

% catch #1: identify cases where #TRs does not match for fMRI and MPE
if( size( epimat,2 ) ~= size(MPE,1) )
    %
    error('ERROR: number of timepoints in fMRI data and MPE data do not match!');
end
% catch #2: identify cases where mask and 4D EPI volumes are not the same size
dimsV = size(VV.img(:,:,:,1));
dimsM = size(MM.img(:,:,:,1));
%
if( (dimsV(1)~=dimsM(1)) || (dimsV(2)~=dimsM(2)) || (dimsV(3)~=dimsM(3)) )
    %
    error('ERROR: dimensions of mask and 4D fMRI data do not match!');
end

%% (1) MOTION OUTLIER TESTING

% singular value decomposition (get the principal components of MPEs)
[u s v] = svd(MPE','econ');
   s2   = s.^2;               % squared eigenvalues
   p    = diag(s2)./trace(s2);% fraction of variance explained by each PC
% project data onto PC space (each timepoint = PC-space column vector)
Qmot    = (v * s)'; 

% run through timepoints (1...Ntime)
for(t=1:Ntime)
    % estimate 15-TR time window
    wind=(t-7):(t+7);
    if(wind( 1 )< 1   ) wind(wind<1)     = []; end
    if(wind(end)>Ntime) wind(wind>Ntime) = []; end
                        wind(wind==t)    = [];
	% get distance estimate at each timepoint, relative to median coordinate
    Dist_mot2(t,1)=sum(( Qmot(:,t) - median( Qmot(:,wind),2 ) ).^2);
end
%
% (1) Chi2 fit of MPEs (using gamma distribution, ML estimator)
Dist_mot2_norm = Dist_mot2 ./ max(Dist_mot2);
par_ab         = gamfit( Dist_mot2_norm );
% get Chi2 probability of each point + outlier threshold
probChi2_mot   = 1-gamcdf( Dist_mot2_norm, par_ab(1), par_ab(2) );
outThr_mot     = gaminv( 0.95, par_ab(1), par_ab(2) );

%% (2) FULL BRAIN-VOLUME OUTLIER TESTING

% singular value decomposition, on full fMRI data matrix
[V S2 temp] = svd( epimat'*epimat ); 
   S   = sqrt(S2);
   U   = epimat*V*inv(S);
   P   = diag(S2) ./ trace( S2 );
% project data onto PC space (each timepoint = PC-space column vector)
Qvol   = (V * S)';

% run through data timepoints (1...Ntime)
for(t=1:Ntime)
    % estimate 15-TR time window
    wind=(t-7):(t+7);
    if(wind( 1 )< 1   ) wind(wind<1)     = []; end
    if(wind(end)>Ntime) wind(wind>Ntime) = []; end
                        wind(wind==t)    = [];
	% get distance estimate at each timepoint, relative to median coordinate
    Dist_vol2(t,1)=sum(( Qvol(:,t) - median( Qvol(:,wind),2 ) ).^2);
    %
    Qtraj(:,t) = mean( Qvol(:,wind),2 );
end
% Chi2 fit full volume fMRI data (using gamma distribution, ML estimator)
Dist_vol2_norm = Dist_vol2 ./ max(Dist_vol2);
par_ab         = gamfit( Dist_vol2_norm );
% get Chi2 probability of each point + outlier threshold
probChi2_vol   = 1-gamcdf( Dist_vol2_norm, par_ab(1), par_ab(2) );
outThr_vol     = gaminv( 0.95, par_ab(1), par_ab(2) );

%% (3) SLICE-BASED OUTLIER TESTING

% get 3D volume dimensions
[Nx Ny Nz] = size( mask );
% get index-vector, where each voxel is labelled with its slice# (z=1...Nz)
ixvol = ones( [Nx Ny Nz] );
for(z=1:Nz) ixvol(:,:,z) = ixvol(:,:,z) .* z; end
ixvect = ixvol(mask>0);

% initialize distance/probability values
outThr_slc     =       ones(   1,  Nz  );    %threshold =1
probChi2_slc   =       ones( Ntime, Nz );%prob of outlier=1
Dist_slc2_norm = 0.001*ones( Ntime, Nz );%zero displacement

% iterate through each brain slice
for( z=1:Nz )
    
    % check for signal in slice
    fullsum = sum(sum(abs(epimat(ixvect==z,:))));
    
    % requires that there be signal in this slice
    if( (fullsum~=0) && isfinite( fullsum ) )
        
        % singular value decomposition on individual fMRI data slice
        [V_tmp S2_tmp temp] = svd( epimat(ixvect==z,:)'*epimat(ixvect==z,:) ); 
         S_tmp  = sqrt(S2_tmp./sum(ixvect==z));
        % project data onto PC space (each timepoint = PC-space column vector)
         Qslc   = (V_tmp * S_tmp)';
        % run through timepoints (1...Ntime)
        for(t=1:Ntime)
            % estimate 15-TR time window
            wind=(t-7):(t+7);
            if(wind( 1 )< 1   ) wind(wind<1)     = []; end
            if(wind(end)>Ntime) wind(wind>Ntime) = []; end
                                wind(wind==t)    = [];
            % get distance estimate at each timepoint, relative to median coordinate
            Dist_tmp2(t,1)=sum(( Qslc(:,t) - median( Qslc(:,wind),2 ) ).^2);
        end

        % (2) Chi2 fit, individual fMRI data slice (using gamma distribution, ML estimator)
        Dist_tmp2_norm = Dist_tmp2 ./ max(Dist_tmp2);
        par_ab         = gamfit( Dist_tmp2_norm );
        % get Chi2 probability of each point + outlier threshold
        probChi2_tmp   = 1-gamcdf( Dist_tmp2_norm, par_ab(1), par_ab(2) );
        outThr_slc(1,z)= gaminv( 0.95, par_ab(1), par_ab(2) );
        %
        Dist_slc2_norm(:,z) = Dist_tmp2_norm;
        probChi2_slc(:,z)   = probChi2_tmp;
    end
end

 % Rescale slice-based distance values, so that =1 indicates significant outlier
 Dist_slc2_rescal = bsxfun(@rdivide,Dist_slc2_norm,outThr_slc);
 % Get significant outliers, based on each metric
 outlier_mot   = double( probChi2_mot <= 0.05 );
 outlier_vol   = double( probChi2_vol <= 0.05 );
 outlier_slc   = double( probChi2_slc <= 0.05 );

 % find outliers in fMRI data, that have motion outlier concurrently, or 1TR before
 % -accounts for potential delay effects
 outlier_delay = double( (outlier_mot + [0; outlier_mot(1:end-1)]) > 0 );
 outlier_volmot= outlier_vol .* outlier_delay;
 outlier_slcmot= bsxfun(@times,outlier_slc,outlier_delay);
 
% outliers based on volume/motion only
  maxout    = round( 0.10*Ntime );
% -------- test for excessive outliers ------- %
if( sum(outlier_mot ) > maxout ) disp( 'WARNING: more than 10% of motion points are outliers!'); end
if( sum(outlier_vol ) > maxout ) disp( 'WARNING: more than 10% of fMRI data points are outliers!'); end

%% (4) PLOTTING RESULTS
gsvect  = zscore( mean(epimat)' ); 

% Only plot output if operating in matlab
if( exist('OCTAVE_VERSION','builtin') )
    %
    disp(['figures not plotted in Octave']);
    % just compute global signal regressor
else %if (usejava('desktop') && usejava('jvm'))
    % ----------- (Fig.0) Display Temporal Variance information ----------- %
    h0=figure('visible','off');
    set(gcf, 'Units', 'normalized');
    set(gcf, 'Position', [0.10 0.15 0.60 0.65]);

    stdvol = std(double(VV.img),0,4);
    slcx = permute( stdvol(round(Nx/2),:,:), [3 2 1] );
    slcy = permute( stdvol(:,round(Ny/2),:), [3 1 2] );
    slcz = stdvol(:,:,round(Nz/2));

    slcall= [slcx(:); slcy(:); slcz(:)];

    subplot(2,3,1); imagesc( flipud(slcx), prctile( slcall(slcall~=0), [2.5 97.5] ));
    subplot(2,3,2); imagesc( flipud(slcy), prctile( slcall(slcall~=0), [2.5 97.5] ));
    title('Temporal Standard Deviation plots');
    subplot(2,3,3); imagesc( slcz, prctile( slcall(slcall~=0), [2.5 97.5] ));

    subplot(2,1,2);
    % pre-get vectors
    fpcvect = zscore( V(:,1) ); fpcvect = fpcvect .* sign( corr(fpcvect, gsvect) );
    mpcvect = zscore( v(:,1) ); mpcvect = mpcvect .* sign( corr(mpcvect,fpcvect) );
    %Current leftout: plots the relationship between fMRI and MPE first PCs
    plot( 1:Ntime,mpcvect,'.-b' , 1:Ntime,fpcvect,'.-k', 1:Ntime, gsvect, '.-r', 'markersize',8, 'linewidth',1.5 ); 
    cc1 = abs(corr( mpcvect, fpcvect )); cc1 = round(cc1*100)/100;
    cc2 = abs(corr(  gsvect, fpcvect )); cc2 = round(cc2*100)/100;
    legend('MPE-PC','fMRI-PC','fMRI-GS');
    text(10, 1.1*max( [mpcvect;fpcvect; gsvect] ), strcat('correlation (fmriPC,motionPC): ',num2str(cc1)));
    text(10, 1.1*min( [mpcvect;fpcvect; gsvect] ), strcat('correlation (fmriPC,fmriGS): ',num2str(cc2)));
    title('Principal Component #1, and fMRI global-mean timecourses');
    xlabel('time (TR)');
    ylabel('z-score units');

    saveas(h0,strcat(outprefix,'_diagnostic_plot0.png'));
    close(h0);

    % ----------- (Fig.1) Display PCA-based information ----------- %

    h0=figure('visible','off');
    set(gcf, 'Units', 'normalized');
    set(gcf, 'Position', [0.10 0.05 0.60 0.85]);

    subplot(2,2,1); % plot1: plot datapoints in pca space (dim1+2)
    hold on;
    pbound = 1.1*max(abs(Qvol(1,:)));
    plot( pbound*[-1 1], [0 0], '-k', [0 0], pbound*[-1 1], '-k', 'color',[0.75 0.75 0.75] );
    plot( Qvol(1,:), Qvol(2,:), 'ok', 'markerfacecolor','b', 'markersize',4);
    plot( Qtraj(1,:), Qtraj(2,:), '--r', 'linewidth',2);
    text( Qtraj(1,1), Qtraj(2,1), ['T=1'],'color','r' );
    text( Qtraj(1,end), Qtraj(2,end), ['T=',num2str(Ntime)],'color','r' );
    title('Scans in PCA-space'); xlabel('PC-1'); ylabel('PC-2');
    %
    xlim([-pbound pbound]);
    ylim([-pbound pbound]);

    subplot(2,2,2); % plot3: eigenspectrum
    plot( 1:6,p,'o-b', 1:10, P(1:10), 'o-k', 'markersize',6, 'linewidth',1.5 ); 
    ylim([0 1]); xlim([0.5 10.5]);
    legend('MPEs','fMRI');
    title('Fractional eigenspectrum'); xlabel('PC#'); ylabel('fraction of total variance');
    set(gca,'Ytick',[0.1:0.1:1.0],'gridlinestyle',':','ygrid','on')

    % ---------------------------

    vertvect = sum(sum(mask,1),2);
    aidx     = find( vertvect(:) > 0 );
    seg      = floor(length(aidx)./5);
    %
    tmp=mask;tmp(tmp>0)=U(:,1); pbound = prctile( abs(U(:,1)), 95 );
    slc = [tmp(:,:,aidx(1)+1*seg) tmp(:,:,aidx(1)+2*seg); tmp(:,:,aidx(1)+3*seg) tmp(:,:,aidx(1)+4*seg)];
    subplot(2,2,3); imagesc( slc, [-pbound pbound] );
    title('PC Eigenimage #1');
    %
    tmp=mask;tmp(tmp>0)=U(:,2); pbound = prctile( abs(U(:,2)), 95 );
    slc = [tmp(:,:,aidx(1)+1*seg) tmp(:,:,aidx(1)+2*seg); tmp(:,:,aidx(1)+3*seg) tmp(:,:,aidx(1)+4*seg)];
    subplot(2,2,4); imagesc( slc, [-pbound pbound] );
    title('PC Eigenimage #2');

    saveas(h0,strcat(outprefix,'_diagnostic_plot1.png'));
    close(h0);

    % ----------- (Fig.2) Display Motion-linked information ----------- %

    h0=figure('visible','off');
    set(gcf, 'Units', 'normalized');
    set(gcf, 'Position', [0.10 0.05 0.60 0.85]);
    % 
    subplot(3,2,1);
    plot( 1:Ntime, sqrt(Dist_mot2_norm), '.-b', [1 Ntime], [1 1]*sqrt(outThr_mot), ':r' );
    title('MPE displacement (fraction of max.)');
    xlabel('time (TR)'); ylabel('amplitude (a.u.)');
    xlim([0 Ntime+1]); ylim([-0.05 1.05]);
    %
    subplot(3,2,3);
    plot( 1:Ntime, sqrt(Dist_vol2_norm), '.-b', [1 Ntime], [1 1]*sqrt(outThr_vol), ':r' );
    title('FMRI volume displacement (fraction of max.)');
    xlabel('time (TR)'); ylabel('amplitude (a.u.)');
    xlim([0 Ntime+1]); ylim([-0.05 1.05]);
    %
    subplot(3,2,5);
    imagesc( [outlier_mot';outlier_vol';2.*outlier_volmot'], [0 1.5] );
    title('significant outliers');
    text(-10, 1,'MPE');
    text(-10, 2,'fMRI');
    text(-10, 3,'Both');
    xlabel('time (TR)');
    set(gca,'YTickLabel',['';'';'']);
    %
    subplot( 2,2,2 );
    imagesc( Dist_slc2_rescal', [0 1] );
    title('FMRI slice displacements (fraction of max.)');
    xlabel('time (TR)');
    ylabel('slice');
    %
    subplot( 2,2,4 );
    imagesc( outlier_slc' + outlier_slcmot', [0 1.5] );
    title('significant outliers');
    xlabel('time (TR)');
    ylabel('slice');

    saveas(h0,strcat(outprefix,'_diagnostic_plot2.png'));
    close(h0);
% else
%     disp(['No display found in Matlab']);
end


% Some output text:

disp(sprintf('\n\nDone. Diagnostic results...\n')),
numsignif = sum( cumsum(P) < 0.95 )+1;
disp(strcat('>  1-',num2str(numsignif),' PCs required to explain 95% of fMRI variance'));

if( sum( outlier_mot ) ==0 ) disp( '(1) No motion outliers identified.' );
else                         disp( '(1) Possible motion outliers at timepoints:' ); disp( find( outlier_mot(:)' > 0 ) );
end
%
if( sum( outlier_volmot ) ==0 ) disp( '(2) No fMRI volume +motion outliers identified.' );
else                         disp( '(2) Possible fMRI volume +motion outliers at timepoints:' );   disp( find( outlier_volmot(:)' > 0 ) );
end
if( sum( sum(outlier_slcmot,2) ) ==0 ) disp( '(3) No fMRI slice +motion outliers identified.' );
else                         disp( '(3) Possible fMRI slice +motion outliers at timepoints:' );   disp( find( sum(outlier_slcmot,2)' > 0 ) );
end

%% (5) Define "output" structure and save

% censor files, where 1=non-outlier, and 0=significant outlier
output.censor_mot    = 1-outlier_mot;
output.censor_vol    = 1-outlier_vol;
output.censor_slc    = 1-outlier_slc;
output.censor_volmot = 1-outlier_volmot;
output.censor_slcmot = 1-outlier_slcmot;
% eigenspectrum data
output.eigfract_fmri  = P(  1:numsignif );
output.eigfract_mot   = p;
% eigenvector data
output.eigimages_fmri= U(:, 1:numsignif);
output.eigvect_fmri  = V(:, 1:numsignif);
output.eigweights_mot= u;
output.eigvect_mot   = v;
% global signal regressor
output.global_signal = gsvect;

% save results to matfile
if( exist('OCTAVE_VERSION','builtin') )
     % matlab-compatible
     save( strcat(outprefix,'_QC_output.mat'), 'output', '-mat7-binary' );    
else save( strcat(outprefix,'_QC_output.mat'), 'output' );
end

%%
function dataMat = convert_nii_to_mat( niiVol, niiMask )
%
% take niiVol and niiMask (nifti) format, and convert
% to matlab vector/matrix structure:
%
% dataMat = nifti_to_mat( niiVol, niiMask )
%
vol = double(niiVol.img);
msk = double(niiMask.img);

dataMat = zeros( sum(msk(:)>0), size(vol,4) );

for(t=1:size(vol,4))
    tmp=vol(:,:,:,t);
    dataMat(:,t) = tmp(msk>0);
end
