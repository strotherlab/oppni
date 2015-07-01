function [windowAvg] = simple_averaging_for_ER( dataMat, onsets, params )
% 
% SIMPLE_AVERAGING_FOR_ER : procedure to (1) estimate stimulus-locked BOLD
% responses, then (2) average across stimuli to get stable split estimates
% of BOLD response.
%
% Syntax:
%           windowAvg = averaging_for_ER( dataMat, onsets, params )
% 
% Input:
%             dataMat : a 2D (voxels x time) data matrix
%              onsets : a 1D vector of stimulus onset times (IN MILLISECONDS!!)
%              params : a structure with the following sub-fields:
%                       
%                           params.TR    : acqusition time (IN MILLISECONDS)
%                           params.WIND  : time window size (integer, #TR)
%                           params.Nblock: number of data splits produced; must be an even number >=4
%                           params.delay : interpolation delay from slice-timing correction (usually TR/2)
%                           params.norm  : integer defining normalization on time-windows before averaging
%                                             { 0=do nothing / 1=subtract mean / 2=normalize variance }
%                                          
%                           ** if there is no good reason, just set norm=0
% Output:
%
%           windowAvg : a (voxels x WIND x Nblock) matrix, of "Nblock" time-averaged, stimulus-locked BOLD response blocks
%
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

% parameter settings in data
TR       = params.TR;     % in ms
WIND     = params.WIND;   % # windows
Nblock   = params.Nblock; % # splits
delay    = params.delay;
norm     = params.norm;
%
[Nvox Ntime] = size( dataMat ); % matrix dims
onsetsTR     = round( (onsets./TR) - (delay/TR) + 1 );
%
onsetsTR( onsetsTR > (Ntime-WIND+1) ) = [];
onsetsTR( onsetsTR <= 0 )         = [];
Nonsets = length(onsetsTR);
%
fullAvgSet = zeros( Nvox, WIND, Nonsets);
%
for(n=1:Nonsets)
    fullAvgSet(:,:,n) = dataMat(:,onsetsTR(n)-1+(1:WIND));
end

%% 

% normalization options
if( norm > 0 )

    fullAvgSet =  bsxfun(@minus, fullAvgSet , fullAvgSet(:,1,:) );
    % normalize: take SD at each voxel / get mean over all volumes
    if( norm==2 )
    fullAvgSet =  bsxfun(@rdivide,fullAvgSet, std( fullAvgSet,0,2 ) );
    fullAvgSet(~isfinite(fullAvgSet)) = 0;
    end
end

if( Nblock == 1 )
    % 
    windowAvg = mean( fullAvgSet, 3 );
else
    % 
    midOn     = onsetsTR + WIND/2;
    onList    = 1:Nonsets;
    Ntrunc    = floor( Ntime / Nblock );
    windowAvg = zeros( Nvox, WIND, Nblock );
    % for each split, now get windowed averaging
    for( kk=1:Nblock )

        CutLo = (kk-1)*Ntrunc + 1;
        CutHi = ( kk )*Ntrunc;
        % assign scans to a split, based on whether >50% of time-window is within the split:
        select = onList(  (midOn > CutLo) & (midOn < CutHi)  );
        % catch in case no stim in this block
        if( isempty( select ) ) 
            % take cut-off stimuli
            select = onList( midOn<CutLo ); 
            % if nothing found, terminate. Otherwise take last one in list
            if( isempty( select ) )
                error('No stimulus onsets found in one of the splits. Terminating.');
            else
                select = select(end); 
            end
        end
        % modified averaging for compatibility with octave
        windowAvg(:,:,kk) = sum( fullAvgSet(:,:,select), 3) ./ length(select);
    end
end
