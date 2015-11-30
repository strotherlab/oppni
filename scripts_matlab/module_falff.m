function output = module_falff( datamat, split_info )
%
% =========================================================================
% MODULE_FALFF: fractional amplitude of low-frequency fluctuations
% =========================================================================
%
%   Syntax:
%           output = module_falff( datamat, split_info )
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
Nhalf        = floor(Ntime/2);
TR           = round(split_info.TR_MSEC./1000);
% parameters for spectral estimation
Fny    = 0.5 * (1/TR);                % nyquist frequency
NFFT   = 2^nextpow2( Nhalf );         % Next power of 2 from length of time-axis
f      = Fny*linspace(0,1,NFFT/2+1);  % fourier data corresponds to these frequency points

kq=0;
for(t=1:10:Nhalf)

    kq=kq+1; %%increment
    % get sum of spectral power values, for f > dataInfo.FreqCut
    powMat  = abs( fft( datamat(:,t:t+Nhalf) , NFFT ,2) ) / Nhalf; 
    powMat  = powMat(:,1:NFFT/2+1);
    % fractional power
    falff(:,kq) = sum( powMat(:, (f>0.015)&(f<0.08)),2) ./ sum(powMat(:, (f>0.0001)),2);
end
% catch for bad-valued voxels
falff(~isfinite(falff))=eps;

% reproducibility estimates on "most independent" splits
[output.metrics.R rSPMZ] = get_rSPM( falff(:,1), falff(:,end), 1 );

if    ( strcmp( split_info.spm, 'zscore' ) ) output.images  = rSPMZ;
elseif( strcmp( split_info.spm, 'raw'   ) )  output.images  =  mean(falff,2);
else  error( 'invalid output type for falff (should be zscore or raw).' );
end

output.temp    = datamat'  * (rSPMZ ./ sqrt(sum(rSPMZ.^2)));
