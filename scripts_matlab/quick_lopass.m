function [ X_filt ] = quick_lopass( X, TR )

% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 169 $';
% CODE_DATE    = '$Date: 2014-12-15 18:09:33 -0500 (Mon, 15 Dec 2014) $';
% ------------------------------------------------------------------------%

[Nvox Ntime] = size(X); % matrix dimensions
% using simple Butterworth filter -- linear phase/ flat frequency, rolloff not great but this is tolerable for fmri
Wp = (2*TR)*0.08; % passband is below 0.08 Hz
Ws = (2*TR)*0.10; % stopband is above 0.10 Hz
% filter design: max passband attn. =50% / min stopband attn =1%
[Nord, Wcut] = buttord( Wp, Ws, 3,10 );
% lowpass butterworth filter with desired cutoff
[B1,A1] = butter(Nord,Wcut,'low');
% zero-phase forward/reverse filter
X_filt  = filtfilt( B1,A1, X' )';


end