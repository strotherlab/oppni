function [ bold_hrf ] = design_to_hrf( design, TR, PSET )
%
% =========================================================================
% DESIGN_TO_HRF:  This function convolves your (binary) task-design matrix with AFNI's 
% standard gamma-variate basis function 'SPMG1' - the 1-parameter estimated function.
% =========================================================================
%
% Syntax:
%               [ bold_hrf ] = design_to_hrf( design, TR, PSET )
%
% Input:
%       design: a (time x K) matrix, where K= number of different stimulus types
%       TR:     rate of volume acquisition (in sec.)
%       PSET:   2-element vector [P1 P2], where  P1 = rise time / P2 = falloff time (in sec.)
%       -->     recommended settings: [P1 = 5, P2 = 15]
%
% Output:
%       bold_hrf: a (time x K) matrix of canonical HRF responses to 'design' stimulus
%
%  SPMG1 design based on:
%  3dDeconvolve: AFNI version=AFNI_2010_10_19_1028 (Feb 22 2011) [32-bit]
%  Authored by:  B. Douglas Ward, et al.
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

[T K]     = size(design); % design matrix size (T=time, K=#stimuli)
bold_hrf  = zeros(T,K);   % initialize the convolved HRF matrix
t = (0:TR:20)';           % timepoints used to specify the HRF (in sec.)
                          % length ~ 20s should be sufficient fall-off time
% --- HRF parameters ---
A1 = 0.0083333;
A2 = 1.274527E-13;
P1 = PSET(1);             % positive lobe rise time
P2 = PSET(2);             % undershoot time
% Canonical HRF transfer function
spmg1 = exp(-t).*(A1.*t.^P1-A2*t.^P2);
% Now, convolve each design vector (columns of matrix) with 'SPMG1'
for( kk = 1:K )
    Resp            = conv( spmg1, design(:,kk) ); 
    bold_hrf(:,kk) = Resp(1:T);          % truncate endpiece (eg convol. tail)
end
