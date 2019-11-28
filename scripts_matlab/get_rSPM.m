function [ rep, rSPM ] = get_rSPM( vect1, vect2, keepMean )
%
%    [ rep, rSPM ] = get_rSPM( vect1, vect2, keepMean )
%
%  * This script takes in 2 vectors, returns reproducibility (rep)
%    and reproducible SPM (rSPM)
%
%  * keepMean: reinsert the mean offset present in vectors
%
%
% ------------------------------------------------------------------------%
% Author: Nathan Churchill, University of Toronto
%  email: nathan.churchill@rotman.baycrest.on.ca
% ------------------------------------------------------------------------%
% version history: March 15 2012
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

inOctave = in_octave();
%if inOctave
% this section was removed as octave uses the corr function
% now, and its corrcoef function does not work completely.
%    for k = 1:size(vect1,2)
%        rep(k) = corrcoef(vect1(:,k), vect2(:,k));
%    end
%else
    for k =1:size(vect1,2)
        rep(k) = corr(vect1(:,k), vect2(:,k));
    end
%end

%(1) getting the mean offsets (normed by SD)
normedMean1 = mean(vect1)./std(vect1);
normedMean2 = mean(vect2)./std(vect2);
%    and rotating means into signal/noise axes
sigMean = (normedMean1 + normedMean2)/sqrt(2);
%noiMean = (normedMean1 - normedMean2)/sqrt(2);
%(2) getting  signal/noise axis projections of (zscored) betamaps
zvect1 = zscore(vect1);
zvect2 = zscore(vect2);
sigProj = ( zvect1 + zvect2)  / sqrt(2);
noiProj = ( zvect1 - zvect2)  / sqrt(2);
% noise-axis SD
noiStd = std(noiProj);
%(3) norming by noise SD:
%     ...getting the (re-normed) mean offsets
sigMean = sigMean./noiStd;
%noiMean = noiMean./noiStd; 
%  getting the normed signal/noise projection maps
sigProj = bsxfun(@rdivide,sigProj , noiStd);
%noiProj = noiProj ./ noiStd;

% Produce the rSPM:
if    ( keepMean == 1 )   rSPM = bsxfun(@plus, sigProj, sigMean);
elseif( keepMean == 0 )   rSPM = sigProj;
end



function inOctave = in_octave()

try
    ver_num = OCTAVE_VERSION;
    inOctave = 1;
    version_str = ['OCTAVE ' ver_num];
catch
    inOctave = 0;
    version_str  = ['MATLAB ' version];
end
