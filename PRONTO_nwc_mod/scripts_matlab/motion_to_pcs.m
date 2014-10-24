function PCset = motion_to_pcs( infile, outfile, fract, record )
%
% MOTION_TO_PCS:  this function takes in motion parameter estimate (MPE)
% text-file, performs PCA, and outputs the first "numpc" component
% timeseries that account for >"fract" fraction of motion variance
%
% Syntax:
%          motion_to_pcs( infile, outfile, fract, record )
% Input:
%          infile : string giving path/name of MPE textfile
%         outfile : string giving path/name of textfile where PCs are recorded
%           fract : determines number of PCs K to keep; PCs 1-K must
%                   account for >= fract of data variance
%          record : indicates 0=no output, 1=record PCs in 'outfile'
%
% Output:
%           PCset : set of PC timecourses, explaining >= fract data variance

% load text file
X = load(infile);
% mean-center timeseries
Y = bsxfun(@minus,X,mean(X));
% singular value decomposition on Y
[u s v] = svd(Y,'econ');
%
fractLine = diag(s.^2)./trace(s.^2);
numpc = sum( cumsum(fractLine) <= fract ) +1;
%
% Get first "npc" principal components:
PCset = u(:,1:numpc) * s(1:numpc,1:numpc);
% print regressors into new text file

% option: write to textfile
if( record>0)
    dlmwrite( outfile, PCset, 'delimiter', '\t', 'precision', '%.4f' );
end
%
disp( sprintf( '\nPCA with %d components, accounts for %.2f%% variance.\n', numpc, 100*sum(fractLine(1:numpc)) ) );


