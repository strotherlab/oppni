function out = Ttest( datamat, design )
%
% . T-test analysis, compares 1-group rel. 0, or between 2 groups
%

if    ( isempty(design) || numel(unique(design))==1 )
    
    disp('T-test, 1-sample...');
    
    % parameters
    n    = size(datamat,2);
    out.tstat    = sqrt(n) .* mean(datamat,2) ./ std(datamat,0,2);
    out.tstat_p  = 2.*tcdf( -abs(out.tstat), n-1 ); %2-tailed likelihood
    
    out.testname = 'ttest_1samp';
    
elseif( numel(unique(design))>1 && size(design,1)>1 && size(design,2)==1 )
    
    ix = unique(design);
    datamat1 = datamat(:,design==ix(1));
    datamat2 = datamat(:,design==ix(2));
    disp('T-test, 2-sample (unpaired)...');

    % parameters
    n1 = size(datamat1,2);
    n2 = size(datamat2,2);
    mu = mean(datamat2,2)-mean(datamat1,2);
    sd = sqrt( ( (n1-1).*var(datamat1,0,2) + (n2-1).*var(datamat2,0,2) )./( n1+n2-2 ) );
    
    out.tstat    = sqrt( n1*n2/(n1+n2) ) .* mu ./ sd;
    out.tstat_p  = 2.*tcdf( -abs(out.tstat), n1+n2-2 ); %2-tailed likelihood

    out.testname = 'ttest_2samp_unpair';

elseif( numel(unique(design))>1 && size(design,1)>1 && size(design,2)==2 )
    
    ix = unique(design(:,1));
    datamat1 = datamat(:,design(:,1)==ix(1)); 
    datamat2 = datamat(:,design(:,1)==ix(2)); 
    
    clear datamat; %% wipe out old datamat
    
    des1     = design(design(:,1)==ix(1),2);
    des2     = design(design(:,1)==ix(2),2);
    
    if( size(datamat1,2) ~= size(datamat2,2) ) error('unequal splits, cannot match'); end
    
    for(i=1:length(des1))
        datamat(:,i) = datamat1(:,i) - datamat2(:,des2==des1(i)); 
    end
    disp('T-test, 2-sample (paired)...');
    % parameters
    n    = size(datamat,2);
    out.tstat    = sqrt(n) .* mean(datamat,2) ./ std(datamat,0,2);
    out.tstat_p  = 2.*tcdf( -abs(out.tstat), n-1 ); %2-tailed likelihood

    out.testname = 'ttest_2samp_paired';
end
