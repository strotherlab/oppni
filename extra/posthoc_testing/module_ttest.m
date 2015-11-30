function out = module_ttest( datamat, type )

%
% 1-group, 2-group, 2-group(paired)
%

if    ( strcmpi(type,'1sample') )
    if(length(datamat)~=1) error('1-sample t-test can only be done on a single group'); end
    %
    n    = size(datamat{1},2);
    out.tmap = sqrt(n) .* mean(datamat{1},2) ./ std(datamat{1},0,2);
    out.pmap = 2.*tcdf( -abs(out.tmap), n-1 ); %2-tailed likelihood
    
elseif( strcmpi(type,'2sample') )
    if(length(datamat)~=2) error('2-sample t-test can only be done on 2 groups'); end
    % parameters
    n1 = size(datamat{1},2);
    n2 = size(datamat{2},2);
    mu = mean(datamat{2},2)-mean(datamat{1},2);
    sd = sqrt( ( (n1-1).*var(datamat{1},0,2) + (n2-1).*var(datamat{2},0,2) )./( n1+n2-2 ) );
    
    out.tmap = sqrt( n1*n2/(n1+n2) ) .* mu ./ sd;
    out.pmap = 2.*tcdf( -abs(out.tmap), n1+n2-2 ); %2-tailed likelihood

elseif( strcmpi(type,'paired') )
    if(length(datamat)~=2) error('paired (2-sample) t-test can only be done on 2 groups'); end
    if(size(datamat{1},2)~=size(datamat{2},2)) error('paired t-test groups do not match!'); end
    %
    n    = size(datamat{1},2);
    dif  = datamat{2} - datamat{1};
    out.tmap = sqrt(n) .* mean(dif,2) ./ std(dif,0,2);
    out.pmap = 2.*tcdf( -abs(out.tmap), n-1 ); %2-tailed likelihood
end
