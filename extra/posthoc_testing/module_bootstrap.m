function out = module_bootstrap( datamat )

if    ( strcmpi(type,'1sample') )
    if(length(datamat)~=1) error('1-sample bootstrap can only be done on a single group'); end
    %
    n    = size(datamat{1},2);
    disp('running bootstrap resampling...');
    bsrmat = zeros( size(datamat{1},1), 500 );
    for(bsr=1:500)
        list = ceil(n*rand(n,1));
        bsrmat(:,bsr) = mean(datamat{1}(:,list),2);
    end
    out.bsrmap = mean(bsrmat,2)./std(bsrmat,0,2);
    out.pmap   = min([sum(bsrmat>0,2) sum(bsrmat<0,2)],[],2)./bsr;
    
elseif( strcmpi(type,'2sample') )
    if(length(datamat)~=2) error('2-sample bootstrap can only be done on 2 groups'); end
    % parameters
    n1 = size(datamat{1},2);
    n2 = size(datamat{2},2);
    disp('running bootstrap resampling...');
    bsrmat = zeros( size(datamat{1},1), 500 );
    for(bsr=1:500)
        list1 = ceil(n1*rand(n1,1));
        list2 = ceil(n2*rand(n2,1));        
        bsrmat(:,bsr) = mean(datamat{2}(:,list2),2) - mean(datamat{1}(:,list1),2);
    end
    out.bsrmap = mean(bsrmat,2)./std(bsrmat,0,2);
    out.pmap   = min([sum(bsrmat>0,2) sum(bsrmat<0,2)],[],2)./bsr;    

elseif( strcmpi(type,'paired') )
    if(length(datamat)~=2) error('paired (2-sample) bootstrap can only be done on 2 groups'); end
    if(size(datamat{1},2)~=size(datamat{2},2)) error('paired bootstrap groups do not match!'); end
    %
    n    = size(datamat{1},2);
    dif  = datamat{2} - datamat{1};
    disp('running bootstrap resampling...');
    bsrmat = zeros( size(dif,1), 500 );
    for(bsr=1:500)
        list = ceil(n*rand(n,1));
        bsrmat(:,bsr) = mean(dif(:,list),2);
    end
    out.bsrmap = mean(bsrmat,2)./std(bsrmat,0,2);
    out.pmap   = min([sum(bsrmat>0,2) sum(bsrmat<0,2)],[],2)./bsr;
end
