function out = Bootstrap( datamat, design )
%
% . Bootstrapped analysis, compares 1-group rel. 0, or between 2 groups
%

NBOOT = 500;


if    ( isempty(design) || numel(unique(design))==1 )
    
    disp('Bootstrap, 1-sample...');
    
    % parameters
    n    = size(datamat,2);
    disp('running resampling...');
    bsrmat = zeros( size(datamat,1), NBOOT );
    for(bsr=1:NBOOT)
        [bsr NBOOT],
        list = ceil(n*rand(n,1));
        bsrmat(:,bsr) = mean(datamat(:,list),2);
    end
    out.bsr     = mean(bsrmat,2)./std(bsrmat,0,2);
    out.bsr_p   = 2*min([sum(bsrmat>0,2) sum(bsrmat<0,2)],[],2)./bsr;
    
    out.testname = 'bootstrap_1samp';

elseif( numel(unique(design))>1 && size(design,1)>1 && size(design,2)==1 )
    
    ix = unique(design);
    datamat1 = datamat(:,design==ix(1));
    datamat2 = datamat(:,design==ix(2));
    disp('Bootstrap, 2-sample (unpaired)...');

    % parameters
    n1 = size(datamat1,2);
    n2 = size(datamat2,2);
    disp('running resampling...');
    bsrmat = zeros( size(datamat1,1), NBOOT );
    for(bsr=1:NBOOT)
        [bsr NBOOT],
        list1 = ceil(n1*rand(n1,1));
        list2 = ceil(n2*rand(n2,1));        
        bsrmat(:,bsr) = mean(datamat2(:,list2),2) - mean(datamat1(:,list1),2);
    end
    out.bsr     = mean(bsrmat,2)./std(bsrmat,0,2);
    out.bsr_p   = 2*min([sum(bsrmat>0,2) sum(bsrmat<0,2)],[],2)./bsr;    

    out.testname = 'bootstrap_2samp_unpair';

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
    disp('Bootstrap, 2-sample (paired)...');
    
    % parameters
    n    = size(datamat,2);
    disp('running resampling...');
    bsrmat = zeros( size(datamat,1), NBOOT );
    for(bsr=1:NBOOT)
        [bsr NBOOT],
        list = ceil(n*rand(n,1));
        bsrmat(:,bsr) = mean(datamat(:,list),2);
    end
    out.bsr    = mean(bsrmat,2)./std(bsrmat,0,2);
    out.bsr_p  = 2*min([sum(bsrmat>0,2) sum(bsrmat<0,2)],[],2)./bsr;
    
    out.testname = 'bootstrap_2samp_paired';

end

