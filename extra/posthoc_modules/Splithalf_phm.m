function out = Splithalf( datamat, design )
%
% . Splithalf Z-scored analysis, compares 1-group rel. 0, or between 2 groups
%

NBOOT = 100;

if    ( isempty(design) || numel(unique(design))==1 )
    
    disp('Split-half analysis, 1-sample...');
    
    % parameters
    n    = size(datamat,2);
    disp('running resampling...');
    sp_av = zeros( size(datamat,1), NBOOT );
    sp_df = zeros( size(datamat,1), NBOOT );
    for(bsr=1:NBOOT)
        [bsr NBOOT],
        list = randperm(n);
        sp1  = mean( datamat(:, list(1:round(n/2)))     ,2);
        sp2  = mean( datamat(:, list(round(n/2)+1:end)) ,2);
        %[r spm] = get_rSPM( sp1,sp2,1 );
        %bsrmat = bsrmat + spm./NBOOT;
        sp_av(:,bsr) = (sp1+sp2);
        sp_df(:,bsr) = (sp1-sp2);
    end
    out.bsr     = mean(sp_av,2)./std(sp_df,0,2);
    out.bsr_p   = 2*normcdf(out.bsr);
    
    out.testname = 'splithalf_1samp';

elseif( numel(unique(design))>1 && size(design,1)>1 && size(design,2)==1 )
    
    ix = unique(design);
    datamat1 = datamat(:,design==ix(1));
    datamat2 = datamat(:,design==ix(2));
    disp('Split-half analysis, 2-sample (unpaired)...');

    % parameters
    n1 = size(datamat1,2);
    n2 = size(datamat2,2);
    disp('running resampling...');
    bsrmat = zeros( size(datamat1,1), 1 );
    for(bsr=1:NBOOT)
        [bsr NBOOT],
        list1 = randperm(n1);
        list2 = randperm(n2);
        sp1  = mean( datamat2(:, list2(1:round(n2/2)))     ,2) - mean( datamat1(:, list1(1:round(n1/2)))     ,2);
        sp2  = mean( datamat2(:, list2(round(n2/2)+1:end)) ,2) - mean( datamat1(:, list1(round(n1/2)+1:end)) ,2);
        [r spm] = get_rSPM( sp1,sp2,1 );
        bsrmat = bsrmat + spm./NBOOT;
    end
    out.bsr     = bsrmat;
    out.bsr_p   = 2*normcdf(bsrmat);

    out.testname = 'splithalf_2samp_unpair';

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
    disp('Split-half analysis, 2-sample (paired)...');
    
    % parameters
    n    = size(datamat,2);
    disp('running resampling...');
    bsrmat = zeros( size(datamat,1), NBOOT );
    for(bsr=1:NBOOT)
        [bsr NBOOT],
        list = randperm(n);
        sp1  = mean( datamat2(:, list(1:round(n/2)))     ,2) - mean( datamat1(:, list(1:round(n/2)))     ,2);
        sp2  = mean( datamat2(:, list(round(n/2)+1:end)) ,2) - mean( datamat1(:, list(round(n/2)+1:end)) ,2);
        [r spm] = get_rSPM( sp1,sp2,1 );
        bsrmat = bsrmat + spm./NBOOT;
    end
    out.bsr     = bsrmat;
    out.bsr_p   = 2*normcdf(bsrmat);
    
    out.testname = 'splithalf_2samp_paired';

end

