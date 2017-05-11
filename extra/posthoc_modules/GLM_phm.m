function out = GLM( datamat, design )

    disp('GLM, t-statistic...');

    n    = size(datamat,2);
    k    = size(design, 2);
    op   = GLM_model_fmri( datamat, 0, [], design, [], [] );
    
    % parameters
    out.tstat    = op.Tmap_signl;
    out.tstat_p  = 2.*tcdf( -abs(out.tstat), n-k-1 ); %2-tailed likelihood
    
    out.testname = 'glm_tstat';
