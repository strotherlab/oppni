function [denoised_data] = apply_glm(data,Regressors)

% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%

if ~iscell(data)  % use when we are running single subject analysis splitting in half
    if ~isfield(Regressors,'DET')
        Regressors.DET = 0;
    end
    if ~isfield(Regressors,'MP')
        Regressors.MP = [];
    end
    if ~isfield(Regressors,'GSPC1')
        Regressors.GSPC1 = [];
    end
    if ~isfield(Regressors,'PHYPLUS')
        Regressors.PHYPLUS = [];
    end
    if ~isfield(Regressors,'Signal')
        Regressors.Signal = [];
    end
    if ~isfield(Regressors,'NOISEROI')
        Regressors.NOISEROI = [];
    end
    
    det_sp1 = get_legendre(Regressors.DET,ceil(size(data,2)/2));
    det_sp2 = get_legendre(Regressors.DET,size(data,2) - ceil(size(data,2)/2));
    xdet_sp1 = [det_sp1;zeros(size(det_sp2,1),size(det_sp1,2))];
    xdet_sp2 = [zeros(size(det_sp1,1),size(det_sp2,2));det_sp2];
    
    Xn = [Regressors.MP xdet_sp1 xdet_sp2];
    Xall = [ Xn Regressors.Signal];
    BetaWeights = data*Xall*inv(Xall'*Xall);
    idxNoise    = 1:size(Xn,2);
    noi_estim   = BetaWeights(:,idxNoise) * Xall(:,idxNoise)';
    denoised_data   = data - noi_estim;
    
    if ~isempty(Regressors.NOISEROI)
        Xn   = [Regressors.NOISEROI];
        Xall = [Xn Regressors.Signal];
        BetaWeights = denoised_data*Xall*inv(Xall'*Xall);
        idxNoise    = 1:size(Xn,2);
        noi_estim   = BetaWeights(:,idxNoise) * Xall(:,idxNoise)';
        denoised_data   = denoised_data - noi_estim;
    end
    
    if ~isempty(Regressors.GSPC1)
        Xn   = [Regressors.GSPC1];
        Xall = [Xn Regressors.Signal];
        BetaWeights = denoised_data*Xall*inv(Xall'*Xall);
        idxNoise    = 1:size(Xn,2);
        noi_estim   = BetaWeights(:,idxNoise) * Xall(:,idxNoise)';
        denoised_data   = denoised_data - noi_estim;
    end
    
    
    if ~isempty(Regressors.PHYPLUS)
        Xall   = [Regressors.PHYPLUS];
        BetaWeights     = denoised_data*Xall*inv(Xall'*Xall);
        noi_estim       = BetaWeights * Xall';
        denoised_data   = denoised_data - noi_estim;
    end
else % use when we are running group analysis or multi-run, no spliting
    for run_counter = 1:numel(data)
        if ~isfield(Regressors{run_counter},'DET')
            Regressors{run_counter}.DET = 0;
        end
        if ~isfield(Regressors{run_counter},'MP')
            Regressors{run_counter}.MP = [];
        end
        if ~isfield(Regressors{run_counter},'GSPC1')
            Regressors{run_counter}.GSPC1 = [];
        end
        if ~isfield(Regressors{run_counter},'PHYPLUS')
            Regressors{run_counter}.PHYPLUS = [];
        end
        if ~isfield(Regressors{run_counter},'Signal')
            Regressors{run_counter}.Signal = [];
        end
        if ~isfield(Regressors{run_counter},'NOISEROI')
            Regressors{run_counter}.NOISEROI = [];
        end
        det = get_legendre(Regressors{run_counter}.DET,size(data{run_counter},2));
        
        Xn = [Regressors{run_counter}.MP det];
        Xall = [ Xn Regressors{run_counter}.Signal];
        BetaWeights = data{run_counter}*Xall*inv(Xall'*Xall);
        idxNoise    = 1:size(Xn,2);
        noi_estim   = BetaWeights(:,idxNoise) * Xall(:,idxNoise)';
        denoised_data{run_counter}   = data{run_counter} - noi_estim;
        
        if ~isempty(Regressors{run_counter}.NOISEROI)
            Xn   = [Regressors{run_counter}.NOISEROI];
            Xall = [Xn Regressors{run_counter}.Signal];
            BetaWeights = denoised_data{run_counter}*Xall*inv(Xall'*Xall);
            idxNoise    = 1:size(Xn,2);
            noi_estim   = BetaWeights(:,idxNoise) * Xall(:,idxNoise)';
            denoised_data{run_counter}   = denoised_data{run_counter} - noi_estim;
        end
           
        
        if ~isempty(Regressors{run_counter}.GSPC1)
            Xn   = [Regressors{run_counter}.GSPC1];
            Xall = [Xn Regressors{run_counter}.Signal];
            BetaWeights = denoised_data{run_counter}*Xall*inv(Xall'*Xall);
            idxNoise    = 1:size(Xn,2);
            noi_estim   = BetaWeights(:,idxNoise) * Xall(:,idxNoise)';
            denoised_data{run_counter}   = denoised_data{run_counter} - noi_estim;
        end
        
        
        if ~isempty(Regressors{run_counter}.PHYPLUS)
            Xall   = [Regressors{run_counter}.PHYPLUS];
            BetaWeights     = denoised_data{run_counter}*Xall*inv(Xall'*Xall);
            noi_estim       = BetaWeights * Xall';
            denoised_data{run_counter}   = denoised_data{run_counter} - noi_estim;
        end
    end
end