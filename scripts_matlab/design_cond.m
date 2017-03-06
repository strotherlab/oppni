function [cond_struc] = design_cond(data,Regressors)
%
% Syntax:
%         denoised_data = apply_glm(data,Regressors)
%
% .internal script that takes 2D matrix "data" and regresses out
% "Regressors" from data. Returns "denoised_data" matrix
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

if ~iscell(data)  % use when we are running single subject analysis splitting in half

    %%%
    nsz(1) = ceil( size(data,2)/2 );
    nsz(2) = size(data,2)-ceil( size(data,2)/2 );

    if ~isfield(Regressors,'DET')
        Regressors_tmp{1}.DET = 0;
        Regressors_tmp{2}.DET = 0;
    else
        Regressors_tmp{1}.DET = Regressors.DET;
        Regressors_tmp{2}.DET = Regressors.DET;
    end
    if ~isfield(Regressors,'MP') || isempty(Regressors.MP)
        Regressors_tmp{1}.MP = [];
        Regressors_tmp{2}.MP = [];
    else
        xx1 = Regressors.MP(1:nsz(1),:);
        xx2 = Regressors.MP(nsz(1)+1:end,:);
        %%%
        Regressors_tmp{1}.MP = xx1(:,sum(abs(xx1),1)>0);
        Regressors_tmp{2}.MP = xx2(:,sum(abs(xx2),1)>0);
    end
    if ~isfield(Regressors,'GSPC1') || isempty(Regressors.GSPC1)
        Regressors_tmp{1}.GSPC1 = [];
        Regressors_tmp{2}.GSPC1 = [];
    else
        xx1 = Regressors.GSPC1(1:nsz(1),:);
        xx2 = Regressors.GSPC1(nsz(1)+1:end,:);
        %%%
        Regressors_tmp{1}.GSPC1 = xx1(:,sum(abs(xx1),1)>0);
        Regressors_tmp{2}.GSPC1 = xx2(:,sum(abs(xx2),1)>0);
    end
    if ~isfield(Regressors,'PHYPLUS') || isempty(Regressors.PHYPLUS)
        Regressors_tmp{1}.PHYPLUS = [];
        Regressors_tmp{2}.PHYPLUS = [];
    else
        xx1 = Regressors.PHYPLUS(1:nsz(1),:);
        xx2 = Regressors.PHYPLUS(nsz(1)+1:end,:);
        %%%
        Regressors_tmp{1}.PHYPLUS = xx1(:,sum(abs(xx1),1)>0);
        Regressors_tmp{2}.PHYPLUS = xx2(:,sum(abs(xx2),1)>0);
    end
    if ~isfield(Regressors,'Signal') || isempty(Regressors.Signal)
        Regressors_tmp{1}.Signal = [];
        Regressors_tmp{2}.Signal = [];
    else
        xx1 = Regressors.Signal(1:nsz(1),:);
        xx2 = Regressors.Signal(nsz(1)+1:end,:);
        %%%
        Regressors_tmp{1}.Signal = xx1(:,sum(abs(xx1),1)>0);
        Regressors_tmp{2}.Signal = xx2(:,sum(abs(xx2),1)>0);
    end
    if ~isfield(Regressors,'NOISEROI') || isempty(Regressors.NOISEROI)
        Regressors_tmp{1}.NOISEROI = [];
        Regressors_tmp{2}.NOISEROI = [];
    else
        xx1 = Regressors.NOISEROI(1:nsz(1),:);
        xx2 = Regressors.NOISEROI(nsz(1)+1:end,:);
        %%%
        Regressors_tmp{1}.NOISEROI = xx1(:,sum(abs(xx1),1)>0);
        Regressors_tmp{2}.NOISEROI = xx2(:,sum(abs(xx2),1)>0);
    end
    
    Regressors = Regressors_tmp; clear Regressors_tmp;
    
    for(run_counter = 1:2 )
        det = get_legendre(Regressors{run_counter}.DET,nsz(run_counter));

        Xall = [det, ...
                Regressors{run_counter}.MP,... 
                Regressors{run_counter}.GSPC1,...
                Regressors{run_counter}.PHYPLUS,...
                Regressors{run_counter}.Signal,...
                Regressors{run_counter}.NOISEROI];
        Xall = bsxfun(@rdivide,Xall,sqrt(sum(Xall.^2)));
        
        cond_struc.XCond(run_counter,1)     = cond( Xall );
        cond_struc.X_num_rank(run_counter,1) = size(Xall,2);
        cond_struc.X_num_rank(run_counter,2) = rank(Xall  );
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

        Xall = [det, ...
                Regressors{run_counter}.MP,... 
                Regressors{run_counter}.GSPC1,...
                Regressors{run_counter}.PHYPLUS,...
                Regressors{run_counter}.Signal,...
                Regressors{run_counter}.NOISEROI];
        Xall = bsxfun(@rdivide,Xall,sqrt(sum(Xall.^2)));
            
        cond_struc.XCond(run_counter,1)     = cond( Xall );
        cond_struc.X_num_rank(run_counter,1) = size(Xall,2);
        cond_struc.X_num_rank(run_counter,2) = rank(Xall  );
    end
end