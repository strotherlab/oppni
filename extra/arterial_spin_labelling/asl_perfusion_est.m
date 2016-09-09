function out = asl_perfusion_est( input_file, M0_ref, mask_file, subtract_type, dataInfo, format, outname )
%
% =========================================================================
% ASL_PERFUSION_EST:  code to compute pasl perfusion estimates in OPPNI 
% code framework. Adapted from "ASLtbx" package by Ze Wang
% (https://cfn.upenn.edu/~zewang/ASLtbx.php) and the associated paper.
% Any results should be carefully checked against their codeset if possible.
% =========================================================================
%

V = load_untouch_nii(input_file);

if(isempty(mask_file))
    % no mask provided
    disp('no mask provided, estimating (a liberal) one from functional data.') 
    disp('Check outputs carefully to make sure it is correct!');
    mask = double(  mean(V.img,4) > ( 0.2*max(reshape(mean(V.img,4),[],1)) )  );
    M    = V;
    M.img= mask;
else
    % load mask
    M = load_untouch_nii(mask_file);
    mask = double(M.img>0);
end
voldims = size(V.img);
for(i=1:3) 
    if( voldims(i) ~= size(M.img,i) ) error(['mask and input volume do not match on dim.',num2str(i)]); end 
end
% to 2D matrix
volmat = nifti_to_mat(V,M);

if( ischar(M0_ref) )
    disp('separate M0 file specified');
    M0vol = load_untouch_nii(M0_ref);
    m0vct = nifti_to_mat(M0vol,M);
    %
    drop  = 0;
elseif( isnumeric(M0_ref) && isfinite(M0_ref) )
    if( M0_ref > size(voldims(4)) ) error('M0_ref index exceeds size of input file'); end
    if( M0_ref < 1 )                error('M0_ref index must be positive integer');   end
    m0vct = volmat(:,M0_ref); %% assign as M0 volume
    volmat(:,M0_ref)=[]; %% drop from run
    %
    drop  = 1;
else
    error('M0_ref is something strange?');
end

Ntime = size( volmat,2 );
if( rem(Ntime,2) > 0 ) 
    error('tag and controls are not matched? Should be an even number'); edit 
end

% auto-checking on label order
delM_init = mean(volmat(:,1:2:end),2) - mean(volmat(:,2:2:end),2);
% based on median perfusion (must be >0)
if( median(delM_init)>0 ) 
    disp('autocheck: label precedes control');
    lab_ctl =  1; %%label precedes control
    labix   =  1:2:Ntime;
    conix   =  2:2:Ntime;
else
    disp('autocheck: control precedes label');
    lab_ctl = -1; %%control precedes label
    conix   =  1:2:Ntime;
    labix   =  2:2:Ntime;
end

% Now setting some parameters for CBF estimation
% model: TCBF = lambda*PERF*6000*1000./(2*m0vct.*exp(-TI/BloodT1)*TI1*labeff*qTI);
%
% pre-initializing par. values
%
par.lambda =  0.9; % - blood/tissue water partition coefficient %*100*60;   %0.9 mL/g
par.BloodT1= 1664; % - T1 for blood (msec); ref Lu 04 and Cavusoglu 09 MRI. This is B0 dependent! (BloodT1=1200 at 1.5T...originally BloodT1=1490 at 3.0T / Wang 03)
par.TR     = 2500; % - repetition time. This is default...may need to adjust
par.TI1    =  700; % - Tagging bolus duration in ms, for the QUIPSS II. This is the default value...may need to adjust
par.TI2    = 1800; % - second inversion time; delay time for labeled spin to enter the imaging slice. This is the default value...may need to adjust
par.labeff = 0.95; % - labeling efficiency, 0.95 for PASL, 0.68 for CASL, 0.85 for pCASL, this should be measured for onsite scanner
par.qTI    = 0.85; % - close to unit, and is set to 0.85 in Warmuth 05

% now modifying according to user
if( (nargin >= 5) && ~isempty(dataInfo) )
   % current available list of parameters to modify
   paramlist = {'lambda','BloodT1','TI1','labeff','qTI'};
   % copying over pre-specified values
   for(i=1:length(paramlist))
        if(isfield(dataInfo,paramlist{i}))
            disp(['updating ',paramlist{i}])
            % updating par. structure
            par.(paramlist{i}) = dataInfo.(paramlist{i});
        end
   end
end

% Slice timing array
%
% offset in acq. per slice (min-TR - labeltime - delaytime)/#slices
Slicetime   = (par.TR - par.TI1 - par.TI2)/voldims(3);
% get per-slice adjustments
tmp         = bsxfun(@times, mask, permute(0:voldims(3)-1,[3 1 2]) ); %% scale voxel values by slice order (start at zero)
slc_timevct = tmp(mask>0) .* Slicetime; %% vectorized, scaled by slice timing


%% COMPUTING PERFUSION WEIGHTED IMAGES + BOLD series

if    ( (isnumeric(subtract_type) && (subtract_type==1)) || strcmpi( subtract_type,'simple') )

    % subtract nearest control from each label
    PERF = volmat(:,labix) - volmat(:,conix);
    % BOLD: simple average of tag-control
    BOLD = (volmat(:,labix) + volmat(:,conix))./2;
    
elseif( (isnumeric(subtract_type) && (subtract_type==2)) || strcmpi( subtract_type,'surround') )
    
    if    ( lab_ctl>0 ) 
        % avg of matched(subseq) and prev
        con_mat = ( volmat(:,conix) + [volmat(:,conix(1)) volmat(:,conix(1:end-1))] )./2;
        %
    elseif( lab_ctl<0 )
        % avg of matched(prior) and subseq
        con_mat = ( volmat(:,conix) + [volmat(:,conix(2:end)) volmat(:,conix(end))] )./2;
        %
    end
    % subtract from tag
    PERF = volmat(:,labix) - con_mat;
    % BOLD: simple average of tag-control
    BOLD = (volmat(:,labix) + con_mat)./2;    
    
elseif( (isnumeric(subtract_type) && (subtract_type==3)) || strcmpi( subtract_type,'sinc') )
    
    % lab_ctl -1  same as  FirstImageType==1 (control before label)
    Timeshift = 0.5;
    for(n=1:round(Ntime/2)) %% step through volumes
       %
       % 6 point sinc interpolation
       if     lab_ctl>0
           idx=n+[-3 -2 -1 0 1 2];
           normloc=3-Timeshift;
       elseif lab_ctl<0
           idx=n+[-2 -1 0 1 2 3];
           normloc=2+Timeshift;
       end
       idx(idx<1)=1;
       idx(idx>round(Ntime/2))=round(Ntime/2);
       con_img=sinc_interpVec( volmat(:,conix(idx)) ,normloc);
       % now subtract interpolated value from tag
       PERF(:,n) = volmat(:,labix(n)) - con_img;
       % BOLD: averaging of tag-control
       BOLD(:,n) = (volmat(:,labix(n)) + con_img)./2;
    end
    
else
    error('subtract_type needs to be (1,2,3) or (simple, surround, sinc)');
end

clear delM_init; %% clear initialized

%% NOW QUANTIFICATION OF CBF

% adjusted TI2, based on TI2 (delay time) & slice-specific delay
TI_adj= (par.TI2) + slc_timevct;
% using the ASLtbx formulation of TCBF (ml/100g/ms) --> (ml/100g/min):
TCBF = 6000*1000 * bsxfun(@rdivide, (par.lambda*PERF), (2*m0vct.*exp(-TI_adj./par.BloodT1)*par.TI1*par.labeff*par.qTI) );

% Storing outputs
out.TCBF = TCBF;
out.PERF = PERF;
out.BOLD = BOLD;
% averages too
out.tcbf_mean = mean(TCBF,2);
out.perf_mean = mean(PERF,2);
out.bold_mean = mean(BOLD,2);
out.m0        = m0vct;

if( format~=0 )
    
    % define path+name of output files
    if(~isempty(outname))
         [apath aname aext] = fileparts( outname );
    else [apath aname aext] = fileparts( input_file );
         aname = [aname,'_output'];
    end
    outname = [apath,'/',aname];

    % predefine header nii
    nii=V;
    nii.hdr.dime.datatype = 16;
    nii.hdr.hist = V.hdr.hist;
    % initialize 4d tmpvol
    TMPVOL = zeros( [size(mask), round(Ntime/2)] ); 

    % --------------------- saving 4d functional data ---------------------
    nii.hdr.dime.dim(5) = size(TMPVOL,4);
    % --tcbf--
    for(t=1:round(Ntime/2)) 
        tmp=mask;tmp(tmp>0)=out.TCBF(:,t); 
        TMPVOL(:,:,:,t) = tmp; 
    end
    nii.img = TMPVOL;
    save_untouch_nii(nii,[outname,'_TCBF.nii']);
    % --perfusion--
    for(t=1:round(Ntime/2)) 
        tmp=mask;tmp(tmp>0)=out.PERF(:,t); 
        TMPVOL(:,:,:,t) = tmp; 
    end
    nii.img = TMPVOL;
    save_untouch_nii(nii,[outname,'_PERF.nii']);
    % --bold--
    for(t=1:round(Ntime/2)) 
        tmp=mask;tmp(tmp>0)=out.BOLD(:,t); 
        TMPVOL(:,:,:,t) = tmp; 
    end
    nii.img = TMPVOL;
    save_untouch_nii(nii,[outname,'_BOLD.nii']);

    % --------------------- saving average volumes ---------------------
    nii.hdr.dime.dim(5) = 1;
    % --mean:tcbf--
    TMPVOL=mask;TMPVOL(TMPVOL>0)=out.tcbf_mean;
    nii.img = TMPVOL;
    save_untouch_nii(nii,[outname,'_TCBF_avg.nii']);
    nii.img = TMPVOL;
    % --mean:perfusion--
    TMPVOL=mask;TMPVOL(TMPVOL>0)=out.perf_mean;
    nii.img = TMPVOL;
    save_untouch_nii(nii,[outname,'_PERF_avg.nii']);
    % --mean:bold--
    TMPVOL=mask;TMPVOL(TMPVOL>0)=out.bold_mean;
    nii.img = TMPVOL;
    save_untouch_nii(nii,[outname,'_BOLD_avg.nii']);
end
if( format <0 ) out=[]; end %clear "out" struct if unwanted

%%
function y = sinc_interpVec(x,u)
% sinc interpolation function, the input can be a data vector or matrix,
% the number of rows is assumed to be the number of data dimension, the
% number of columns is assumed to be the time points. The interpolation is
% applied column wise.
% Ze Wang @ 6-09-2004 Upenn
[dim,len]=size(x);
[dim2,ulen]=size(u);
if dim2==1
    u=repmat(u,dim,1);
else
    if dim2~=dim disp('We can'' figure out what you want to do\n');return; end;
end
m = 0:len-1;
m=repmat(m,[dim,1]);
for i=1:ulen
    weight=sinc(m- repmat(u(:,i),1,len));
    swei=sum(weight,2);
    if abs(swei(1)-1)>0.1
        weight=weight./repmat(swei,1,len);
    end
  
  y(:,i) = sum(x.*weight, 2);
end
