function output = Analysis_Combine_Multirun( out, runlist, maskname, onsetlist, durationlist, onsetcens, durationcens, TR, DROP, contrast, wmspat )
%
% Script to run analysis that combines multiple runs, in cases where
% they cannot be easily analyzed independently
%

[opath,oname,oext] = fileparts(out);

mkdir([opath, '/analysis_multirun_intermed']);

if(nargin<10) contrast=[]; end
if(nargin<11) wmspat  =[]; end

nruns = length(runlist);
ncond = length(onsetlist);

% initialized list of onsets + durations
VAL_ONSET=cell(nruns,ncond);
VAL_DURAT=cell(nruns,ncond);

for(c=1:ncond) % per condition:
    
    disp(['checking file...', onsetlist{c}]);
    
    % list of onsets
    fid    = fopen(onsetlist{c});
    tline=fgetl(fid);
    klin=0;
    while(ischar(tline))
       klin=klin+1;
       VAL_ONSET{klin,c} = str2num(tline);
       tline=fgetl(fid);
    end
    fclose(fid);
    
    %^ list of durations
    fid    = fopen(durationlist{c});
    tline=fgetl(fid);
    klin=0;
    while(ischar(tline))
       klin=klin+1;
       VAL_DURAT{klin,c} = str2num(tline);
       tline=fgetl(fid);
    end
    fclose(fid);
end

if( ~isempty(onsetcens) )
    %^ list of durations
    VAL_ONSET_C=cell(nruns,1);
    fid    = fopen(onsetcens);
    tline=fgetl(fid);
    klin=0;
    while(ischar(tline))
       klin=klin+1;
       VAL_ONSET_C{klin} = str2num(tline);
       tline=fgetl(fid);
    end
    fclose(fid);

    %^ list of durations
    VAL_DURAT_C=cell(nruns,1);
    fid    = fopen(durationcens);
    tline=fgetl(fid);
    klin=0;
    while(ischar(tline))
       klin=klin+1;
       VAL_DURAT_C{klin} = str2num(tline);
       tline=fgetl(fid);
    end
    fclose(fid);
else
    VAL_ONSET_C=[];
    VAL_DURAT_C=[];
end

%% NOW LOAD FMRI DATA

M=load_untouch_nii(maskname);
mask=double(M.img);

if( exist([opath, '/analysis_multirun_intermed/',oname,'_dat.mat'],'file') ==0 )

    disp('concatenating...');
    
    for(c=1:ncond)
    ONSET_1LIST{c} = [];
    DURAT_1LIST{c} = [];
    end
    ONSETc1LIST = [];
    DURATc1LIST = [];
    epicat      = [];

    for(n=1:nruns)
        
        disp(strcat(['run ',num2str(n),' of ', num2str(nruns)]));

        % into 2d volumes
        V=load_untouch_nii(runlist{n});
        voltmp     = nifti_to_mat( V,M );
        epicat     = [epicat voltmp];
        ntime(n,1) = size(voltmp,2);

        % adjusted onsets
        for(c=1:ncond)
        if(~isempty(VAL_ONSET{n,c})) 

            tmp1 = VAL_ONSET{n,c} - DROP(1)*TR;               % adjust time to first un-dropped scan 
            % keep negative onsets to subtract from duration
            dursub = tmp1.*double(tmp1<0);
            % zero out the negative onsets;
            tmp1(tmp1<=0)=eps;
            % put into the list
            ONSET_1LIST{c} = [ONSET_1LIST{c}, tmp1+sum(ntime(1:n-1))*TR];
            % trim durations to fit
            tmp2     = VAL_DURAT{n,c} + dursub;
            ix       = find( (tmp1 + tmp2)>(ntime(n)*TR) ); % if exceeds run length...
            tmp2(ix) = (ntime(n)*TR) - tmp1(ix);            % trim back to =run length
            DURAT_1LIST{c} = [DURAT_1LIST{c} tmp2];
        end    
        end
        if(~isempty(VAL_ONSET_C) && ~isempty(VAL_ONSET_C{n})) 

            tmp1 = VAL_ONSET_C{n} - DROP(1)*TR;               % adjust time to first un-dropped scan 
            % keep negative onsets to subtract from duration
            dursub = tmp1.*double(tmp1<0);
            % zero out the negative onsets;
            tmp1(tmp1<=0)=eps;
            % put into the list
            ONSETc1LIST = [ONSETc1LIST, tmp1+sum(ntime(1:n-1))*TR];
            % trim durations to fit
            tmp2     = VAL_DURAT_C{n} + dursub;
            ix       = find( (tmp1 + tmp2)>(ntime(n)*TR) ); % if exceeds run length...
            tmp2(ix) = (ntime(n)*TR) - tmp1(ix);            % trim back to =run length
            DURATc1LIST = [DURATc1LIST tmp2];
        end        
    end

    save([opath, '/analysis_multirun_intermed/',oname,'_dat.mat'],'epicat','ONSET_1LIST','DURAT_1LIST','ONSETc1LIST','DURATc1LIST');
end

%% RELOAD CONCATENATED DATA

load([opath, '/analysis_multirun_intermed/',oname,'_dat.mat']);

% design matrix (time x condition+cens)
design=zeros( size(epicat,2), ncond );
for(c=1:ncond)
    ONSET_TR =  ceil(ONSET_1LIST{c}./TR);
    DURAT_TR = round(DURAT_1LIST{c}./TR);
    for(j=1:length(ONSET_TR))
        design( ONSET_TR(j):ONSET_TR(j)+DURAT_TR(j), c)=1;
    end
end
if(~isempty(ONSETc1LIST))
    censor=zeros(size(epicat,2),1);
    ONSET_TR =  ceil(ONSETc1LIST./TR);
    DURAT_TR = round(DURATc1LIST./TR);    
    for(j=1:length(ONSETc1LIST))
        censor( ONSET_TR(j):ONSET_TR(j)+DURAT_TR(j), 1)=1;
    end 
else
    censor=[];
end

% truncate
design = design( 1:size(epicat,2), :);
if( ~isempty(censor) )
censor = censor( 1:size(epicat,2) );
end

figure,imagesc([design censor]);

%% CHOOSING ANALYSIS METHOD

if( isempty(contrast) ) % analyze all contrasts in GLM
    
    ntot = size(epicat,2);
    Xsignl = design_to_hrf( design,TR,[5 15]);
    
%     o1 = GLM_model_fmri( epicat(:,1:ceil(ntot/2)    ), 0, censor(1:ceil(ntot/2)    ), Xsignl(1:ceil(ntot/2),:    ), [] );
%     o2 = GLM_model_fmri( epicat(:,ceil(ntot/2)+1:end), 0, censor(ceil(ntot/2)+1:end), Xsignl(ceil(ntot/2)+1:end,:), [] );
%     
%     for(c=1:ncond)
%         [output.rep(c,1) output.spm(:,c)] = get_rSPM( o1.Beta_signl(:,c), o2.Beta_signl(:,c), 1 );
%     end

    %% WM regression
    %epicat = detrend_matrix( epicat, 0, mean(epicat(wmspat>0.5,:))',[]);

    if( ~isempty(censor))
    oo = GLM_model_fmri( epicat, 0, censor, Xsignl(1:size(epicat,2),:), [] );
    else
    oo = GLM_model_fmri( epicat, 0, [], Xsignl(1:size(epicat,2),:), [] );
    end
    output.beta = oo.Beta_signl;
    output.tmap = oo.Tmap_signl;

    % now save files
    nii = M;
    nii.img = zeros([size(M.img) size(output.tmap,2)]);
    Z = zeros(size(M.img));
    for nvol = 1:size(output.tmap,2)
        Z(M.img~=0) = output.tmap(:,nvol);
        nii.img(:,:,:,nvol) = Z;
    end
    nii.hdr.dime.datatype = 16;
    nii.hdr.dime.dim([1 5]) = [4 nvol];
    save_untouch_nii(nii,[opath,'/',oname,'.nii']);

    
end

