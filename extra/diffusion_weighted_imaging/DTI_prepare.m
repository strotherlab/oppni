function DTI_prepare( InputStruct, EC_type )

% parse input structure if not already done
if( ~isstruct(InputStruct) )
    [InputStruct] = Read_Input_DTI(InputStruct);
end

disp('now running dti preparation...');

if(strcmpi(EC_type,'eddy_correct'))

    % iterate through all subjects/runs
    for ksub = 1:numel(InputStruct)

            outname = [InputStruct(ksub).run(1).Output_nifti_file_path,'/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix];
            mkdir_r(outname)

            %% ---- STANDARD DTI ACQUISITION, PROCESSING AND ANALYSIS ---- %%

            inname_30dir = [InputStruct(ksub).run(1).Input_nifti_file_path,'/',InputStruct(ksub).run(1).Input_nifti_file_prefix{1}];

            disp('DTI image analysis');

            % recopying data data
            unix(['cp ',inname_30dir, '.nii ', outname, '/DTI_rawdat.nii']);
            unix(['cp ',inname_30dir, '.bval ', outname, '/DTI_ec.bval']);
            unix(['cp ',inname_30dir, '.bvec ', outname, '/DTI_raw.bvec']);

            if( ~exist([outname, '/DTI_ec.nii'],'file' ) && ~exist([outname, '/DTI_ec.nii.gz'],'file') )
                disp('running motion+eddy correction...');
                % run eddy corrections
                unix(['eddy_correct ', outname, '/DTI_rawdat.nii.gz ', outname, '/DTI_ec.nii 0']); % create eddy corrected
                % rotate bvec directions as well (this overwrites the unrotated bvec)
                unix(['fdt_rotate_bvecs ', outname,'/DTI_raw.bvec ',outname,'/DTI_ec.bvec ',outname,'/DTI_ec.ecclog']); % update + duplicated
                disp('done corrections.');
            else
                disp('skipping EC');
            end

            if( ~exist([outname,'/DTI_fit_FA.nii'],'file' ) && ~exist([outname,'/DTI_fit_FA.nii.gz'],'file' ) )
                disp('running DTI fitting...');
                % quick brain mask
                unix(['bet2 ',outname,'/DTI_ec.nii.gz ',outname,'/DTI_bet -m']); % create mask
                % running dti parameter fits
                unix(['dtifit -k ',outname,'/DTI_ec.nii.gz -o ',outname,'/DTI_fit -m ',outname,'/DTI_bet_mask.nii.gz -r ',outname,'/DTI_ec.bvec -b ',outname,'/DTI_ec.bval']); %create output parameters            
            else
                disp('skipping DTI param. fitting 1');
            end
            % unzipping files...
            if( exist([outname,'/DTI_fit_FA.nii.gz'],'file' ) || exist([outname,'/DTI_fit_S0.nii.gz'],'file' ) )
               unix(['gunzip ', outname,'/DTI_fit_*.nii.gz']); 
            end

            %-- generating additional diffusivity parameters
            if( ~exist([outname,'/DTI_fit_RDx.nii'],'file' ) && ~exist([outname,'/DTI_fit_RDx.nii.gz'],'file' ) )

                disp('generating additional DTI measures...');

                L1=load_untouch_nii([outname,'/DTI_fit_L1.nii']);
                L2=load_untouch_nii([outname,'/DTI_fit_L2.nii']);
                L3=load_untouch_nii([outname,'/DTI_fit_L3.nii']);
                VV=L1;
                % additional measures of axial diffusivity and radial diffusivity
                VV.img = L1.img;
                save_untouch_nii(VV,[outname,'/DTI_fit_ADx.nii']);
                VV.img = (L2.img + L3.img)./2;
                save_untouch_nii(VV,[outname,'/DTI_fit_RDx.nii']);
            else
                disp('skipping DTI param. fitting 2');
            end

    %%
            if( InputStruct(ksub).run(1).MULTI_RUN_NODDI )

                disp('Multiple DTI acqusitions...preparing for (optional) noddi analysis');

                % path to dti data
                dti_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix];   
                % noddi output directory
                mkdir_r([dti_path,'/noddi_out']);

                N_runs = length(InputStruct(ksub).run(1).Input_nifti_file_prefix);

                % collecting catlist for functional data
                catlist = [];
                for(i=1:N_runs)
                    inname_multi_list{i} = [InputStruct(ksub).run(1).Input_nifti_file_path,'/',InputStruct(ksub).run(1).Input_nifti_file_prefix{i}];
                    catlist = [catlist, ' ', inname_multi_list{i},'.nii'];
                end

                % --- checking on concatt'd data
                if ~exist([outname, '/noddi_out/DTI_Multi_cat.nii.gz'],'file') || ~exist([outname,'/noddi_out/DTI_Multi_ec.bval'],'file')  || ~exist([outname,'/noddi_out/DTI_Multi_cat.bvec'],'file')       

                    disp('concatenating data...');
                    % concat functional data
                    unix(['fslmerge -t ', outname, '/noddi_out/DTI_Multi_cat.nii ', catlist]);

                    % concat the bval files now (automatically "EC'd"
                    bvcat=cell(1);
                    for(i=1:N_runs)
                        fidbv  = fopen([inname_multi_list{i}, '.bval'],'rt'); % open file
                        tline  = fgetl(fidbv); kq=0; % get 1st line
                        while ischar(tline) 
                            kq=kq+1;
                            bvcat{kq} = [bvcat{kq} tline]; % read in lines, concat to cells (resp each. line)
                            tline     = fgetl(fidbv);
                        end
                        fclose(fidbv);
                    end
                    % store to concatt'd file
                    fidOUT = fopen([outname,'/noddi_out/DTI_Multi_ec.bval'],'wt'); % create bval
                    for(i=1:length(bvcat))
                        fprintf(fidOUT, '%s\n', bvcat{i});
                    end
                    fclose(fidOUT);

                    % concat the bvec files now (need to EC)
                    bvcat=cell(3);
                    for(i=1:N_runs) 
                        fidbv  = fopen([inname_multi_list{i}, '.bvec'],'rt'); % open file
                        tline  = fgetl(fidbv); kq=0; % get 1st line
                        while ischar(tline) 
                            kq=kq+1;
                            bvcat{kq} = [bvcat{kq} tline]; % read in lines, concat to cells (resp each. line)
                            tline     = fgetl(fidbv);
                        end
                        fclose(fidbv);
                    end            
                    % store to concatt'd file
                    fidOUT = fopen([outname,'/noddi_out/DTI_Multi_cat.bvec'],'wt'); % create bval
                    for(i=1:length(bvcat))
                        fprintf(fidOUT, '%s\n', bvcat{i});
                    end
                    fclose(fidOUT);

                end

                % --- checking on EC'd data
                if ~exist([outname, '/noddi_out/DTI_Multi_ec.nii'],'file') && ~exist([outname, '/noddi_out/DTI_Multi_ec.nii.gz'],'file')      

                    disp('running motion+eddy correction...');
                    % run eddy corrections
                    unix(['eddy_correct ', outname, '/noddi_out/DTI_Multi_cat.nii.gz ', outname, '/noddi_out/DTI_Multi_ec.nii 0']); % create eddy corrected
                    % rotate bvec directions as well (this overwrites the unrotated bvec)
                    unix(['fdt_rotate_bvecs ', outname,'/noddi_out/DTI_Multi_cat.bvec ',outname,'/noddi_out/DTI_Multi_ec.bvec ',outname,'/noddi_out/DTI_Multi_ec.ecclog']); % update + duplicated
                    disp('done corrections.');
                end
                % catch when multiple files
                if( exist([outname, '/noddi_out/DTI_Multi_ec.nii'],'file') && exist([outname, '/noddi_out/DTI_Multi_ec.nii.gz'],'file') )
                    unix(['rm ',outname, '/noddi_out/DTI_Multi_ec.nii']);
                end
            else
                disp('single DTI run, initial processing complete.');
            end 
    end

elseif(strcmpi(EC_type,'eddy'))

    % iterate through all subjects/runs
    for ksub = 1:numel(InputStruct)

            outname = [InputStruct(ksub).run(1).Output_nifti_file_path,'/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix];
            mkdir_r(outname)

            %% ---- STANDARD DTI ACQUISITION, PROCESSING AND ANALYSIS ---- %%

            inname_30dir = [InputStruct(ksub).run(1).Input_nifti_file_path,'/',InputStruct(ksub).run(1).Input_nifti_file_prefix{1}];

            disp('DTI image analysis');

            % recopying data data
            unix(['cp ',inname_30dir, '.nii ', outname, '/DTI_rawdat.nii']);
            unix(['cp ',inname_30dir, '.bval ', outname, '/DTI_eddy.bval']);
            unix(['cp ',inname_30dir, '.bvec ', outname, '/DTI_raw.bvec']);

            if( ~exist([outname, '/DTI_eddy.nii'],'file' ) && ~exist([outname, '/DTI_eddy.nii.gz'],'file') )
                disp('running motion+eddy correction...');
                % premasking
                unix(['bet2 ',outname,'/DTI_rawdat.nii ',outname,'/DTI_preEddyBet -m']); % create mask
                % run eddy corrections
                unix(['eddy --imain=', outname, '/DTI_rawdat.nii --mask=',outname,'/DTI_preEddyBet_mask.nii.gz --index=',InputStruct(ksub).run(1).idx_file{1},' --acqp=',InputStruct(ksub).run(1).acq_file{1},' --bvecs=',outname,'/DTI_raw.bvec --bvals=',outname, '/DTI_eddy.bval --out=',outname,'/DTI_eddy.nii']); % create eddy corrected
                % delete rough masks
                unix(['rm ',outname,'/DTI_preEddyBet.nii.gz']);
                unix(['rm ',outname,'/DTI_preEddyBet_mask.nii.gz']);
                disp('done corrections.');
            else
                disp('skipping Eddy');
            end

            if( ~exist([outname,'/DTI_fit_FA.nii'],'file' ) && ~exist([outname,'/DTI_fit_FA.nii.gz'],'file' ) )
                disp('running DTI fitting...');
                % quick brain mask
                unix(['bet2 ',outname,'/DTI_eddy.nii.gz ',outname,'/DTI_bet -m']); % create mask
                % running dti parameter fits
                unix(['dtifit -k ',outname,'/DTI_eddy.nii.gz -o ',outname,'/DTI_fit -m ',outname,'/DTI_bet_mask.nii.gz -r ',outname,'/DTI_eddy.nii.eddy_rotated_bvecs -b ',outname,'/DTI_eddy.bval']); %create output parameters            
            else
                disp('skipping DTI param. fitting 1');
            end
            % unzipping files...
            if( exist([outname,'/DTI_fit_FA.nii.gz'],'file' ) || exist([outname,'/DTI_fit_S0.nii.gz'],'file' ) )
               unix(['gunzip ', outname,'/DTI_fit_*.nii.gz']); 
            end

            %-- generating additional diffusivity parameters
            if( ~exist([outname,'/DTI_fit_RDx.nii'],'file' ) && ~exist([outname,'/DTI_fit_RDx.nii.gz'],'file' ) )

                disp('generating additional DTI measures...');

                L1=load_untouch_nii([outname,'/DTI_fit_L1.nii']);
                L2=load_untouch_nii([outname,'/DTI_fit_L2.nii']);
                L3=load_untouch_nii([outname,'/DTI_fit_L3.nii']);
                VV=L1;
                % additional measures of axial diffusivity and radial diffusivity
                VV.img = L1.img;
                save_untouch_nii(VV,[outname,'/DTI_fit_ADx.nii']);
                VV.img = (L2.img + L3.img)./2;
                save_untouch_nii(VV,[outname,'/DTI_fit_RDx.nii']);
            else
                disp('skipping DTI param. fitting 2');
            end

    %%
            if( InputStruct(ksub).run(1).MULTI_RUN_NODDI )

                disp('Multiple DTI acqusitions...preparing for (optional) noddi analysis');

                % path to dti data
                dti_path = [InputStruct(ksub).run(1).Output_nifti_file_path '/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix];   
                % noddi output directory
                mkdir_r([dti_path,'/noddi_out']);

                N_runs = length(InputStruct(ksub).run(1).Input_nifti_file_prefix);

                % collecting catlist for functional data
                catlist = [];
                for(i=1:N_runs)
                    inname_multi_list{i} = [InputStruct(ksub).run(1).Input_nifti_file_path,'/',InputStruct(ksub).run(1).Input_nifti_file_prefix{i}];
                    catlist = [catlist, ' ', inname_multi_list{i},'.nii'];
                end

                % --- checking on concatt'd data
                if ~exist([outname, '/noddi_out/DTI_Multi_cat.nii.gz'],'file') || ~exist([outname,'/noddi_out/DTI_Multi_eddy.bval'],'file')  || ~exist([outname,'/noddi_out/DTI_Multi_cat.bvec'],'file')       

                    disp('concatenating data...');
                    % concat functional data
                    unix(['fslmerge -t ', outname, '/noddi_out/DTI_Multi_cat.nii ', catlist]);

                    % concat the bval files now (automatically "EC'd"
                    bvcat=cell(1);
                    for(i=1:N_runs)
                        fidbv  = fopen([inname_multi_list{i}, '.bval'],'rt'); % open file
                        tline  = fgetl(fidbv); kq=0; % get 1st line
                        while ischar(tline) 
                            kq=kq+1;
                            bvcat{kq} = [bvcat{kq} tline]; % read in lines, concat to cells (resp each. line)
                            tline     = fgetl(fidbv);
                        end
                        fclose(fidbv);
                    end
                    % store to concatt'd file
                    fidOUT = fopen([outname,'/noddi_out/DTI_Multi_eddy.bval'],'wt'); % create bval
                    for(i=1:length(bvcat))
                        fprintf(fidOUT, '%s\n', bvcat{i});
                    end
                    fclose(fidOUT);

                    % concat the bvec files now (need to EC)
                    bvcat=cell(3);
                    for(i=1:N_runs) 
                        fidbv  = fopen([inname_multi_list{i}, '.bvec'],'rt'); % open file
                        tline  = fgetl(fidbv); kq=0; % get 1st line
                        while ischar(tline) 
                            kq=kq+1;
                            bvcat{kq} = [bvcat{kq} tline]; % read in lines, concat to cells (resp each. line)
                            tline     = fgetl(fidbv);
                        end
                        fclose(fidbv);
                    end            
                    % store to concatt'd file
                    fidOUT = fopen([outname,'/noddi_out/DTI_Multi_cat.bvec'],'wt'); % create bval
                    for(i=1:length(bvcat))
                        fprintf(fidOUT, '%s\n', bvcat{i});
                    end
                    fclose(fidOUT);

                end

                % --- checking on EC'd data
                if ~exist([outname, '/noddi_out/DTI_Multi_eddy.nii'],'file') && ~exist([outname, '/noddi_out/DTI_Multi_eddy.nii.gz'],'file')      

                    disp('running motion+eddy correction...');
                    % premasking
                    unix(['bet2 ',outname,'/noddi_out/DTI_Multi_cat.nii.gz ',outname,'/noddi_out/DTI_preEddyBet -m']); % create mask
                    % run eddy corrections
                    unix(['eddy --imain=', outname, '/noddi_out/DTI_Multi_cat.nii.gz --mask=',outname,'/noddi_out/DTI_preEddyBet_mask.nii.gz --index=',InputStruct(ksub).run(1).idx_file{2},' --acqp=',InputStruct(ksub).run(1).acq_file{2},' --bvecs=',outname,'/noddi_out/DTI_Multi_cat.bvec --bvals=',outname, '/noddi_out/DTI_Multi_eddy.bval --out=',outname,'/noddi_out/DTI_Multi_eddy.nii']); % create eddy corrected
                    % delete rough masks
                    unix(['rm ',outname,'/noddi_out/DTI_preEddyBet.nii.gz']);
                    unix(['rm ',outname,'/noddi_out/DTI_preEddyBet_mask.nii.gz']);
                    disp('done corrections.');
                end
                
                % catch when multiple files
                if( exist([outname, '/noddi_out/DTI_Multi_eddy.nii'],'file') && exist([outname, '/noddi_out/DTI_Multi_eddy.nii.gz'],'file') )
                    unix(['rm ',outname, '/noddi_out/DTI_Multi_eddy.nii']);
                end
            else
                disp('single DTI run, initial processing complete.');
            end 
    end

else
    error('unrecognized eddy correction method');
end
