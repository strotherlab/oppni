function DTI_prepare( InputStruct )

% parse input structure if not already done
if( ~isstruct(InputStruct) )
    [InputStruct] = Read_Input_DTI(InputStruct);
end

disp('now running dti preparation...');

% iterate through all subjects/runs
for ksub = 1:numel(InputStruct)

        outname = [InputStruct(ksub).run(1).Output_nifti_file_path,'/dti_processed/',InputStruct(ksub).run(1).Output_nifti_file_prefix];
        mkdir_r(outname)

%%      OPTION 1: multi-run DWI acquisition --> concatenate for HARDI analysis
        if( InputStruct(ksub).run(1).MULTI_RUN_HARDI )
            
            inname_30dir = [InputStruct(ksub).run(1).Input_nifti_file_path,'/',InputStruct(ksub).run(1).Input_nifti_file_prefix{1}];
            inname_64dir = [InputStruct(ksub).run(1).Input_nifti_file_path,'/',InputStruct(ksub).run(1).Input_nifti_file_prefix{2}];
            
            disp('HARDI image analysis.');
            disp('concatenating data...');

            % concat functional data
            unix(['fslmerge -t ', outname, '/HARDI_concat.nii ', inname_30dir, '.nii ', inname_64dir,'.nii']);

            % concat bval
            fid30  = fopen([inname_30dir, '.bval'],'rt');
            fid64  = fopen([inname_64dir, '.bval'],'rt');
            fidOUT = fopen([outname,'/HARDI_ec.bval'],'wt'); % create bval
            tline30 = fgetl(fid30);
            tline64 = fgetl(fid64);
            while ischar(tline30) 
                tline30_64 = [tline30 tline64];
                fprintf(fidOUT, '%s\n', tline30_64);
                tline30 = fgetl(fid30);
                tline64 = fgetl(fid64);
            end
            fclose(fid30); 
            fclose(fid64);
            fclose(fidOUT);

            % concat bvec
            fid30  = fopen([inname_30dir, '.bvec'],'rt');
            fid64  = fopen([inname_64dir, '.bvec'],'rt');
            fidOUT = fopen([outname,'/HARDI_cat.bvec'],'wt'); % create bvec
            tline30 = fgetl(fid30);
            tline64 = fgetl(fid64);
            while ischar(tline30) 
                tline30_64 = [tline30 tline64];
                fprintf(fidOUT, '%s\n', tline30_64);
                tline30 = fgetl(fid30);
                tline64 = fgetl(fid64);
            end
            fclose(fid30); 
            fclose(fid64);
            fclose(fidOUT);

            disp('done concatenating.')

            disp('running motion+eddy correction...');
            % run eddy corrections
            unix(['eddy_correct ', outname, '/HARDI_concat.nii.gz ', outname, '/HARDI_ec.nii 0']); % create eddy corrected
            % rotate bvec directions as well (this overwrites the unrotated bvec)
            unix(['fdt_rotate_bvecs ', outname,'/HARDI_cat.bvec ',outname,'/HARDI_ec.bvec ',outname,'/HARDI_ec.ecclog']); % update + duplicated
            disp('done corrections.');

            disp('running DTI fitting...');
            % quick brain mask
            unix(['bet2 ',outname,'/HARDI_ec.nii.gz ',outname,'/HARDI_bet -m']); % create mask
            % running dti parameter fits
            unix(['dtifit -k ',outname,'/HARDI_ec.nii.gz -o ',outname,'/HARDI_fit -m ',outname,'/HARDI_bet_mask.nii.gz -r ',outname,'/HARDI_ec.bvec -b ',outname,'/HARDI_ec.bval']); %create output parameters

            disp('done fitting');

%%      OPTION 2: single-run DWI acquisition --> perform standard DTI analysis
        else
            
            inname_30dir = [InputStruct(ksub).run(1).Input_nifti_file_path,'/',InputStruct(ksub).run(1).Input_nifti_file_prefix{1}];

            disp('DTI image analysis');

            % recopying data data
            unix(['cp ',inname_30dir, '.nii ', outname, '/DTI_rawdat.nii']);
            unix(['cp ',inname_30dir, '.bval ', outname, '/DTI_ec.bval']);
            unix(['cp ',inname_30dir, '.bvec ', outname, '/DTI_raw.bvec']);

            disp('done concatenating.')

            disp('running motion+eddy correction...');
            % run eddy corrections
            unix(['eddy_correct ', outname, '/DTI_rawdat.nii.gz ', outname, '/DTI_ec.nii 0']); % create eddy corrected
            % rotate bvec directions as well (this overwrites the unrotated bvec)
            unix(['fdt_rotate_bvecs ', outname,'/DTI_raw.bvec ',outname,'/DTI_ec.bvec ',outname,'/DTI_ec.ecclog']); % update + duplicated
            disp('done corrections.');

            disp('running DTI fitting...');
            % quick brain mask
            unix(['bet2 ',outname,'/DTI_ec.nii.gz ',outname,'/DTI_bet -m']); % create mask
            % running dti parameter fits
            unix(['dtifit -k ',outname,'/DTI_ec.nii.gz -o ',outname,'/DTI_fit -m ',outname,'/DTI_bet_mask.nii.gz -r ',outname,'/DTI_ec.bvec -b ',outname,'/DTI_ec.bval']); %create output parameters

            disp('done fitting');
        end

end

%%