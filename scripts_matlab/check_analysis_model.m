function analysis_model_out = check_analysis_model( analysis_model )
%
% Syntax:
%         analysis_model_out = check_analysis_model( analysis_model )
%
% .internal script checks if requested analysis module exists, returns
%  model attributes in "analysis_model_out"
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

if( strcmpi(analysis_model,'NONE') )
    %-- no analysis being performed
    analysis_model_out.model_name  = 'NONE';
    analysis_model_out.design_type = 'nocontrast';
else

    global CODE_PATH
    if isempty(CODE_PATH)
        CODE_PATH = fileparts(which('Pipeline_PART1_afni_steps.m'));
        if CODE_PATH(end)~='/'
            CODE_PATH = [CODE_PATH '/'];
            addpath(CODE_PATH);
        end
    end
    addpath([CODE_PATH '/analysis_modules'])        

    e = dir( [CODE_PATH 'analysis_modules'] );
    kq=0;
    for(i=1:length(e))
       if( ~isempty(strfind(e(i).name,'.m')) )
          kq=kq+1;
          [path module_list{kq} ext] = fileparts( e(i).name );
       end
    end

    ix = find( strcmpi( analysis_model, module_list ) );

    if( isempty(ix) )
        display(sprintf('ERROR: The analysis model %s not found among:\n',analysis_model));
        display(module_list);    
        sge_exit(100);    
    else
        analysis_model = module_list{ix};
    end

    p   = str2func(analysis_model);
    tmp = p(); %% store attributes
    analysis_model_out          = tmp.attributes;
    analysis_model_out.filepath = [CODE_PATH 'analysis_modules']; %% store actual filename

end