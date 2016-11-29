function QC_wrapper(step,  inputfile, analysis_model, contrast_list_str, reference_file, WARP_TYPE, newmaskname, Npcs )

if nargin < 5
	error('too few args!! must specify: step,  inputfile, analysis_model, contrast_list_str')
end


if nargin < 5
	WARP_TYPE = 'affine';
    newmaskname = [];
    Npcs = 10;
end

if strcmpi(newmaskname,'None')
    newmaskname = [];
end
if strcmpi(Npcs,'None')
    Npcs = 10;
end

if ischar(step)
    step = str2num(step);
end
if ischar(Npcs)
    Npcs = str2num(Npcs);
end

if step==0 || step==1
    Pipeline_QC1( inputfile, analysis_model, contrast_list_str )
end
if step==0 || step==2
    Pipeline_QC2(  inputfile, analysis_model, contrast_list_str, reference_file, WARP_TYPE, newmaskname, Npcs )
end

