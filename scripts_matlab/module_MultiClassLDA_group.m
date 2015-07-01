function output = module_MultiClassLDA_group( datamat, split_info, Resampling_Index )
%
% =========================================================================
% MODULE_LDA: module that performs linear discriminant analysis in
% split-half NPAIRS framework, given 2 task blocks per condition.
% =========================================================================
%
%   Syntax:
%           output = module_LDA( datamat, split_info )
%
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

Subject_Num           = length(split_info);
Classification_Number = size(split_info{1}.idx_cond,1);

if( ~isfield(split_info,'drf') || isempty(split_info.drf) )
    disp('LDA uses default data reduction drf=0.5');
    for(subject_counter = 1:Subject_Num)
    split_info{subject_counter}.drf = 0.5;
    end
end

output.metrics.R   = [];
output.metrics.P   = [];
output.metrics.dPR = [];
output.temp.CV_alt = [];
output.temp.CV_alt_varfract = [];
output.images  = [];
for class_count = 1:Classification_Number
    
    for subject_counter = 1:Subject_Num
        split_info_out{subject_counter}.idx_cond1 = [split_info{subject_counter}.idx_cond(class_count,1).sp1 split_info{subject_counter}.idx_cond(class_count,1).sp2];
        split_info_out{subject_counter}.idx_cond2 = [split_info{subject_counter}.idx_cond(class_count,2).sp1 split_info{subject_counter}.idx_cond(class_count,2).sp2];
        split_info_out{subject_counter}.drf  = split_info{subject_counter}.drf;
        split_info_out{subject_counter}.type =  split_info{subject_counter}.type;
    end
    output_class = module_LDA_group( datamat, split_info_out, Resampling_Index);
    
    output.metrics.R            = [output.metrics.R output_class.metrics.R];
    output.metrics.P            = [output.metrics.P output_class.metrics.P];
    output.metrics.dPR          = [output.metrics.dPR output_class.metrics.dPR];
    output.temp.CV_alt          = [output.temp.CV_alt output_class.temp.CV_alt];
    output.temp.CV_alt_varfract = [output.temp.CV_alt_varfract output_class.temp.CV_alt_varfract];
    output.images  = [output.images output_class.images];
end
output.metrics.dPR_avg = - sqrt((1-mean(output.metrics.R)).^2  + (1-mean(output.metrics.P)).^2);
