function output = module_MultiClassLDA( datamat, split_info )
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

Classification_Number = size(split_info.idx_task,1);

output.metrics.R = [];
output.metrics.P  = [];
output.metrics.Dneg = [];
output.temp.CV_alt = [];
output.temp.CV_alt_varfract = [];
output.images  = [];
for class_count = 1:Classification_Number

    split_info_out.idx_task1_sp1 = split_info.idx_task(class_count,1).sp1;
    split_info_out.idx_task1_sp2 = split_info.idx_task(class_count,1).sp2;
    split_info_out.idx_task2_sp1 = split_info.idx_task(class_count,2).sp1;
    split_info_out.idx_task2_sp2 = split_info.idx_task(class_count,2).sp2;
    
    split_info_out.drf = split_info.drf;
    split_info_out.type =  split_info.type;
    
    output_class = module_LDA( datamat, split_info_out );
    
    output.metrics.R            = [output.metrics.R output_class.metrics.R];
    output.metrics.P            = [output.metrics.P output_class.metrics.P];
    output.metrics.Dneg         = [output.metrics.Dneg output_class.metrics.Dneg];
    output.temp.CV_alt          = [output.temp.CV_alt output_class.temp.CV_alt];
    output.temp.CV_alt_varfract = [output.temp.CV_alt_varfract output_class.temp.CV_alt_varfract];
    output.images  = [output.images output_class.images];
end
output.metrics.Dneg_avg = - sqrt((1-mean(output.metrics.R)).^2  + (1-mean(output.metrics.P)).^2);