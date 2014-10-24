function [Resampling_Index] = generate_split_half_list(num_subject,num_resample)

divider = min(16,num_subject);
step = ceil(num_subject/divider);
list_subject = 1:num_subject;
A = [];
for i = 1:step-1
    xnum_subject = list_subject((i-1)*divider+1:min(i*divider,num_subject));
    A  =   [A nchoosek(xnum_subject,round(length(xnum_subject)/2))];
end
xnum_subject = list_subject((step-1)*divider+1:min(step*divider,num_subject));
temp  = nchoosek(xnum_subject,round(length(xnum_subject)/2));
if ~isempty(A)
    temp2 = repmat(temp,ceil(size(A,1)/size(temp,1)),1);
    A = [A temp2(1:size(A,1),:)];
else
    A = temp;
end
if ((num_subject/2) == fix(num_subject/2))
    index1     = randperm(fix(size(A,1)/2));
else
    if (num_resample<fix(size(A,1)/2))
        index1 = randperm(fix(size(A,1)/2));
    else
        index1 = randperm(size(A,1));
    end
end
num_resample = min(num_resample,length(index1));
index1 = index1(1:num_resample);
Resampling_Index     = A(index1,:);
