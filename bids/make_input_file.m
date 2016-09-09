function make_input_file( filename, argcell, TASK )

% 1 IN path
% 2 OUT path
% 3 STRUCT path
% 4 PHYSIO path
% 5 DROP value1
% 6 DROP value2

fieldcell = {'IN','OUT','STRUCT','PHYSIO'};
fin = fopen(filename,'wt');
for(i=1:4)
    if(strcmpi(fieldcell{i},'PHYSIO') && strcmpi(argcell{i},'None'))
    disp('skipping physio.');  
    else
    fprintf(fin,'%s=%s ',fieldcell{i},argcell{i});
    end
end

fprintf(fin,'DROP=[%u,%u] TASK=%s\n',argcell{5},argcell{6},TASK); 

fclose(fin);
