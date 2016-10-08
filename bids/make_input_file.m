function make_input_file( filename, argcell, TASK )

% 1 IN path
% 2 OUT path
% 3 STRUCT path
% 4 PHYSIO path
% 5 DROP value1
% 6 DROP value2

fieldcell = {'IN','OUT','STRUCT','PHYSIO'};
fin = fopen(filename,'wt');

nsub=numel( TASK    );
nrun=numel( TASK{1} );

for(i=1:nsub)
for(j=1:nrun)

    for(k=1:4)
        if(strcmpi(fieldcell{k},'PHYSIO') && isempty(argcell{k}))
        disp('skipping physio.');  
        else
        fprintf(fin,'%s=%s ',fieldcell{k},argcell{k}{i}{j});
        end
    end
    fprintf(fin,'DROP=[%u,%u] TASK=%s\n',argcell{5},argcell{6},TASK{i}{j}); 

end
end

fclose(fin);
