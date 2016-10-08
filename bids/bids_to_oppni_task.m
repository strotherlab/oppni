function bids_to_oppni_task(jsonfile, tsvfile, task_type, newtaskfile)

%% 1. overal design specifications

fid = fopen(jsonfile);
C_text=textscan(fid,'%s','delimiter','\n');
fclose(fid);

C_str=[];
for(i=1:length(C_text{1})) C_str=[C_str C_text{1}{i}]; end
strucj = p_json(C_str);

if( ~isfield(strucj,'RepetitionTime') ) error('need RepetitionTime field'); end
if(~isnumeric(strucj.RepetitionTime)) strucj.RepetitionTime=str2num(strucj.RepetitionTime); end

nsub=numel( newtaskfile    );
nrun=numel( newtaskfile{1} );

for(i=1:nsub)
for(j=1:nrun)

    fin = fopen(newtaskfile{i}{j},'wt');

    fprintf(fin,'TR_MSEC=[%s]\n',num2str(1000*strucj.RepetitionTime));  
    fprintf(fin,'UNIT=[sec]\n');    
    fprintf(fin,'TYPE=[%s]\n\n',task_type);    

    %% 2. task conditions, onsets, durations

    fid = fopen(tsvfile{i}{j});
    tsvfile{i}{j},
    C_text=textscan(fid,'%s','delimiter','\n');
    fclose(fid);

    headr = C_text{1}{1};
    body  = C_text{1}(2:end);
    colnames=regexp(headr,'\t','split');
    nrow = length(colnames);
    ncol = length(body);

    if( ~strcmp(colnames{1},'onset') ) error('onset should be 1st column'); end
    if( ~strcmp(colnames{2},'duration') ) error('duration should be 2nd column'); end

    ixtrial=0;
    for(k=1:nrow) if(strcmp(colnames{k},'trial_type')) ixtrial=k; end; end
    if(ixtrial==0) error('need trial_type specification'); end

    for(k=1:ncol)
       bodycell(k,:) = regexp(body{k},'\t','split'); % per col
    end
    d=bodycell(:,ixtrial);
    tasklist = unique(d);

    for(k=1:length(tasklist))

        % collect all onsets, durations affiliated with task
        onslist='';
        durlist='';
        for(l=1:size(body,1))
            if(strcmp(d{l},tasklist{k}))
                if( all(ismember(bodycell{l,1},'.0123456789 ')) &&  all(ismember(bodycell{l,2},'.0123456789 ')) )
                    if(isempty(onslist) && isempty(durlist))
                    onslist = [bodycell{l,1}];
                    durlist = [bodycell{l,2}];                        
                    else
                    onslist = [onslist ',' bodycell{l,1}];
                    durlist = [durlist ',' bodycell{l,2}];    
                    end
                end
            end
        end

        if( ~isempty(onslist) && ~isempty(durlist) )

            % print to file
            fprintf(fin,'NAME=[%s]\n',tasklist{k});    
            fprintf(fin,'ONSETS=[%s]\n',onslist);    
            fprintf(fin,'DURATION=[%s]\n\n',durlist);  
        end
    end

    fclose(fin);

end
end