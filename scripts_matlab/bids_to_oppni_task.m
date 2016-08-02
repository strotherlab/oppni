function bids_to_oppni_task(jsonfile, tsvfile, task_type, newtaskfile)

fin = fopen(newtaskfile,'wt');

%% 1. overal design specifications

fid = fopen(jsonfile);
C_text=textscan(fid,'%s','delimiter','\n');
fclose(fid);

C_str=[];
for(i=1:length(C_text{1})) C_str=[C_str C_text{1}{i}]; end
strucj = p_json(C_str);

if( ~isfield(strucj,'RepetitionTime') ) error('need RepetitionTime field'); end

if(~isnumeric(strucj.RepetitionTime)) strucj.RepetitionTime=str2num(strucj.RepetitionTime); end

fprintf(fin,'TR_MSEC=[%s]\n',num2str(1000*strucj.RepetitionTime));  
fprintf(fin,'UNIT=[sec]\n');    
fprintf(fin,'TYPE=[%s]\n\n',task_type);    

%% 2. task conditions, onsets, durations

fid = fopen(tsvfile);
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
for(i=1:nrow) if(strcmp(colnames{i},'trial_type')) ixtrial=i; end; end
if(ixtrial==0) error('need trial_type specification'); end
    
for(i=1:ncol)
   bodycell(i,:) = regexp(body{i},'\t','split'); % per col
end
d=bodycell(:,ixtrial);
tasklist = unique(d);

for(i=1:length(tasklist))
    
    % collect all onsets, durations affiliated with task
    onslist='';
    durlist='';
    for(j=1:size(body,1))
        if(strcmp(d{j},tasklist{i}))
            if( all(ismember(bodycell{j,1},'.0123456789 ')) &&  all(ismember(bodycell{j,2},'.0123456789 ')) )
                if(isempty(onslist) && isempty(durlist))
                onslist = [bodycell{j,1}];
                durlist = [bodycell{j,2}];                        
                else
                onslist = [onslist ' ' bodycell{j,1}];
                durlist = [durlist ' ' bodycell{j,2}];    
                end
            end
        end
    end
    
    if( ~isempty(onslist) && ~isempty(durlist) )

        % print to file
        fprintf(fin,'NAME=[%s]\n',tasklist{i});    
        fprintf(fin,'ONSETS=[%s]\n',onslist);    
        fprintf(fin,'DURATION=[%s]\n\n',durlist);  
    end
end

fclose(fin);
