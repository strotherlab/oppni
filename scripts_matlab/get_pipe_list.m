function [pipeset_half, detSet, mprSet, tskSet, phySet, gsSet, lpSet, Nhalf, Nfull, censorType] = get_pipe_list( filename )

% ------------------------------------------------------------------------%
% Authors: Nathan Churchill, University of Toronto
%          email: nathan.churchill@rotman.baycrest.on.ca
%          Babak Afshin-Pour, Rotman reseach institute
%          email: bafshinpour@research.baycrest.org
% ------------------------------------------------------------------------%
% CODE_VERSION = '$Revision: 158 $';
% CODE_DATE    = '$Date: 2014-12-02 18:11:11 -0500 (Tue, 02 Dec 2014) $';
% ------------------------------------------------------------------------%
% LMP (2021-02-09) - modified to allow censorType parameter from pipeline file
% 
% reads in the inputfile
fid = fopen(filename);
newline = fgetl(fid);
newline(newline==' ') = [];
pipelinelist=[];
while ischar(newline)
    %
    pipelinelist = [pipelinelist ' ' newline];
    newline      = fgetl(fid);
    newline(newline==' ') = [];
end
fclose(fid);
substr = cell(11,1);
proclist = {'MOTCOR='; 'CENSOR='; 'RETROICOR='; 'TIMECOR=';'SMOOTH='; 'DETREND=';'MOTREG=';'TASK=';'PHYPLUS='; 'GSPC1='; 'CUSTOMREG='; 'LOWPASS'};
ileft  = strfind( pipelinelist, '[' );
iright = strfind( pipelinelist, ']' );
for(s=1:length(proclist))
    iStep = strfind( upper(pipelinelist), proclist{s} );
    if isempty(iStep)
        display(sprintf('WARNING: Preprocessing step %s has not been defined, please check the spelling!',proclist{s}));
        display('This step will be turned off for current results');
        substr{s,1} = '[0]';
    else
        bleft       = ileft (ileft >iStep);   bleft= bleft(1);
        bright      = iright(iright>iStep);  bright=bright(1);
        substr{s,1} = pipelinelist(bleft:bright);
    end
end

% PIPE-1:Motion --------------------------------------------
pipeset_old=[];
pipeset_new=[]; K=1;
%
if( ~isempty( strfind(substr{K},'0') ) )  pipeset_new = [pipeset_new; 0];   end
if( ~isempty( strfind(substr{K},'1') ) )  pipeset_new = [pipeset_new; 1];   end

% PIPE-2:Censor --------------------------------------------
pipeset_old = pipeset_new;
pipeset_new = []; K=2;

%
%  LMP (2021-02-09) --- set alternate censor type if it is specified as 3rd arg ----
%  LMP (2021-02-25) --- make censorType last parameter in option list
%
%  'motion'       : outlier in motion parameter estimates (MPEs)
%  'volume'       : outlier in full-volume fMRI data
%  'volume+motion': outlier in both full-volume fMRI and MPEs
%  'slice'        : outlier in individual fMRI axial slices
%  'slice+motion' : outlier in both fMRI slice and MPEs

censortypelist = {'motion'; 'volume'; 'volume+motion'; 'slice';'slice+motion'};
censorType = 'volume+motion';
if((size(strfind(substr{K},','),2) > 0 ) )
    args = split(substr{K},',');
    lastindx = size(args,1);
    ct = split(args{lastindx},"'"); %check for quoted string
    if (size(ct,1) == 3)
        censorType = ct{2};
        if ~any(strcmp(censortypelist, censorType))
            censorType = 'volume+motion';
        end
    end
end

%
if( ~isempty( strfind(substr{K},'0') ) )
    tmpset     =[ pipeset_old, 0*ones(size(pipeset_old,1),1) ];
    pipeset_new=[pipeset_new; tmpset];
end
if( ~isempty( strfind(substr{K},'1') ) )
    tmpset     =[ pipeset_old, 1*ones(size(pipeset_old,1),1) ];
    pipeset_new=[pipeset_new; tmpset];
end
if( ~isempty( strfind(substr{K},'2') ) )
    tmpset     =[ pipeset_old, 2*ones(size(pipeset_old,1),1) ];
    pipeset_new=[pipeset_new; tmpset];
end
if( ~isempty( strfind(substr{K},'3') ) )
    tmpset     =[ pipeset_old, 3*ones(size(pipeset_old,1),1) ];
    pipeset_new=[pipeset_new; tmpset];
end

% PIPE-3:Retroicor --------------------------------------------
pipeset_old = pipeset_new;
pipeset_new = []; K=3;
%
if( ~isempty( strfind(substr{K},'0') ) )
    tmpset     =[ pipeset_old, 0*ones(size(pipeset_old,1),1) ];
    pipeset_new=[pipeset_new; tmpset];
end
if( ~isempty( strfind(substr{K},'1') ) )
    tmpset     =[ pipeset_old, 1*ones(size(pipeset_old,1),1) ];
    pipeset_new=[pipeset_new; tmpset];
end

% PIPE-4:Timecor --------------------------------------------
pipeset_old = pipeset_new;
pipeset_new = []; K=4;
%
if( ~isempty( strfind(substr{K},'0') ) )
    tmpset     =[ pipeset_old, 0*ones(size(pipeset_old,1),1) ];
    pipeset_new=[pipeset_new; tmpset];
end
if( ~isempty( strfind(substr{K},'1') ) )
    tmpset     =[ pipeset_old, 1*ones(size(pipeset_old,1),1) ];
    pipeset_new=[pipeset_new; tmpset];
end

% PIPE-5:Smooth --------------------------------------------
pipeset_old = pipeset_new;
pipeset_new = []; K=5;
fulidx  = [1  strfind(substr{K},',')  length(substr{K})];
numscal = length(fulidx)-1;
for(i=1:numscal)
    scaltok = substr{K}( fulidx(i)+1:fulidx(i+1)-1 );
    tmpset  =[ pipeset_old, str2num(scaltok)*ones(size(pipeset_old,1),1) ];
    pipeset_new=[pipeset_new; tmpset];
end

pipeset_half = pipeset_new; % everything that was already done
Nhalf = size( pipeset_half, 1 );

% ===========================================================

% PIPE-6:Detrend --------------------------------------------
detSet = [];
K=6;
if( ~isempty(strfind(substr{K},'-1')) || ~isempty(strfind(substr{K},'A')) )
    if( ~isempty(strfind(substr{K},',')) ) error('automatic detrending detected - cannot include other options'); end
    detSet = -1;
else
    fulidx  = [1  strfind(substr{K},',')  length(substr{K})];
    numscal = length(fulidx)-1;
    for(i=1:numscal)
        scaltok = substr{K}( fulidx(i)+1:fulidx(i+1)-1 );
        detSet = [detSet str2num(scaltok)];
    end
end

% PIPE-7:Motreg --------------------------------------------
mprSet = [];
K=7;
%
if( ~isempty( strfind(substr{K},'0') ) )
    mprSet     =[ mprSet, 0];
end
if( ~isempty( strfind(substr{K},'1') ) )
    mprSet     =[ mprSet, 1];
end

% PIPE-8:Taskreg --------------------------------------------
tskSet = [];
K=8;
%
if( ~isempty( strfind(substr{K},'0') ) )
    tskSet     =[ tskSet, 0];
end
if( ~isempty( strfind(substr{K},'1') ) )
    tskSet     =[ tskSet, 1];
end

% PIPE-9:phycaa --------------------------------------------
phySet = [];
K=9;
%
if( ~isempty( strfind(substr{K},'0') ) )
    phySet     =[ phySet, 0];
end
if( ~isempty( strfind(substr{K},'1') ) )
    phySet     =[ phySet, 1];
end

% PIPE-10+11:globalsig --------------------------------------------
gsSet = [];
K=10;
%
if( ~isempty( strfind(substr{K},'0') ) ) 
    if( ~isempty( strfind(substr{K+1},'0') ) )
        gsSet     =[0];
    end
    if( ~isempty( strfind(substr{K+1},'1') ) )
        gsSet     =[gsSet 2];
    end
end
if( ~isempty( strfind(substr{K},'1') ) )
    if( ~isempty( strfind(substr{K+1},'0') ) )
        gsSet     =[gsSet 1];
    end
    if( ~isempty( strfind(substr{K+1},'1') ) )
        gsSet     =[gsSet 3];
    end
end

% PIPE-12:lowpass filter --------------------------------------------
lpSet = [];
K=12;
%
if( ~isempty( strfind(substr{K},'0') ) )
    lpSet     =[ lpSet, 0];
end
if( ~isempty( strfind(substr{K},'1') ) )
    lpSet     =[ lpSet, 1];
end




Nfull = Nhalf * length(detSet) * length(mprSet) * length(tskSet) * length(gsSet) * length(lpSet);
