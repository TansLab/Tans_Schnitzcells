
function snew = renumberschnitzes(p,s);


N0 = length(s);
keepit = zeros(N0,1);
for i = 1:N0,
    keepit(i)=0;
    if ~isempty(s(i).frames),
        if length(s(i).frames)>0,
            keepit(i) = 1;
        end;
        if length(s(i).frames)==1,
            if (s(i).P<1),
                if s(i).D<1
                    if s(i).E<1,

                        % let's get rid of schnitzes that are totally disconnected from everything
                        % and have only a single frame
                        keepit(i) = 0;
                        disp(['found an isolated schnitz #',num2str(i),'--snipping it...']);

                    end;
                end;
            end;
        end;

    end;
end;




j = 1;
for i = 1:N0,
    if keepit(i),
        snew(j) = s(i);
        newnum(i) = j;
        j = j + 1;
    end;
end;
f = find(keepit==0);
disp('deleting schnitzes: ');
f'
% newnum,

for i = 1:length(snew),
    if snew(i).D>0,
        snew(i).D = newnum(snew(i).D);
    end;
    if snew(i).E>0,
        snew(i).E = newnum(snew(i).E);
    end;
    if snew(i).P>0,
        snew(i).P = newnum(snew(i).P);
    end;
end;

% now fix guys who have an E but no D...
for i = 1:length(snew),
    if (snew(i).E>0) & ~(snew(i).D>0),
        disp(['daughters swapped...fixing, new schnitz #',num2str(i)]);
        snew(i).D = snew(i).E;
        snew(i).E = 0;
    end;
end;

% finally, find everyone who lacks cenx/ceny definitions and fill those in:
for i = 1:length(snew),
    if length(s(i).frames)>0,
        if length(s(i).cenx)==0,
            disp(['adding cenx/ceny info to: ',num2str(i)]);
            snew(i) = addcenxceny(p,snew(i));
        end;
    end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newschnitz = addcenxceny(p,schnitz);

newschnitz = schnitz;
for i = 1:length(newschnitz.frames),
    [Lc] = loadseg(p,newschnitz.frames(i),'Lc');
    [x,y] = find(Lc==schnitz.cellno(i));
    newschnitz.cenx(i) = mean(y);
    newschnitz.ceny(i) = mean(x);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [varargout] = loadseg(p,fr,varargin);
filename = [p.segmentationDir,p.movieName,'seg',str3(fr-1),'.mat'];
if exist(filename),
    for i = 1:length(varargin),
        if ~isempty(who('-file',filename,char(varargin{i}))),
            x = load(filename,char(varargin{i}));
            varargout(i) = {x.(char(varargin{i}))};
        % added by nitzan 2005June24
        elseif strcmp(char(varargin{i}),'Lc') & p.trackUnCheckedFrames
            x = load(filename,'LNsub');
            varargout(i) = {x.LNsub};
        % end addition by nitzan 2005June24
        else
            disp(['cannot find ',char(varargin{i}),' in ',filename]);
            for k = 1:length(varargin),
                varargout(k)={0};
            end;
        end;
    end;
else
    disp([filename, ' does not exist']);
    for i = 1:length(varargin),
        varargout(i) = {0};
    end;
end;
