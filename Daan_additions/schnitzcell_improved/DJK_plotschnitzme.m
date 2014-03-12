
function [xout,yout,xavgout,yavgout] = DJK_plotschnitzme(schnitzcells, xfield, yfield, whichones, varargin)

% function DJK_plotschnitzme(schnitzcells, xfield, yfield, whichones, option1, option2,...)
%
% same as plotschnitzme, but with few changes:
%   1. removed style option, always same style, different color for each
%   schnitz
%
%
%
%
% EXAMPLES:
%  - DJK_plotschnitzme(schnitzcells,'FYsmins','MYs',[]);
%  - DJK_plotschnitzme(schnitzcells,'FYsmins','MYs',[1:63],'DERIV FYsmins');
%  - DJK_plotschnitzme(schnitzcells,'FYsmins','MYs',[],'AVG');
%
 
if isfield(schnitzcells,'approved')
    if nargin < 4,
        whichones = find([schnitzcells.approved]);
        if length(whichones)
            disp('Showing approved schnitzes (and their ancestors) only.');
        end
    end;
    if isempty(whichones),
        whichones = find([schnitzcells.approved]);
        if length(whichones)
            disp('Showing approved schnitzes (and their ancestors) only.');
        end
    end;
end
if isempty(whichones),
    whichones = find([schnitzcells.D]==0);
    if length(whichones)
        disp('Showing daughterless schnitzes (and their ancestors) only.');
    end
end;
if ~length(whichones)
    disp('No schnitzces selected. See "help plotschnitzme"');
end

xavgout = [];
yavgout = [];

derivon = 0; % if you use "deriv", then we plot the dyfield/dxfield versus xfield.
norm0on = 0;
normon = 0;
bfilton = 0;
normbothon = 0;
invxon = 0;
collapseon = 0;
indmaxon = 0;
finalmedon = 0;
avgon = 0;
subtracton = 0;
subtractval = 0;
gethandles = 0;
keepnans = 0;
if nargin>=5, %DJK 080805
    for i = 1:length(varargin),
        
        if ~isempty(findstr(upper(varargin{i}),'SUBTRACT')),
            s = upper(varargin{i});
            subtracton = 1;
            [T,R] = strtok(s);
            [T,R] = strtok(R);
            subtractval = str2num(T);
            disp(['subtracting ',num2str(subtractval)]);
        end;

        
        if ~isempty(findstr(upper(varargin{i}),'COLLAPSE')),
            collapseon = 1;
        end;
        
        if ~isempty(findstr(upper(varargin{i}),'FINALMED')),
            finalmedon = 1;
        end;
        
        if ~isempty(findstr(upper(varargin{i}),'INDMAX')),
            indmaxon = 1;
        end;
        
        if ~isempty(findstr(upper(varargin{i}),'BFILT')),
            bfilton = 1;
        end;
        
        if ~isempty(findstr(upper(varargin{i}),'INVX')),
            invxon = 1;
        end;
        
        if ~isempty(findstr(upper(varargin{i}),'NORM')),
            normon = 1;
        end;
        
        if ~isempty(findstr(upper(varargin{i}),'AVG')),
            avgon = 1;
        end;
        if ~isempty(findstr(upper(varargin{i}),'HANDLEBACK')),
            gethandles = 1;
        end;
        if ~isempty(findstr(upper(varargin{i}),'KEEPNANS')),
            keepnans = 1;
        end;
        if ~isempty(findstr(upper(varargin{i}),'DERIV')),
            derivon = 1;

            s = (varargin{i});

            [T,R] = strtok(s);
            [T,R] = strtok(R);
            derivfield = T;
            disp(['derivative with respect to ',num2str(derivfield)]);
            
        end;
        if ~isempty(findstr(upper(varargin{i}),'NZERO')),
            norm0on = 1;
            normon = 0;
            disp('norm 0 on');
        elseif ~isempty(findstr(upper(varargin{i}),'NBOTH')),
            normbothon = 1;
            normon = 0;
            disp('norm both on');
        end;
        
        if ~isempty(findstr(upper(varargin{i}),'')),
            bfilton = 1;
        end;
    end;
    
    if normon,
        disp('norm on!');
        eval(['allvals = [schnitzcells.',yfield,'];']);
        normfac = max(allvals);
    end;
end;

%disp(['map: ' num2str(length(whichones))]);
map = DJK_hsv(length(whichones)); % color map for lines
for i = 1:length(whichones),
    thisone = whichones(i);
    
    X = [];
    Y = [];
    D = [];
    done = 0;
    while ~done, 
        x = getfield(schnitzcells(thisone),xfield);
        y = getfield(schnitzcells(thisone),yfield);
        
        x = x(end:-1:1);
        y = y(end:-1:1);
        
        if subtracton,
            y = y - subtractval;
        end;
        
        if derivon,
            d = getfield(schnitzcells(thisone),derivfield);
            d = d(end:-1:1);
        end;
            
        X = cat(2,X,x);
        Y = cat(2,Y,y);
        if derivon,
            D = cat(2,D,d);
        end;
        
        thisone = schnitzcells(thisone).P; %DJK 080606
        %done = (thisone <=0); %DJK 080606
        
        done = 1; %DJK 080606
    end;
    X = X(end:-1:1);
    Y = Y(end:-1:1);
    
    
    if derivon,
        D = D(end:-1:1);
        
        
        
        dY = diff(Y);
        dD = diff(D);
        dYdD = dY./dD;
        Y = dYdD;
        X2 = X(2:end);
        X1 = X(1:end-1);
        X12 = 0.5*(X1+X2);
        X = X12;
    end;
    
    if normon
        Y = Y/normfac;
    elseif norm0on | normbothon,
        Y = Y/Y(1);
    elseif collapseon,
        Y = Y/Y(end);
    elseif indmaxon,
        Y = Y/max(Y);
    end;
    
    if normbothon,
        X = X/X(1);
    end;
    
    fnotnan = find(~isnan(X) & ~isnan(Y));
    
    if ~keepnans
        X2 = X(fnotnan);
        Y2 = Y(fnotnan);
    else
        X2 = X;
        Y2 = Y;
    end
    
    if bfilton & (length(Y2)>5),
        Y2 = bfilt(Y2);
    end;
    
    
    if invxon,
        X2 = 1./X2;
    end;
    
    if ~avgon & length(X2),
        hout(i)=plot(X2,Y2,'Color',map(i,:),'LineWidth',3);  %DJK 080606
        hold on;
    end;
    yout{i} = Y2;
    xout{i} = X2;
    
end;


if avgon,
    
    allxout = [xout{:}];
    allyout = [yout{:}];
    avgX = [];
    avgY = [];
    
    ux = unique(allxout);
    for i = 1:length(ux),
        uu = ux(i);
        f = find(allxout==uu);
        avgX(i) = uu;
        avgY(i) = mean(allyout(f));
        stdY(i) = std(allyout(f));
    end;
    
    if normon,
        avgY = avgY / max(avgY);
    end;

    h = errorbar(avgX,avgY,stdY,'b-'); %DJK 080805
    xavgout = avgX;
    yavgout = avgY;
    
end;
if gethandles
    xout = hout;
end;


fs=findstr(yfield,'_');
for i=fs(end:-1:1);
    yfield(i+1:end+1)=yfield(i:end);yfield(i)='\';end
fs=findstr(xfield,'_');
for i=fs(end:-1:1);
    xfield(i+1:end+1)=xfield(i:end);xfield(i)='\';end
if derivon,
    ylabel(['derivative: \partial(',yfield,')/\partial(',derivfield,')']);
elseif normon,
    ylabel(['normalized ',yfield]);
else
    ylabel(yfield);
end;
xlabel(xfield);
hold off;

