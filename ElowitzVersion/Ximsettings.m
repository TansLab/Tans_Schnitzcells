
function [exptimestr, gainstr,exptime,cube,datenumber] = imsettings(pname, color);

exptimestr = 'empty';
gainstr = 'empty';
exptime = 'empty';
cube = 'empty';
datenumber = 'empty';

if isempty(findstr('.tif', pname)),
    imx = [pname,'.tif'];
end;
if nargin > 1,
    pos = findstr('-p-', pname);
    pname(pos+1) = color;
end;
if exist(pname)==2,
    iminfo = imfinfo(pname);
    if isfield(iminfo,'ImageDescription'),
        descrip = iminfo.ImageDescription;
    else
        descrip = [];
    end;
    % ST 4/10/05:
    if isfield(iminfo,'Software'),
        software = iminfo.Software;
    else
        software = [];
    end;
    if isfield(iminfo,'DateTime'),
        datetime = iminfo.DateTime;
    elseif isfield(iminfo,'FileModDate'),     %%%ADDED SJT
        timestr = iminfo.FileModDate(13:20);
        hour = str2num(timestr(1:2));
        minute = str2num(timestr(4:5));
        second = str2num(timestr(7:8));
        year = 2000;
        month = 10;
        day = 10;
        datenumber = datenum(year,month,day,hour,minute,second);
    else
        datetime = [];
    end;
    % ST end
else
    exptimestr = '';
    gainstr = '';
    exptime = 0;
    cube=-1;
    return;
end;
if isempty(descrip),
    exptimestr = '';
    gainstr = '';
    exptime = 0;
    cube=-1;
    return;
end;
if strmatch('Image-Pro',software)
    exptimepos = findstr('Exptime=',descrip) + length('Exptime=');
    exptime = sscanf(descrip(exptimepos:end),'%f');
    exptimestr = num2str(exptime);
    cubestrpos = findstr('Cube=',descrip) + length('Cube=');
    cube = sscanf(descrip(cubestrpos:end),'%d');

    timesetupstr = 'Acquired at ';
    timestrpos = findstr(timesetupstr,descrip)+length(timesetupstr);
    datetimestr = sscanf(descrip(timestrpos:exptimepos-1),'%s');
    datestr = datetimestr(1:10);
    timestr = datetimestr(11:18);

    year = str2num(datestr(1:4));
    month = str2num(datestr(6:7));
    day = str2num(datestr(9:10));

    hour = str2num(timestr(1:2));
    minute = str2num(timestr(4:5));
    second = str2num(timestr(7:8));

elseif strmatch('MetaMorph',software)     %%%%%%%%%%%%%%%ADDED SJT
    exptimepos = findstr('Exposure: ',descrip) + length('Exposure: ');
    exptime = sscanf(descrip(exptimepos:end),'%f');
    exptimestr = num2str(exptime);
    %     cubestrpos = findstr('Cube=',descrip) + length('Cube=');
    %     cube = sscanf(descrip(cubestrpos:end),'%d');

    datestr = datetime(1:10);
    timestr = datetime(12:19);

    year = str2num(datestr(1:4));
    month = str2num(datestr(6:7));
    day = str2num(datestr(9:10));

    hour = str2num(timestr(1:2));
    minute = str2num(timestr(4:5));
    second = str2num(timestr(7:8));
elseif strmatch('Exposure:',descrip)     %%%%%%%%%%%%%%%ADDED SJT not so nice..if can not add extra fields
    exptimepos = findstr('Exposure: ',descrip) + length('Exposure: ');
    exptime = sscanf(descrip(exptimepos:end),'%f');
    exptimestr = num2str(exptime);
    datetimepos = findstr('DateTime: ',descrip) + length('DateTime: ');
    datetimestr = descrip(datetimepos:end);
    %datetimestr = num2str(datetime);    
    datestr = datetimestr(1:10);
    timestr = datetimestr(12:19);

    year = str2num(datestr(1:4));
    month = str2num(datestr(6:7));
    day = str2num(datestr(9:10));

    hour = str2num(timestr(1:2));
    minute = str2num(timestr(4:5));
    second = str2num(timestr(7:8));
%    keyboard;
end

%end
datenumber = datenum(year,month,day,hour,minute,second);
% keyboard;
% if exptime < 1,
%     exptimestr = ['0',exptimestr];
% end;
exptimestr(exptimestr=='.')=[];

% gainpos = findstr('Gain: Gain ',descrip) + length('Gain: Gain ');
% ORIGINAL LINE - MW
gainpos = strfind('Gain: Gain ',descrip) + length('Gain: Gain ');

% EDITED -MW
if isempty(gainpos)
    % Use high setting if not found.    
    disp('can''t find gain setting -- using high');
    gainstr = 'high';
else
    % Mark beginning of gain value
    gain = sscanf(descrip(gainpos:end), '%f');

    % And convert numerical value to string value
    switch(gain),
        case 1,
            gainstr = 'low';
        case 2,
            gainstr = 'med';
        case 3,
            gainstr = 'high';        
        % Or otherwise not found
        otherwise,
            disp('gain setting not recognized -- using high');
            gainstr = 'high';
    end;
end;

