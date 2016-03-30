
function [exptimestr, gainstr,exptime,cube,datenumber] = imsettings(pname, color);

%%
exptimestr = 'empty';
gainstr = 'empty';
exptime = 'empty';
cube = 'empty';
datenumber = 'empty';

%% 
if isempty(findstr('.tif', pname)),
    imx = [pname,'.tif'];
end;
if nargin > 1,
    pos = findstr('-p-', pname);
    pname(pos+1) = color;
end;

%%
if exist(pname)==2,
    iminfo = imfinfo(pname);
    
    if isfield(iminfo,'ImageDescription'),
        descrip = iminfo.ImageDescription;
    else
        descrip = [];
    end;
    
    
    % Read Software
    if isfield(iminfo,'Software'),  %(metamorph/metaseries: orig images)
        software = iminfo.Software;
    elseif ~isempty(descrip) %DE 2013/11/05 (for micromanager: orig&crop, for metamorph/metaseries: crop images)
        descrip = iminfo.ImageDescription;
        ps1=(strfind(descrip,'Software: '))+length('Software: ');
        software = descrip(ps1:end);
    else
        software = [];
    end;
    
    
    % Read Acquisition Date
    if isfield(iminfo,'DateTime'),  %(metamorph/metaseries: orig images)
        datetime = iminfo.DateTime;
    elseif ~isempty(strfind(descrip,'DateTime')) & strmatch('MetaMorph',software)   % (metamorph 7.1: crop images)
        ps2=(strfind(descrip,'DateTime: '))+length('DateTime: ');
        datetime = descrip(ps2:ps2+18); %to be tested! NW 2014-11
        
    elseif ~isempty(strfind(descrip,'DateTime')) & strmatch('MetaSeries',software)   % (metamorph 7.8 (=metaseries): crop images)
        ps2=(strfind(descrip,'DateTime: '))+length('DateTime: ');
        datetime = descrip(ps2:ps2+20);
    
    elseif ~isempty(strfind(descrip,'DateTime')) %DE 2013/11/05  (micromanager and hopefully no other cases)
         ps2=(strfind(descrip,'DateTime: '))+length('DateTime: ');
         datetime=descrip(ps2:ps2+18);
         
         year=str2num(datetime(1:4));
         month=str2num(datetime(6:7));
         day=str2num(datetime(9:10));
         hour=str2num(datetime(12:13));
         minute=str2num(datetime(15:16));
         second=str2num(datetime(18:19));
         
         datenumber = datenum(year,month,day,hour,minute,second);
         
         %DE 2014-04-11: exposure time is read from descrip too:
         exptimepos = findstr('Exposure: ',descrip) + length('Exposure: ');
         exptime = sscanf(descrip(exptimepos:end),'%f');
         exptimestr = num2str(exptime);
         
         %     elseif isfield(iminfo,'FileModDate'),     %%%ADDED SJT
         %         timestr = iminfo.FileModDate(13:20);
         %         hour = str2num(timestr(1:2));
         %         minute = str2num(timestr(4:5));
         %         second = str2num(timestr(7:8));
         %         year = 2000;
         %         month = 10;
         %         day = 10;
         %         datenumber = datenum(year,month,day,hour,minute,second);
    else
        datetime = [];
        disp(['Current ''descrip'' value:' descrip '.']);
        warning('Could not find date (yet).');
            % MW note: I think this function could be organized a bit
            % better; date is determined either in these if-statements
            % (micromanager) or later in code (metamorph).
    end;
    % ST end
else
    exptimestr = '';
    gainstr = '';
    exptime = 0;
    cube=-1;
    return;
end;

%%
if isempty(descrip),
    exptimestr = '';
    gainstr = '';
    exptime = 0;
    cube=-1;
    return;
end;

if strmatch('Image-Pro',software)   % unused software in Amolf
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

elseif strmatch('MetaMorph',software)     % old Metamorph version 7.1. Used at Amolf until 2014-10 
    exptimepos = findstr('Exposure: ',descrip) + length('Exposure: ');
    exptime = sscanf(descrip(exptimepos:end),'%f');
    exptimestr = num2str(exptime);
    %     cubestrpos = findstr('Cube=',descrip) + length('Cube=');
    %     cube = sscanf(descrip(cubestrpos:end),'%d');

    datestr = datetime(1:10);   % example: 2013:12:14
    timestr = datetime(12:19);

    year = str2num(datestr(1:4));
    month = str2num(datestr(6:7));
    day = str2num(datestr(9:10));

    hour = str2num(timestr(1:2));
    minute = str2num(timestr(4:5));
    second = str2num(timestr(7:8));
    
elseif strmatch('MetaSeries',software)      % new Metamorph version 7.8. Used at Amolf from 2014-10 on
    exptimepos = findstr('Exposure: ',descrip) + length('Exposure: ');
    exptime = sscanf(descrip(exptimepos:end),'%f');
    exptimestr = num2str(exptime);
    
    datestr = datetime(1:8);    % example: 20141111   (different in old metamorph)
    timestr = datetime(10:end);  %
    
    year = str2num(datestr(1:4));
    month = str2num(datestr(5:6));
    day = str2num(datestr(7:8));

    hour = str2num(timestr(1:2));
    minute = str2num(timestr(4:5));
    second = round(str2num(timestr(7:end)));    % round to full seconds (to adjust to older version)

% elsif micromanager
% note that when software is micromanager, timestamp data is parsed
% earlier.
end

%%

datenumber = datenum(year,month,day,hour,minute,second);

% keyboard;
% if exptime < 1,
%     exptimestr = ['0',exptimestr];
% end;
exptimestr(exptimestr=='.')=[];

%%

% gainpos = findstr('Gain: Gain ',descrip) + length('Gain: Gain ');
% ORIGINAL LINE - MW
gainpos = strfind(descrip, 'Gain: Gain ') + length('Gain: Gain ');

% EDITED -MW
if isempty(gainpos)
    % Use high setting if not found.    
    disp('can''t find gain setting -- set to n/a'); %NW2014-11
    gainstr = 'n/a';
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
            disp('gain setting not recognized -- set to n/a'); %NW2014-11
            gainstr = 'set to n/a';
    end;
end;

%%
end
