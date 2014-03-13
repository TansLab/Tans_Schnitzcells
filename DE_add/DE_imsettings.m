
function [exptimestr, gainstr,exptime,cube,datenumber] = DE_imsettings(p, pname, color)



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
    
    %Image info (timestamp and exposure) I get from Micromanager (p.micromanager==1)
    %This info is actually written in txt metadata file; I parse it DE_adjustinfo
    descrip = iminfo.ImageDescription;
    pos_software=(findstr(descrip,'Software: '))+length('Software: ');
    software = descrip(pos_software:end);
    
    pos_datetime=(findstr(descrip,'DateTime: '))+length('DateTime: ');
    datetime=descrip(pos_datetime:pos_datetime+18);
    
    year=str2num(datetime(1:4));
    month=str2num(datetime(6:7));
    day=str2num(datetime(9:10));
    hour=str2num(datetime(12:13));
    minute=str2num(datetime(15:16));
    second=str2num(datetime(18:19));
        
    datenumber = datenum(year,month,day,hour,minute,second);
        
    %Here is exposure
    pos_exposure=(findstr(descrip,'Exposure: '))+length('Exposure: ');
    exptimestr=descrip(pos_exposure:end);
    exptime=str2num(exptimestr);

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
% gain = sscanf(descrip(gainpos:end), '%f');
% switch(gain),
%     case 1,
%         gainstr = 'low';
%     case 2,
%         gainstr = 'med';
%     case 3,
%         gainstr = 'high';
%     otherwise,
%         disp('can''t find gain setting -- using high');
%         gainstr = 'high';
% end;

%DE: I impose gainstr. BTW, it is not used anyways.
gainstr = 'high';
