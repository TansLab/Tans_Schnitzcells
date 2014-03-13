function im_description=DE_adjustiminfo(p, DphaseRange_name)

%% If the data is taken with micromanager, there has to be some
% readjustment of timestring and other iminfo.

% Inputs from the superfunction DJK_cropImages_3colors
%(p, DphaseRange(i).name)
fullpath=[p.movieDir 'MetaData\' DphaseRange_name(1:end-4) '.txt'];

% disp(['DE_adjustiminfo || fullpath = ''' fullpath '''']) % DEBUG MW

metaData=importdata(fullpath );

year=metaData{1}(end-3:end);
month=metaData{1}(5:7);
day=metaData{1}(9:10);
time=metaData{1}(12:19);

%exposuretime
exposure=metaData{2}(10:end);

% monthN=DE_month_convert(month); % ORIGINAL LINE -MW
% Convert to string which holds number of months, add prefix zero if req.
monthN= sprintf('%02d',  DE_month_convert(month)  ); % EDITED LINE -MW

% disp(['DE_adjustiminfo || month = ''' month '''']) % DEBUG MW
% disp(['DE_adjustiminfo || monthN = ' num2str(monthN) '']) % DEBUG MW

%DateTime=[year ':' str2(monthN) ':' day ' ' time]; % ORIGINAL LINE -MW
DateTime=[year ':' monthN ':' day ' ' time]; % EDIT MW

% disp(['DE_adjustiminfo || DateTime = ''' DateTime '''']) % DEBUG MW

%disp([year ' ' month ' ' day ' ' time] );
year=str2num(DateTime(1:4));
month=str2num(DateTime(6:7));
day=str2num(DateTime(9:10));
hour=str2num(DateTime(12:13));
minute=str2num(DateTime(15:16));
second=str2num(DateTime(18:19));

datenumber = datenum(year,month,day,hour,minute,second);


im_info = imfinfo([p.imageDir DphaseRange_name]);

%software info was already there:
Software=im_info.ImageDescription(1:12);
 

 % this image info will be added to crop
 im_description = ['DateTime: ' DateTime  char(10) 'Software: ' Software char(10) 'Exposure: ' exposure];
 
 
