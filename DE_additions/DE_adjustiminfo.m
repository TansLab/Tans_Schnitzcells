function im_description=DE_adjustiminfo(p, DphaseRange_name)

%% If the data is taken with micromanager, there has to be some
% readjustment of timestring and other iminfo.

% Inputs from the superfunction DJK_cropImages_3colors
%(p, DphaseRange(i).name)
fullpath=[p.movieDir 'MetaData\' DphaseRange_name(1:end-4) '.txt'];

metaData=importdata(fullpath );

year=metaData{1}(end-3:end);
month=metaData{1}(5:7);
day=metaData{1}(9:10);
time=metaData{1}(12:19);

%exposuretime
exposure=metaData{2}(10:end);


monthN = sprintf('%02d',  DE_month_convert(month));

DateTime=[year ':' monthN ':' day ' ' time];


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
 
 
