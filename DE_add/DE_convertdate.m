function [datenumber]=DE_convertdate(p,pname,color)

MetaDataPath= [p.movieDir(1:end-5),'\MetaData\']; 
%MetaDataPath= [pname(1:findstr(pname,'images')-1),'MetaData\'];

%get the index of the image:
%image_indx=str2num(pname(findstr(pname,'.tif')-4:findstr(pname,'.tif')-1));


image_indx1=strrep(p.movieName,'crop',[]);
image_indx2=strrep(image_indx1,'pos',[]);

%get the name of the metafile:
metafilename=[ image_indx1,'-',color,'-', image_indx2,'.txt'];
fullpath=[MetaDataPath,metafilename];
metaData=importdata(fullpath );

year=metaData{1}(end-3:end);
month=metaData{1}(5:7);
day=metaData{1}(9:10);
time=metaData{1}(12:19);

switch month
    case 'Nov'
    monthN=11;        
end

datetime=[year ':' num2str(monthN) ':' day ' ' time];

year=str2num(datetime(1:4));
month=str2num(datetime(6:7));
day=str2num(datetime(9:10));
hour=str2num(datetime(12:13));
minute=str2num(datetime(15:16));
second=str2num(datetime(18:19));


datenumber = datenum(year,month,day,hour,minute,second);
