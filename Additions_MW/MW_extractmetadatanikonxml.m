
%% Extract metadata test
% To be expanded later to import Nikon image in schnitzcells software.
% - MW 2016-06

% Give directory where your file is
ORIGINDIR = 'G:\EXPERIMENTAL_DATA_2016\a_incoming\2016-06-24_Petra_test_data\';

% File name
FILENAME = 'pos1-p-1-005.tif';
% (Remove the extension from the filename)
FILENAMEROOT = FILENAME(1:end-4);

% Extract the meta image information with imfinfo function
myiminfo=imfinfo([ORIGINDIR FILENAME]);
% Obtain the xml code (which also contains the timestamp)
thexmlcode = myiminfo.ImageDescription;

% In the xml code, find where the date+time are mentioned.
dateIndices = strfind(thexmlcode,'AcquisitionDate');
% Obtain that string
myDateRaw = thexmlcode(dateIndices(1)+numel('AcquisitionDate>'):dateIndices(2)-numel('</A'));

% And now polish that string a  little bit
numbersIndices = find(ismember(myDateRaw,'1234567890'));
myDateClean = myDateRaw(numbersIndices(1):numbersIndices(end));

%%

% Result in string
datetime=myDateClean

% We now the format, so we can actually extract all information separately
year=str2num(datetime(1:4));
month=str2num(datetime(6:7));
day=str2num(datetime(9:10));
hour=str2num(datetime(12:13));
minute=str2num(datetime(15:16));
second=str2num(datetime(18:19));

% Generate a timestamp (for further matlab calculation)
datenumber = datenum(year,month,day,hour,minute,second);

% Let user know what we found
disp(['For ' FILENAME ', date according to matlab is: ' num2str(year) '-' num2str(month) '-' num2str(day) ',' num2str(hour) ':' num2str(minute) ':' num2str(second) '.']);

