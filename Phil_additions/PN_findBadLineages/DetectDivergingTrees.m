%5.3 time points per cell cyle on average
clear all
close all
p = DJK_initschnitz('pos5crop','2013-08-14','e.coli.AMOLF','rootDir','D:\ExperimentalDataTodo\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040],'fluor1','r','fluor2','g','fluor3','none');
[p,schnitzcells] = DJK_compileSchnitzImproved_3colors(p,'quickMode',1);

%%
munb = '11';
ls = length(schnitzcells);
fieldName = 'dR5_cycCor'; %'G6_mean';
fieldTime = 'dR5_time';  %'R_time';
smoothField = [fieldName '_smoothed'];
% DON't USE EXTENSIVE =1 YET -> FUNCTIONS ARE NOT READY YET.
%extensive=1;  % =1: prod. rate, volume, total fluor.    =0  (-intensive) :concentration, growthrate
extensive=0;
% extensive=0 is not so perfect for e.g. prod rate because the two daughers
% have half the prod rate of the parent (on average)


if extensive==0  % conc etc
    schnitzcells = SmoothTree(schnitzcells,fieldName,31); %smooth each schnitz, taking into account the parents and the children tree
else % rate etc
    schnitzcells = SmoothTree_ext(schnitzcells,fieldName,31);
end

% timeaxis = [schnitzcells(1).(fieldTime)(:)  DescendAverage(schnitzcells,fieldTime,1)];  (%works with one or more init
% schnitzes)
%population = [schnitzcells(1).(smoothField)(:) ; DescendAverage(schnitzcells,smoothField,1)]; %compute the overall population average trend
% I think this assumes that the colony starts of with one single schnitz
% (1) !!
% more general solution with more initial schnitzes:
alltime=[schnitzcells.(fieldTime)];
allfield=[schnitzcells.(fieldName)];
% remove NaN fields  % be careful with the order of steps!
allfield=allfield(~isnan(alltime));
alltime=alltime(~isnan(alltime));
alltime=alltime(~isnan(allfield));
allfield=allfield(~isnan(allfield));
% check the NaN removal
idx1=find(isnan(alltime)); idx2=find(isnan(allfield));
if ~isempty(idx1) | ~isempty(idx2)
    disp('error: Time or Field has NaN values. Remove them!')
end

timaxis=unique(alltime);
usedpoints=0;
population=[];
for i=1:length(timeaxis)
    idxfield=find(alltime>timeaxis(i)-0.2 & alltime<timeaxis(i)+0.2);  % some tolerance
    population=[population; mean(allfield(idxfield))];
    usedpoints=usedpoints+length(idxfield);
end

%%
%find diverging subtrees
distance = zeros(1,ls);
for i = 1:ls
    s = schnitzcells(i);
    parentsAverage = AscendField(schnitzcells,smoothField,i);
    kidsAverage = DescendAverage(schnitzcells,smoothField,i);
    
    %reconstitutes a lineage passing through a particular schnitz, by
    %concatenating the signal of the parent, of the considered schnitz and the average of the children tree
    line = [parentsAverage(:) ; s.(smoothField)(:) ; kidsAverage(:)]'; 
    
    %compute the distance of this reconstituted lineage to the whole
    %population average
    distance(i) = sqrt(sum((line(:)-population(:)).^2));
end

[m im] = max(distance(:));

%output the list of schnitzes order by the divergence of their lineage
[Y,I] = sort(distance);
[Y',I']  % output