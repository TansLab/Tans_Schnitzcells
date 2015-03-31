% DJK_addToSchnitzes_cycleCor takes a schnitzcells data set, and corrects a
% parameter according to cell phase
%
% Procedure:
% Find the mean trend of any parameter against the cell phase (i.e. where
% the cell is in the cell cycle), and then correct that
% parameter for the cell phase (substraction).
%
% OUTPUT
% 'schnitzcells'      schnitzcells
% A field will be added to schnitzcells with the suffix _cyccor, which is
% the cell cycle corrected version of the given paramter.
%
% REQUIRED ARGUMENTS:
% 'schnitzcells'      schnitzcells
% 'field'             name of field to be corrected
% 'phaseField'        name of cell cycle phase field ('phase')
%
function [schnitzcells] = DJK_addToSchnitzes_cycleCor(schnitzcells, field, phaseField) 
PLOTTING = false;


%--------------------------------------------------------------------------
% INITIALIZE NEW FIELD
%--------------------------------------------------------------------------
newField                   = [field '_cycCor'];
schnitzcells(1).(newField) = [];
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% DETERMINE CORRECTION
%--------------------------------------------------------------------------
% field = 'dY5_sum_dt';
% phaseField = 'phase_atY';

data =  [schnitzcells(:).(field)];
phase = [schnitzcells(:).(phaseField)];
%if length(data) ~= length(phase)  % moved to further below, otherwise error
%                                  % for growth rate (NW 2012-09-11)
%  error(['length of fields ' field ' and ' phaseField ' is not the same']);
%end
idx = ~isnan(data);
data = data(idx);
phase = phase(idx);
idx = ~isnan(phase);
data = data(idx);
phase = phase(idx);
meanData = mean(data);
if length(data) ~= length(phase)  
  error(['length of fields ' field ' and ' phaseField ' is not the same']);
end


% % Fit with binning
% binWidth = 0.04;
% [N,bin] = histc(phase,[0:binWidth:1]);
% for i = 1:length(N)-1
%   meanOfBin(i) = mean(data(find(bin==i)));
% end
% 
% % Plot binning
% if PLOTTING
%   figure; plot(phase,data,'.'); hold on;
%   x = [binWidth:binWidth:1-binWidth];
%   for i=1:length(x)
%     x2(2*i-1) = x(i);
%     x2(2*i  ) = x(i);
%   end
%   x2 = [0 x2 1];
%   for i=1:length(meanOfBin)
%     meanOfBin2(2*i-1) = meanOfBin(i);
%     meanOfBin2(2*i  ) = meanOfBin(i);
%   end
%   plot(x2,meanOfBin2,'r-'); hold on;
%    %plot([0, x],meanOfBin,'.r','MarkerSize',15)
% end


% Fit with 3rd degree polynomial
fitCoef3 = DJK_polyfit(phase,data,3);
f = polyval(fitCoef3,[0:0.1:1]);
if PLOTTING
  figure;
  plot(phase,data,'.'); hold on;
  plot([0:0.1:1],f,'r-'); hold on;
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
  s = schnitzcells(i);
  
  %------------------------------------------------------------------------
  % LOOP OVER PHASE FIELD
  %------------------------------------------------------------------------
  if length(s.(phaseField)) > 0
    for f = 1:length(s.(phaseField))
      if isnan(s.(phaseField)(f))
        s.(newField)(f) = NaN;
      else
%         cycleCor = meanOfBin( ceil(s.(phaseField)(f)/binWidth) );

        % Determine extrapolated value at current cell phase
        cycleCor = polyval(fitCoef3,s.(phaseField)(f));
        % subtract that value, but keep same mean
        s.(newField)(f) = s.(field)(f) + meanData - cycleCor;
      end
    end
  end
  %------------------------------------------------------------------------
  
  schnitzcells(i) = s;
end 

if PLOTTING
  figure;
  data =  [schnitzcells(:).(newField)];
  phase = [schnitzcells(:).(phaseField)];
  if length(data) ~= length(phase)
    error(['length of fields ' newField ' and ' phaseField ' is not the same']);
  end
  idx = ~isnan(data);
  data = data(idx);
  phase = phase(idx);
  idx = ~isnan(phase);
  data = data(idx);
  phase = phase(idx);
  plot(phase,data,'g.'); hold on;
end



% function doubleData = DJK_doubleArray(data)
% for i=1:length(data)
%   doubleData(2*i-1) = data(i);
%   doubleData(2*i  ) = data(i);
% end


