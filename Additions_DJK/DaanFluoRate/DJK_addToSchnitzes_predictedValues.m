% DJK_addToSchnitzes_predictedValues replaces NaN values with predicted
% values, based on how field depends on 3rd degree polynomial fit of predictionField
%
% OUTPUT
% 'schnitzcells'      schnitzcells
%
% REQUIRED ARGUMENTS:
% 'schnitzcells'      schnitzcells
% 'field'             name of field whose NaN need to be predicted
% 'predictionfield'   name of field that will be used to predict NaN
% 'newField'          name of new field
% 'minMax'            min and max that values can become
%
function [schnitzcells] = DJK_addToSchnitzes_predictedValues(schnitzcells, field, predictionfield, newField, minMax) 
PLOTTING = false;


%--------------------------------------------------------------------------
% INITIALIZE NEW FIELD
%--------------------------------------------------------------------------
schnitzcells(1).(newField) = [];
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% DETERMINE TREND
%--------------------------------------------------------------------------
x_data = [schnitzcells(:).(predictionfield)];
y_data =  [schnitzcells(:).(field)];
if length(x_data) ~= length(y_data)
  error(['length of fields ' predictionfield ' and ' field ' is not the same']);
end
idx = ~isnan(x_data); x_data = x_data(idx); y_data = y_data(idx);
idx = ~isnan(y_data); x_data = x_data(idx); y_data = y_data(idx);

% Fit with 3rd degree polynomial
fitCoef3 = DJK_polyfit(x_data,y_data,3);

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% LOOP OVER SCHNITZES
%--------------------------------------------------------------------------
for i = 1:length(schnitzcells)
  s = schnitzcells(i);
  
  %------------------------------------------------------------------------
  % LOOP OVER FIELD VALUES
  %------------------------------------------------------------------------
  for f = 1:length(s.(field))
    if isnan(s.(field)(f))
      s.(newField)(f) = polyval(fitCoef3,s.(predictionfield)(f));
      if s.(newField)(f) < minMax(1), s.(newField)(f) = minMax(1); end
      if s.(newField)(f) > minMax(2), s.(newField)(f) = minMax(2); end
    else
      s.(newField)(f) = s.(field)(f);
    end
  end
  %------------------------------------------------------------------------
  
  schnitzcells(i) = s;
end 
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% PLOTTING
%--------------------------------------------------------------------------
if PLOTTING
  figure;
  plot(x_data,y_data,'.'); hold on;
  range = [min(x_data):(max(x_data)-min(x_data))/50:max(x_data)];
  plot(range,polyval(fitCoef3,range),'r-');
end
%--------------------------------------------------------------------------
