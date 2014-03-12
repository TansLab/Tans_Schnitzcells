% DJK_scaleRange scales a range of the input data to another range
% Values that end up outside the target range are set to the limit
% 
% OUTPUT:
% 'output'        output data
%
% REQUIRED ARGUMENTS:
% 'input'         input data
% 'source_range'  source range, eg [0 100]
% 'target_range'  source range, eg [0 100]
%
% formula, according to http://www.vias.org/tmdatanaleng/cc_scaling.html
 
function output = DJK_scaleRange(input, source_range, target_range);

% correct for error, when source range is not a range but a point
if source_range(1) == source_range(2)
  source_range = source_range + [-1 1];
end
  
scale_factor = (target_range(2) - target_range(1)) / (source_range(2) - source_range(1)) ;
shift_factor = ( target_range(1)*source_range(2) - target_range(2)*source_range(1) ) / (source_range(2) - source_range(1));

output = input * scale_factor + shift_factor;
output(find(output>target_range(2))) = target_range(2);
output(find(output<target_range(1))) = target_range(1);
