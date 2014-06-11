% [schnitz cellid]=recalc_schnitz(P,D,E,G,varargout)
% recalculates cellid correspondence cell array
% and schnitz structure with fields P,D,E,frames,.. and anything else
% being in the format {frame}(id); 

%     [schnitz cellid]=recalc_schnitz(P,D,E,G,varargin)
% 
% recalculates the 'schnitz' structure array (with cell ID indexing) 
% and the 'cellid'  {frame}(label) -> cell number correspondence array 
% out of the variables P,D,E,G, and optional 
% arguments in the format   ...'field_name', field_data... which will  be added
% into the schnitz.'field_name' field of the structure  array.
%
%   All input data and  must be cell arrays with indexing:  data{frame}(label_id) ;
% P,D,E,G -are  correspondencies between label IDs 
% on consequtive frames mbatching same cells:
% P{nframe}(j) corresponds to the label from (nframe+1) frame 
% matching  the cell labeled j on the nframes-th frame, if it did not
% split.
% D,E -same for the case where there was a split into two other cells,
% G-logical index array indicating labels of 'ghost' cells which don't have
% a match in previous frames.
%
% The P,D,E,G arrays can obtained from raw 'softassign output' with
% the data_treat function.
% T.B. 04/05

% JCR Modified: added trackRange argument 

function [schnitz cellid]=recalc_schnitz(P,D,E,G,trackRange,varargin)

[cellid PP DD EE]=cellids(P,D,E,G);
schnitz=[];




schnitz=toschnitz(cellid,schnitz,'P',PP);
schnitz=toschnitz(cellid,schnitz,'E',EE);
schnitz=toschnitz(cellid,schnitz,'D',DD);


for i=1:size(cellid,   2)
for j=1:max(size(cellid{i}))    
    % JCR Modified
    % frames{i}(j)=i;
    % frames{i}(j)=trackRange(i)+1;  % JCR Hack- add 1 for schntizedit  % MW 2014/06/11 removal JCR hack, aka removal N+1 bug
    
    frames{i}(j)=trackRange(i);  % MW 2014/06/11 removal JCR hack, aka removal N+1 bug
end
end
    
schnitz=toschnitz(cellid,schnitz,'frame_nrs',frames);

for i=1:2: floor(length(varargin)/2)*2;
    schnitz=toschnitz(cellid,schnitz,varargin{i},varargin{i+1});
end


