function [schnitz cellid]=recalc_schnitz(P,D,E,G,trackRange,varargin)
% [schnitz cellid]=recalc_schnitz(P,D,E,G,varargout)
% recalculates cellid correspondence cell array
% and schnitz structure with fields P,D,E,frames,.. and anything else
% being in the format {frame}(id); 
%
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

% Hierarchy of calling:
% >> MW_makeSchnitzFileFromTracking 
%       >> DJK_data_treat
%       >> recalc_schnitz 
%           >> cellids

% JCR Modified: added trackRange argument 

%%

% cellid now contains a list of member cells per frames
[cellid PP DD EE]=cellids(P,D,E,G);
% Create empty schnitz
schnitz=[];

%%

% Add ancestry information to schnitz
schnitz=toschnitz(cellid,schnitz,'P',PP);
schnitz=toschnitz(cellid,schnitz,'E',EE);
schnitz=toschnitz(cellid,schnitz,'D',DD);

%% now make lookup table for which frame belongs to which schnitz, note 
% that the lookup table  labels the first frame in the range as "1"
for frameIdx=1:size(cellid,   2) % looping over frames
for j=1:max(size(cellid{frameIdx})) % looping over 
    % JCR Modified
    % frames{i}(j)=i;
    % frames{i}(j)=trackRange(i)+1;  % JCR Hack- add 1 for schntizedit  % MW 2014/06/11 removal JCR hack, aka removal N+1 bug
    
    frames{frameIdx}(j)=trackRange(frameIdx);  % MW 2014/06/11 removal JCR hack, aka removal N+1 bug
end
end

%%
% Compiling main structure - MW comment
schnitz=toschnitz(cellid,schnitz,'frame_nrs',frames);

%%
% Iterating over additional arguments that are added - MW comment
for i=1:2:floor(length(varargin)/2)*2;
    disp(['Adding [' varargin{i} '].']);
    schnitz=toschnitz(cellid,schnitz,varargin{i},varargin{i+1});
end



