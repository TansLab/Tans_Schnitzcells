function [p, schnitzcells] = NDL_addToSchnitzes_skeletonLength(p)
framerange=[54];
%% % load schnitzcells
load([p.tracksDir p.movieName '-Schnitz.mat'])

% calculate all lengths
[p, allLengthsOfBacteriaInPixels, allLengthsOfBacteriaInMicrons] = NDL_lengthforfillamentedcells(p, framerange)
%% Assign calculated pixel-length to schnitzcells
schnitzcells(1).pixLength_filamented=[];
A=0;
for h=1:max(size(schnitzcells))
    for i=framerange
        for j=1:length(schnitzcells(h).frame_nrs)
            if (schnitzcells(h).frame_nrs(j)==i)==1
                A=schnitzcells(h).cellno(j);
                schnitzcells(h).pixLength_filamented(j)=allLengthsOfBacteriaInPixels{i}(A);
            end
        end
    end
end
%% Assign calculated length in micron to schnitzcells
schnitzcells(1).length_filamented=[];
B=0;
for k=1:max(size(schnitzcells))
    for l=framerange
        for m=1:length(schnitzcells(k).frame_nrs)
            if (schnitzcells(k).frame_nrs(m)==l)==1
                B=schnitzcells(k).cellno(m);
                schnitzcells(k).length_filamented(m)=allLengthsOfBacteriaInMicrons{l}(B);
            end
        end
    end
end

% use schnitzcells.frame_nrs and schnitzcells.cellno
% schnitzcells(s).length_fitNew(age) = quad(func_length, schnitzcells(s).fitNew_x_rot_left(age),schnitzcells(s).fitNew_x_rot_right(age)) * p.micronsPerPixel;
% schnitzcells(Schnitzcell_number).pixLength_filamented will look like 1x#frames double (1, frame_no)
% schnitzcells(1).pixLength_filamented=allLengthsOfBacteriaInPixels;
%% save schnitzcells
save([p.tracksDir p.movieName '-Schnitz.mat'])
end