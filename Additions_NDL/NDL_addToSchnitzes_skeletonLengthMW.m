function [p, schnitzcells] = NDL_addToSchnitzes_skeletonLengthMW(p)

%% % load schnitzcells
load([p.tracksDir p.movieName '-Schnitz.mat'])
% Determine frame range
frameRange=unique([schnitzcells.frame_nrs]);
frameRange=sort(frameRange);

% calculate all lengths
[p, allLengthsOfBacteriaInPixels, allLengthsOfBacteriaInMicrons,allLengthsOfBacteriaInPixelsMW, allLengthsOfBacteriaInMicronsMW] = NDL_lengthforfillamentedcellsMW(p, frameRange)

%% Assign calculated pixel-length to schnitzcells
schnitzcells(1).pixLength_skeleton=[];
A=0;
for h=1:max(size(schnitzcells))
    for i=frameRange
        for j=1:length(schnitzcells(h).frame_nrs)
            if (schnitzcells(h).frame_nrs(j)==i)==1
                A=schnitzcells(h).cellno(j);
                schnitzcells(h).pixLength_skeleton(j)=allLengthsOfBacteriaInPixelsMW{i}(A);
                schnitzcells(h).length_skeleton(j)=allLengthsOfBacteriaInMicronsMW{i}(A);
            end
        end
    end
end
% %% Assign calculated length in micron to schnitzcells
% schnitzcells(1).length_skeleton=[];
% B=0;
% for k=1:max(size(schnitzcells))
%     for l=frameRange
%         for m=1:length(schnitzcells(k).frame_nrs)
%             if (schnitzcells(k).frame_nrs(m)==l)==1
%                 B=schnitzcells(k).cellno(m);
%                 schnitzcells(k).length_skeleton(m)=allLengthsOfBacteriaInMicrons{l}(B);
%             end
%         end
%     end
% end

% use schnitzcells.frame_nrs and schnitzcells.cellno
% schnitzcells(s).length_fitNew(age) = quad(func_length, schnitzcells(s).fitNew_x_rot_left(age),schnitzcells(s).fitNew_x_rot_right(age)) * p.micronsPerPixel;
% schnitzcells(Schnitzcell_number).pixLength_skeleton will look like 1x#frames double (1, frame_no)
% schnitzcells(1).pixLength_skeleton=allLengthsOfBacteriaInPixels;
%% save schnitzcells
save([p.tracksDir p.movieName '-Schnitz.mat'],'schnitzcells')
disp(['Schnitzcells file saved at ' p.tracksDir p.movieName '-Schnitz.mat'])
end