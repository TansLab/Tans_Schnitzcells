function MW_makeSchnitzFileFromTracking(p, opts)
% Function that is shared among all trackers.
% 
% See DJK_tracker_djk, MW_tracker, NW_tracker_centroid_vs_area.
% Updates schnitzcells .mat file with lineage based on tracking files.

disp('Updating schnitzcells lineage file..');

%%
[P, D, E, G] = DJK_data_treat(p);

%%
[schnitzcells, cellid] = recalc_schnitz(P,D,E,G,p.manualRange,'',opts); %
schnitzcells = renumberschnitzes(p,schnitzcells); % MW TODO, is this necessary step?

%%
disp(['saving schnitzcells lineage structure to ' p.lineageName]);
save(p.lineageName,'schnitzcells');

end