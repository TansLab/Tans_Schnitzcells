function []=NW_Backup_SegTrack(p)
% Creates a backup of 
% 1) segmentation files (in segmentation folder), e.g. pos5cropseg003.m
% 2) tracking files (in data folder), e.g. pos5crop-djk-output-003-to-004.txt
%
% A backup of the lineage file is not created but it can easily be
% recovered by running DJK_tracker_djk with the backupped tracking files.
%
%
% ********** IMPORTANT ***************
% - The function should be run directly after correction of
%   tracking (DJK_tracker_djk.m).
% - The timestamp (last modification time) of segmentation files is not
%   altered during backup creation.
% - The timestep of tracking files is updated to the current date to ensure
%   that they are always newer than segmentation files.
% 
% This implies: 
% - By reperforming tracking, the DJK_tracker does not (accidently) overwrite
%   the tracking files because they are newer than the segmentation files.
% - If you wish them to be overwritten you need to save the segmentation
%   file again or set ('override',1) in DJK_tracker.
%
% - Since the function updates the tracking timestamp it can also still be
%   applied once the fluo-addition (-> update of segmentation timestamp, 
%   NW_initializeFluoData) has been performed.
% ************************************
%
% ************Reason for backup***************
% After fluorescence analysis segmentation files are
% saved with a newer date (because fluorescence data is added). This date is
% newer than the date of the tracking file. If now an error in segmentation
% or tracking is discovered and tracking is redone all the manually
% corrected errors within the tracking files will be overwritten by the
% automated tracking.
% To avoid this, copy the backuped segmentation and tracking-output files
% back into the standard folders and correct then the error.
% ********************************************
%
% Backup is only created if the files exist, otherwise the file type is
% skipped.
%
% A unique subfolder in /data/ is created to store the backup. The name of
% the folder contains the date and a number if more backups are created on
% the same day.
% e.g. '/AutomaticBackup_2013-11-28/'
%      '/AutomaticBackup_2013-11-28_2/'  (newer version of the same day)
%
%

% -------------------------------------------------------------------------
% Create Save Directory for backup
% -------------------------------------------------------------------------
mySubDir=['AutomaticBackup_' datestr(date,'yyyy-mm-dd')];
mySaveDir=[p.tracksDir mySubDir filesep];
% check if directory already exists
if exist(mySaveDir)~=7 % if it does not exist -> create it
  [status,msg,id] = mymkdir([mySaveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' mySaveDir ' : ' msg]);
    return;
  end
else % if it exists then create a new Backup folder
    diridx=2;
        while exist(mySaveDir)==7  % loop over folder names (_2, _3, etc) until a folder name does not exist yet
            if diridx==2 %first round -> only remove one charachter (filesep)
                mySaveDir=[mySaveDir(1:end-1) '_' num2str(diridx) filesep]; 
            else
                mySaveDir=[mySaveDir(1:end-3) '_' num2str(diridx) filesep]; 
            end
        diridx=diridx+1;
    end
    % create new directory
    [status,msg,id] = mymkdir([mySaveDir]);
    if status == 0
        disp(['Warning: unable to mkdir ' mySaveDir ' : ' msg]);
        return;
    end
end


% -------------------------------------------------------------------------
% Copy files to the backup folder
% -------------------------------------------------------------------------
%  ------------------Segmentation files -----------------------
mysegfiles=dir([p.segmentationDir p.movieName 'seg*.mat']);
if ~isempty(mysegfiles)
    segSourceDir=[p.segmentationDir p.movieName 'seg*.mat'];
    [status,msg,id]=copyfile(segSourceDir,mySaveDir);
    if status==0
        disp(['Warning: unable to copy segmentation files : ' msg]);
    return;
    end
    disp(['Segmentation files back-upped.']);
else
    disp('No segmentation files found. Will continue with tracking files.')
end


% ------------------Tracking Files----------------------------

% ------ don't update modification date of tracking files ------
%mytrackfiles=dir([p.tracksDir p.movieName '*output*.txt']);
%if ~isempty(mytrackfiles)
%    trackSourceDir=[p.tracksDir p.movieName '*output*.txt'];
%    [status,msg,id]=copyfile(trackSourceDir,mySaveDir);
%    if status==0
%        disp(['Warning: unable to copy tracking files : ' msg]);
%    return;
%    end
%else
%    disp('No tracking files found.')
%end
% -------------------

% ------ update modification date of tracking files to current date (newer than segfiles) ------
mytrackfiles=dir([p.tracksDir p.movieName '*output*.txt']);
if ~isempty(mytrackfiles)
    % copy files to new folder
    trackSourceDir=[p.tracksDir p.movieName '*output*.txt'];
    [status,msg,id]=copyfile(trackSourceDir,mySaveDir);
    if status==0
        disp(['Warning: unable to copy tracking files : ' msg]);
    return;
    end
    % open the copied files , 'touch' them and save again -> modification
    % date updated to today
    for i=1:length(mytrackfiles)
        currentfile=[mySaveDir mytrackfiles(i).name];
        fid = fopen(currentfile, 'r+');   %'r+' open for reading and writing (don't create new file)
        byte = fread(fid, 1); %reads first byte (character or number or blank) of a binary file
        fseek(fid, 0, 'bof'); %move to beginning of file
        fwrite(fid, byte); %write first byte back into the file
        fclose(fid); %close file again
    end
    disp(['Tracking files back-upped.']);
else
    disp('No tracking files found.')
end
% ------------------


disp('-------------------------------------------')
disp('Finished copying files to backup folder.')
disp(['(' mySaveDir ')']);
disp('-------------------------------------------')

