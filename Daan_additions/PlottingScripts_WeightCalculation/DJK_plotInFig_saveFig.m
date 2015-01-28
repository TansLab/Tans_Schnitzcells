% DJK_plotInFig_saveFig saves figure depending on settings:
% * figu.figure.save   = 0 no save
%                      = 1 save
%                      = 2 ask to save
%
% REQUIRED INPUTS:
%  * figu.fig{}
%
% OPTIONAL INPUTS
%  * figu.figure.filename   - 'unnamed figure'
%  * figu.figure.filetypes  - 'pdf'
%  * figu.figure.save       - 0
%  * figu.figure.show       - 1
%
%
% Code written by Daan Kiviet

function DJK_plotInFig_saveFig(figu);
% SET OPTIONAL INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([' * ' mfilename]);
figu = DJK_setOptionalInput(figu, ['figu.figure.filename']  , 'unnamed figure');
figu = DJK_setOptionalInput(figu, ['figu.figure.filetypes'] , 'pdf');
figu = DJK_setOptionalInput(figu, ['figu.figure.save']      , 0);
figu = DJK_setOptionalInput(figu, ['figu.figure.show']      , 1);

figu.figure = DJK_fillCellsToSize(figu.figure, length(figu.fig));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% SAVE FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(figu.fig)
  % DETERMINE WHAT TO DO 
  switch figu.figure.save{i}
    case 0
      fig_save = 0;
    case 1
      fig_save = 1;
    case 2
      saveFigInput = questdlg('Save Figure?','Save Figure?','Yes','No','Yes');
      pause(0.2);
      if (upper(saveFigInput(1))=='Y')
        fig_save = 1;
      else
        fig_save = 0;
      end;
    otherwise
      fig_save = 0;
  end

  % SAVE FIGURE 
  if fig_save
    for j = 1:length(figu.figure.filetypes{i})
      figureFilename = char(strcat(figu.figure.filename{i} , '.', figu.figure.filetypes{i}(j)));
      saveas(figu.fig{i},figureFilename);
      disp([' * Saved plot in ' figureFilename]);
    end
  end

  % CLOSE FIGURE 
  if ~figu.figure.show{i}
    close(figu.fig{i});
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
