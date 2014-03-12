function dispAndWrite(fid,text);
%
%

disp(text); 
fprintf(fid, '%s\n', text); 
