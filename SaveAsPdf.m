clear all 
close all
%%
% A script to save all fig files as Pdf.
%%
folderStr = 'Output\New\';

Subfiles = dir( [folderStr,'\*.fig'] );

N_CASE = numel(Subfiles);

% set(0, 'DefaultAxesFontSize', 28);
     
for i_case = 1:N_CASE
    openfig([folderStr,'\',Subfiles(i_case).name]);
    
    fg = gcf;
    fg.Units = 'points';
    
    fg.PaperPositionMode = 'manual';
    fg.PaperOrientation  = 'landscape';
    fg.PaperUnits        = 'points';
    fg.PaperSize         = fg.Position(3:4)+1;
    fg.PaperPosition     = [1,1,fg.Position(3:4)];
    fg.PaperPositionMode = 'auto';

%     print([folderStr,'\',Subfiles(i_case).name(1:end-4),'.pdf'], '-dpdf')
    print([folderStr,'\',Subfiles(i_case).name(1:end-4),'.jpeg'], '-djpeg')
%     print([folderStr,'\',Subfiles(i_case).name(1:end-4),'.pdf'], '-dpdf', '-opengl')
    close gcf
end