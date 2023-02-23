% This script should either be run while the pwd is Gannet_GUI_V0.3 (on
% startup of the GUI), or the folder containing Gannet-main, spm12, and
% Gannet_GUI_V0.3

currentFolder = pwd;
cd ../

if isfolder('Gannet-main')
    addpath('Gannet_GUI_V0.3');
    addpath('spm12');
    addpath('Gannet-main');
else
    cd(currentFolder)
end

clear currentFolder

GannetGUI
