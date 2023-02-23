GannetGUI is a graphical user interface designed to provide a front-end to the Gannet software package. To run Gannet with this GUI follow the steps below:

1. Select a BIDS-organized directory containing your Input Data.
2. Select a Derivatives folder that will contain Gannet outputs and associated files for the GUI. If you have Read and Write access to the Input Data Directory, the Derivatives Folder should be in there in accordance with BIDS.
3. If your Derivatives Folder has saved Workspaces, you will be prompted to load one. If not, you will work in a New Workspace.
4. Your Workspace shows which edited MRS scan sequences you want to process. Each row in the workspace is a different scan sequence found in your Input Data. For each row, select which Metabolite File sequence you want to process, along with its associated Water Reference File and Anatomical Image File.
5. If you do not have saved Configuration Files in your Derivatives Folder, make one for each type of Metabolite File you want to process.
6. Once each row has the correct Metabolite File, Water Reference File, Anatomical Image File, and Configuration File, you can select 'Save Current Workspace' to save the these rows and load them in for future uses of the GUI.
7. When ready, select 'Continue' to go to the Subject and Function Selector GUI. Note that you will not be able to select the 'Continue' button if any row is missing a Configuration File, or if the same scan sequence was selected multiple times.
8. On the Subject and Function Selector GUI, each tab shows a different scan sequence from the previous GUI. Each row shows a subject from the Input Data and whether that subject's Metabolite File, Reference File, and Anatomical Image File exist.
9. On each tab, choose at which step you would like to stop Gannet.
10. On each tab, select which subjects you would like to process for that particular scan sequence.
11. When all tabs are ready, select 'Run Gannet' to begin processing on all tabs.
12. Wait until Gannet finishes running, then close the GUI.

Quick descriptions of the buttons on the main Gannet GUI:
1. 'Reset Workspace' - Resets the screen to a New Workspace
2. 'Load Existing Workspace' - Allows you to load an existing Workspace found in the current Derivatives Folder
3. 'Save Current Workspace' - Saves the current Workspace to the Derivatives Folder
4. 'Create New Configuration' - Opens the Configuration File Editor, allowing you to make a new Configuration File and save it to the Derivatives Folder
5. 'View/Edit Existing Configuration' - Prompts which Configuration you would like to view or edit, then opens the Configuration File Editor and allows you to make changes to the file
6. 'Import Configuration' - Creates a copy of an existing Configuration File from somewhere other than your Derivatives Folder, allowing it to be used for the study being processed

For information on Gannet, please see 'README.md' in the 'Gannet-main' directory. For more information, see www.gabamrs.com/gannet.
