function DICOMDeIdentify
% DICOMDeIdentify
%   1. Create a top-level directory and then create subdirectories for each
%      of your subjects and name them, e.g., 'S01', 'S02', 'S03', ... etc.
%   2. Copy each subject's DICOM (.dcm) files into their respective folder
%      that you created in step 1.
%   3. Run this script from within the top-level directory

%   Author: Mark Mikkelsen (Johns Hopkins University, 2018)

% Additional fields to remove
attrbs.StudyDate = '';
attrbs.SeriesDate = '';
attrbs.AcquisitionDate = '';
attrbs.ContentDate = '';
attrbs.StudyTime = '';
attrbs.SeriesTime = '';
attrbs.AcquisitionTime = '';
attrbs.ContentTime = '';
attrbs.ClinicalTrialProtocolName = '';
attrbs.PerformedStationName = '';
attrbs.PerformedLocation = '';
attrbs.PerformedProcedureStepDescription = '';
attrbs.PerformedProcedureStepStartDate = '';
attrbs.PerformedProcedureStepStartTime = '';
attrbs.PerformedProcedureStepID = '';
attrbs.PerformedProcedureStepDescription = '';

homeDir = pwd;
A = dir;
A = A(~ismember({A.name},{'.','..'}));
reverseStr = '';
for ii = 1:length(A)
    if isdir(A(ii).name) %#ok<ISDIR>
        if ~exist([A(ii).name '_anon'],'dir')
            mkdir([A(ii).name '_anon']);
        end
        msg = sprintf('\n\nAnonymizing: %s',A(ii).name);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
        cd(A(ii).name);
        B = dir;
        B = B(~ismember({B.name}, {'.','..','.DS_Store'}));
        fprintf('\n\n');
        for jj = 1:length(B)
            if exist(B(jj).name,'file') == 2
                dicomanon(B(jj).name, ['../' A(ii).name '_anon/' B(jj).name], 'update', attrbs);
                msg2 = sprintf([B(jj).name ' ---> '  A(ii).name '_anon/' B(jj).name '\n']);
                fprintf(msg2);
            end
        end
    end
    cd(homeDir);
end
cd(homeDir);
fprintf('\nDone!\n');



