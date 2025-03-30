function DeIdentifyDisclaimer
% Display a disclaimer when a user runs any of the de-identification
% functions

persistent lastDisclaimTime

if isempty(lastDisclaimTime) || (datetime('now') - lastDisclaimTime) > days(1)

    fprintf(['\nDISCLAIMER: The Gannet developers provide NO GUARANTEE ' ...
             'that the de-identification\nfunctions will remove all ' ...
             'protected health information (PHI) and personally\nidentifiable ' ...
             'information (PII) from data files.\n\n']);

    lastDisclaimTime = datetime('now');

end
