classdef EndToEndTest < GannetTestCase
    % Run GannetLoad and GannetFit on a Big GABA Philips MEGA-PRESS dataset
    % (S01, TE = 68 ms, with water reference) and check the outputs against
    % reference values from Gannet 3.5.4 with default settings
    % (metabolite ECC on; MATLAB R2026b, macOS).

    properties (Constant)
        % Tolerance on reference values, to allow for differences across
        % MATLAB releases and platforms
        RelTol = 0.10
        Ref = struct('GABA_ConcCr', 0.1185, 'GABA_ConcIU', 2.1277)
    end

    methods (Test)
        function loadAndFitPhilipsData(testCase)
            f  = testCase.dataFile('Philips_data', 'P3_MP', 'S01', ...
                                   'S01_GABA_68_act.sdat');
            fw = testCase.dataFile('Philips_data', 'P3_MP', 'S01', ...
                                   'S01_GABA_68_ref.sdat');

            % Gannet writes its output folders to the current folder
            testCase.applyFixture( ...
                matlab.unittest.fixtures.WorkingFolderFixture);

            MRS_struct = testCase.quietly(@() GannetLoad({f}, {fw}));
            testCase.assertTrue(isfield(MRS_struct, 'spec'), ...
                'GannetLoad did not produce spectra.');
            testCase.verifyTrue(isfile(fullfile('GannetLoad_output', ...
                'S01_GABA_68_act_vox1_load.pdf')), ...
                'GannetLoad did not write its output PDF.');

            MRS_struct = testCase.quietly(@() GannetFit(MRS_struct));
            testCase.verifyTrue(isfile(fullfile('GannetFit_output', ...
                'S01_GABA_68_act_vox1_fit.pdf')), ...
                'GannetFit did not write its output PDF.');

            out = MRS_struct.out.vox1;
            testCase.verifyEqual(out.GABA.ConcCr, ...
                testCase.Ref.GABA_ConcCr, 'RelTol', testCase.RelTol, ...
                'GABA+/Cr differs from the reference value.');
            testCase.verifyEqual(out.GABA.ConcIU, ...
                testCase.Ref.GABA_ConcIU, 'RelTol', testCase.RelTol, ...
                'GABA+ (i.u.) differs from the reference value.');

            % Data-quality sanity checks
            testCase.verifyLessThan(out.GABA.FitError, 10, ...
                'GABA+ fit error is too high (%).');
            testCase.verifyGreaterThan(out.GABA.SNR, 10, ...
                'GABA+ SNR is too low.');
            testCase.verifyLessThan(out.Cr.FWHM, 12, ...
                'Cr linewidth is too broad (Hz).');
            testCase.verifyLessThan(out.water.FWHM, 12, ...
                'Water linewidth is too broad (Hz).');
        end
    end

end
