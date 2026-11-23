classdef DataReaderTest < GannetTestCase
    % Each vendor reader should load a Big GABA MEGA-PRESS dataset
    % (TE = 68 ms) with the expected dimensions and acquisition parameters.

    methods (Test)
        function readsGEPfile(testCase)
            f = testCase.dataFile('GE_data', 'G6_MP', 'S01', ...
                                  'S01_GABA_68.7');
            MRS_struct = GERead(testCase.baseStruct(), f);

            testCase.verifyMEGAPRESS(MRS_struct, 2048, 2000);
            testCase.verifySize(MRS_struct.fids.data, [2048 160]);
            testCase.verifyEqual(MRS_struct.p.Navg(1), 320);
            testCase.verifyEqual(MRS_struct.p.LarmorFreq(1), ...
                                 127.721838, 'AbsTol', 1e-3);
        end

        function readsPhilipsSDAT(testCase)
            f  = testCase.dataFile('Philips_data', 'P3_MP', 'S01', ...
                                   'S01_GABA_68_act.sdat');
            fw = testCase.dataFile('Philips_data', 'P3_MP', 'S01', ...
                                   'S01_GABA_68_ref.sdat');
            MRS_struct = PhilipsRead(testCase.baseStruct(), f, fw);

            testCase.verifyMEGAPRESS(MRS_struct, 2048, 2000);
            testCase.verifySize(MRS_struct.fids.data, [2048 320]);
            testCase.verifyEqual(MRS_struct.p.Navg(1), 320);
            testCase.verifyEqual(MRS_struct.p.LarmorFreq(1), ...
                                 127.75092, 'AbsTol', 1e-3);
        end

        function readsSiemensTWIX(testCase)
            f  = testCase.dataFile('Siemens_data', 'S5_MP', 'S01', ...
                                   'S01_GABA_68.dat');
            fw = testCase.dataFile('Siemens_data', 'S5_MP', 'S01', ...
                                   'S01_GABA_68_H2O.dat');
            % SiemensTWIXRead prints progress; capture it to keep the
            % test log tidy
            MRS_struct = testCase.quietly(@() ...
                SiemensTWIXRead(testCase.baseStruct(), f, fw));

            testCase.verifyMEGAPRESS(MRS_struct, 4096, 4000);
            testCase.verifySize(MRS_struct.fids.data, [4096 320]);
            testCase.verifyEqual(MRS_struct.p.Navg(1), 320);
            testCase.verifyEqual(MRS_struct.p.LarmorFreq(1), ...
                                 123.240321, 'AbsTol', 1e-3);
        end
    end

    methods (Access = private)
        function verifyMEGAPRESS(testCase, MRS_struct, npoints, sw)
            p = MRS_struct.p;
            testCase.verifyEqual(p.npoints(1), npoints);
            testCase.verifyEqual(p.sw(1), sw, 'AbsTol', 1e-6);
            testCase.verifyEqual(p.TE(1), 68, 'AbsTol', 1e-6);
            testCase.verifyEqual(p.TR(1), 2000, 'AbsTol', 1e-6);
            testCase.verifyEqual(p.voxdim(1,:), [30 30 30], ...
                                 'AbsTol', 1e-6);

            data = MRS_struct.fids.data;
            testCase.verifyEqual(size(data,1), npoints);
            testCase.verifyEqual(mod(size(data,2), 2), 0, ...
                ['MEGA-PRESS data should contain an even number of ' ...
                 'transients (ON/OFF pairs).']);
            testCase.verifyTrue(all(isfinite(data(:))), ...
                'Metabolite data contain non-finite values.');
            testCase.verifyGreaterThan(max(abs(data(:))), 0, ...
                'Metabolite data are all zero.');

            testCase.assertTrue(isfield(MRS_struct.fids, 'data_water'), ...
                'Water reference was not loaded.');
            water = MRS_struct.fids.data_water;
            testCase.verifyEqual(size(water,1), npoints);
            testCase.verifyTrue(all(isfinite(water(:))), ...
                'Water data contain non-finite values.');
            % Unsuppressed water should be much larger than the
            % water-suppressed signal
            testCase.verifyGreaterThan(max(abs(water(:,1))), ...
                                       5 * max(abs(mean(data,2))));
        end
    end

end
