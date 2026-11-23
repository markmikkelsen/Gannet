classdef (Abstract) GannetTestCase < matlab.unittest.TestCase
    % GANNETTESTCASE Shared setup for Gannet unit tests.
    %
    %   Adds Gannet to the MATLAB path, hides figures, and locates the test
    %   data. The test data root is read from the GANNET_TEST_DATA
    %   environment variable; if unset, tests/data is used. Tests that need
    %   data are skipped (marked Incomplete) when it is missing, unless the
    %   environment variable GANNET_TEST_DATA_REQUIRED is set to 1 (as in
    %   CI), in which case they fail.

    properties (SetAccess = private)
        GannetRoot
        DataRoot
    end

    methods (TestClassSetup)
        function addGannetToPath(testCase)
            testCase.GannetRoot = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                testCase.GannetRoot, 'IncludingSubfolders', true));
        end

        function hideFigures(testCase)
            origVisible = get(groot, 'DefaultFigureVisible');
            set(groot, 'DefaultFigureVisible', 'off');
            testCase.addTeardown(@() set(groot, ...
                'DefaultFigureVisible', origVisible));
            testCase.addTeardown(@() close('all', 'force'));
        end

        function locateData(testCase)
            root = getenv('GANNET_TEST_DATA');
            if isempty(root)
                root = fullfile(fileparts(mfilename('fullpath')), 'data');
            end
            testCase.DataRoot = root;
        end
    end

    methods (Access = protected)
        function f = dataFile(testCase, varargin)
            % Return the full path to a test-data file, or skip/fail if
            % it does not exist
            f = fullfile(testCase.DataRoot, varargin{:});
            msg = sprintf('Test data file not found: %s', f);
            if strcmp(getenv('GANNET_TEST_DATA_REQUIRED'), '1')
                testCase.fatalAssertTrue(isfile(f), msg);
            else
                testCase.assumeTrue(isfile(f), msg);
            end
        end

        function MRS_struct = baseStruct(~)
            % Minimal MRS_struct, as GannetLoad builds before calling a
            % reader
            MRS_struct = GannetPreInitialise(struct);
            MRS_struct.ii = 1;
        end

        function out = quietly(~, fn) %#ok<INUSD> used in evalc
            % Call fn and return its output, suppressing anything it
            % prints to the Command Window
            [~, out] = evalc('fn()');
        end
    end

end
