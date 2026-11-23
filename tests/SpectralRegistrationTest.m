classdef SpectralRegistrationTest < GannetTestCase
    % Spectral registration (Near et al., 2015) should recover known
    % frequency and phase offsets applied to synthetic transients.

    properties (Constant)
        NPoints     = 2048
        SW          = 2000      % Hz
        LarmorFreq  = 127.8     % MHz (3 T)
        NTrans      = 32
        MaxFreqHz   = 5         % largest applied freq. offset (Hz)
        MaxPhaseDeg = 20        % largest applied phase offset (deg)
    end

    methods (Test)
        function recoversFrequencyAndPhaseOffsets(testCase)
            rng(42, 'twister'); % for reproducibility
            [MRS_struct, df, dphi, fidTrue] = testCase.makeSyntheticData();

            evalc('[~, MRS_struct] = SpectralRegistration(MRS_struct, 0);');

            fEst   = MRS_struct.out.SpecReg.freq{1}(:);
            phiEst = MRS_struct.out.SpecReg.phase{1}(:);

            % The estimated corrections should cancel the applied
            % offsets, up to a common offset (the reference is the
            % median transient)
            fResid   = df + fEst;
            phiResid = dphi + phiEst;
            testCase.verifyLessThan(std(fResid), 0.05, ...
                ['Residual frequency error after registration is ' ...
                 'too large (Hz).']);
            testCase.verifyLessThan(std(phiResid), 1, ...
                ['Residual phase error after registration is ' ...
                 'too large (degrees).']);

            % Aligned transients should agree far better than unaligned
            % ones
            t = (0:testCase.NPoints-1).' / testCase.SW;
            ref = fidTrue .* exp(1i*2*pi*mean(fResid)*t) ...
                          .* exp(1i*pi/180*mean(phiResid));
            errBefore = mean(vecnorm(MRS_struct.fids.data - ref));
            errAfter  = mean(vecnorm(MRS_struct.fids.data_align - ref));
            testCase.verifyLessThan(errAfter, 0.1 * errBefore, ...
                'Alignment did not reduce the spread of transients by 10x.');
        end
    end

    methods (Access = private)
        function [MRS_struct, df, dphi, fid] = makeSyntheticData(testCase)
            N  = testCase.NPoints;
            sw = testCase.SW;
            LF = testCase.LarmorFreq;
            t  = (0:N-1).' / sw;

            % NAA, Cr, Cho, and Cr CH2 singlets (ppm, amplitude, T2* in s)
            peaks = [2.01 1.0 0.08
                     3.02 0.8 0.08
                     3.20 0.6 0.08
                     3.92 0.5 0.06];
            fid = zeros(N,1);
            for k = 1:size(peaks,1)
                f = (peaks(k,1) - 4.68) * LF; % Hz, relative to water
                fid = fid + peaks(k,2) * exp(1i*2*pi*f*t) ...
                                       .* exp(-t/peaks(k,3));
            end

            n    = testCase.NTrans;
            df   = (2*rand(n,1) - 1) * testCase.MaxFreqHz;
            dphi = (2*rand(n,1) - 1) * testCase.MaxPhaseDeg;
            df   = df - median(df);
            dphi = dphi - median(dphi);

            noiseSD = 0.002 * max(abs(fid));
            data = zeros(N, n);
            for k = 1:n
                data(:,k) = fid .* exp(1i*2*pi*df(k)*t) ...
                                .* exp(1i*pi/180*dphi(k)) ...
                            + noiseSD * complex(randn(N,1), randn(N,1));
            end

            MRS_struct = testCase.baseStruct();
            MRS_struct.p.npoints    = N;
            MRS_struct.p.sw         = sw;
            MRS_struct.p.LarmorFreq = LF;
            MRS_struct.p.ZeroFillTo = round(32768 / 2000 * sw);
            % Skip the in vivo Cr-based re-referencing step
            MRS_struct.p.phantom    = 1;
            MRS_struct.fids.data    = data;
        end
    end

end
