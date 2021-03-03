function data = EddyCurrentCorrection(data, water_data)

% Based on Klose (1990), MRM, 14(1):26-30. The equation was taken from
% Jiru (2008), EJR, 67(2):202-217

K           = abs(data);
phase_data  = unwrap(angle(data));
phase_water = unwrap(angle(water_data));
phase_corr  = phase_data - repmat(phase_water, [1 size(data,2)]);
data        = complex(K .* exp(1i*phase_corr)); % make sure water data remain complex