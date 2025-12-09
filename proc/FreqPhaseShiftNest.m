function output = FreqPhaseShiftNest(pars,input)
% Heavily copied from jn_freqPhaseShiftNest of Jamie Near, McGill

f   = pars(1); % frequency shift [Hz]
phi = pars(2); % phase shift [deg]

t = 0:input.dwelltime:(length(input.data)/2-1)*input.dwelltime;
input.data = reshape(input.data, [length(input.data)/2 2]);
fid = complex(input.data(:,1), input.data(:,2));

y = fid .* exp(1i*pi*(t'*f*2+phi/180));
output = [real(y); imag(y)];

end

