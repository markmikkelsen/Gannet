function out = SpecReg(data, target, t, x)

f   = x(1);
phi = x(2);

data = reshape(data, [length(data)/2 2]);
fid  = complex(data(:,1), data(:,2));

fid = fid .* exp(1i * pi * (t' * f * 2 + phi/180));
fid = [real(fid); imag(fid)];

out = target - fid;

end
