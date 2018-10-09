function g = makeGElength(g)
% Input:
%   g    waveform array [nt ...]
% if length not divisible by 4, pad with zeroes at end

if mod(size(g,1),4)
	[nrows,ncols] = size(g);
	g = [g; zeros(4-mod(length(g),4),ncols)];
end

% EOF
