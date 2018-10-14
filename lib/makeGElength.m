function g = makeGElength(g)
% if length not divisible by 4, pad with zeroes at end

if mod(length(g),4)
	[nrows,ncols] = size(g);
	if nrows > ncols
		g = [g; zeros(4-mod(length(g),4),ncols)];
	else
		g = [g  zeros(ncols,4-mod(length(g),4))];
	end
end

% EOF
