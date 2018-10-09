function out = tryread(fn, arg)
try
	out = fn(arg);
catch ME
	error('Failed to read %s\n', arg);
end
