function G = makegirf
% Play gradients proposed in Vannesjo et al, MRM 69:583-593 (2013):
% 12 triangular gradients with same slew rate but different time-to-peak.
%
% $Id: makegirf.m,v 1.1 2017/09/19 16:20:57 jfnielse Exp $

mxg = 5;               % max gradient (Gauss/cm)
mxs = 18e3;            % max gradient slew rate (Gauss/cm/sec)

TimeToPeakUs = [50:10:160];  % microseconds

% construct triangular gradients (using 1 us dwell time)
for ip = numel(TimeToPeakUs):-1:1
	dwell = 1e-6;
	np = TimeToPeakUs(ip);       % number of 1 us samples to peak
	dg = mxs*dwell;              % change in gradient per 1 microsecond sample (Gauss/cm)
	gtmp = (0:dg:(dg*np))';
	gtmp = [zeros(4,1); gtmp];
	gtmp = [gtmp; flipud(gtmp)];

	if ip < numel(TimeToPeakUs)
		dn = N - numel(gtmp);
		gtmp = [zeros(dn/2,1); gtmp; zeros(dn/2,1)];
		G(:,ip) = gtmp;
	else
		G(:,ip) = gtmp;
		N = size(G,1);
	end
end

T = dwell:dwell:(dwell*N);
plot(1e3*T,G);
xlabel('time (msec)');
ylabel('Gauss/cm');
hold on;

% interpolate to 4us dwell time
dwell2 = 4e-6;          % gradient dwell (update) time, sec
T2 = dwell2/2:dwell2:(dwell*N-dwell/2);
G = interp1(dwell:dwell:(dwell*N), G, T2);
plot(1e3*T2,G,'ro')

% write to .mod file
N = size(G,1);
paramsint16 = [0 N 0 N];
mat2mod(0.002*ones(size(G)), 0*G, G, 0*G, 0*G, 90, 'readout.mod', 'GIRF calibration waveforms');

return;

% EOF
