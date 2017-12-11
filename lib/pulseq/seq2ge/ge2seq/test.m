function test

% load sequence waveforms, freq/phase offsets, and rasterTime
load testseq

n = numel(gxwav);   % number of Blocks to be created

seq=mr.Sequence();  % initialize Pulseq sequence object

% create seq object blocks
for ii = 1:n
	rf = makeArbitraryRf(rf2hz(rho{ii}.*exp(1i*th{ii}),rasterTime),pi/4);
	gx = makeArbitraryGrad('x',g2hzcm(gxwav{ii},rasterTime));
	gy = makeArbitraryGrad('y',g2hzcm(gywav{ii},rasterTime));
	gz = makeArbitraryGrad('z',g2hzcm(gzwav{ii}/2,rasterTime));
	
	if hasRF{ii}
		seq.addBlock(rf,gx,gy,gz);
	elseif hasDAQ{ii}
		[numel(gxwav{ii}) numel(gx.waveform)]
		adc = makeAdc(numel(gxwav{ii}),'Dwell',rasterTime,'delay',0,'freqOffset',freqOffset{ii}, ...
									'phaseOffset',phaseOffset{ii});
		seq.addBlock(gx,gy,gz,adc);
	else
		seq.addBlock(gx,gy,gz);
	end

	%{
	if delay > 0
		del.delay = delay/1000;             % microseconds
		del.type = 'delay';
		%seq.addBlock(del);
	end
	%}
end

seq.plot();
seq.write('test.seq');

return;


function g = g2hzcm(g,rasterTime)
% convert gradient Gauss/cm to Hz/cm, and interpolate to GradRasterTime
gamma = 42.576e6;       % Hz/T
g = 1e-2 * g * gamma;   % Hz/cm
GradRasterTime = 1e-5;  % Pulseq gradient raster time
n = numel(g);
g = interp1(1:n,g,floor(GradRasterTime/rasterTime):GradRasterTime/rasterTime:n);
return;

function rf = rf2hz(rf,rasterTime) %,rasterTime)
% convert rf units from Gauss to Hz, and interpolate to RfRasterTime
gamma = 42.576e6;    % Hz/T
rf = 1e-4*rf*gamma;  % Hz
RfRasterTime = 1e-6;    % Pulseq RF raster time
rf = interp(rf,rasterTime/RfRasterTime);     % upsample from 4us to 1us
return;

	if 0
		% display waveforms
		subplot(221); plot(abs(rf.signal),'r'); title(sprintf('max = %f', max(abs(rf.signal))));
		subplot(222); plot(angle(rf.signal),'r');  title(sprintf('max = %f', max(angle(rf.signal))));
		subplot(2,2,[3 4]); plot(gx.waveform,'r'); hold on; plot(gy.waveform,'g'); plot(gz.waveform,'b'); hold off;
		title(sprintf('max = %f', max([gx.waveform gy.waveform gz.waveform])));
		input('press any key to continue');
		%keyboard
	end
