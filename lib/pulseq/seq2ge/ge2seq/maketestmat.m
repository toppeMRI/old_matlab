function maketestmat

% get TOPPE timing CVs
fid = fopen('timing.txt', 'r', 'ieee-be');
s = fgets(fid);  % skip line
s = fscanf(fid, '%s ', 1);
start_core = fscanf(fid, '%d\n', 1);
s = fscanf(fid, '%s ', 1);
myrfdel = fscanf(fid, '%d\n', 1);
s = fscanf(fid, '%s ', 1);
daqdel = fscanf(fid, '%d\n', 1);
s = fscanf(fid, '%s ', 1);
timetrwait = fscanf(fid, '%d\n', 1);
s = fscanf(fid, '%s ', 1);
timessi = fscanf(fid, '%d\n', 1);
fclose(fid);
TPARAMS = [start_core myrfdel daqdel timetrwait timessi];

% Get module information.
fid = fopen('modules.txt', 'r', 'ieee-be');
s = fgets(fid);  % skip line
nmodules = fscanf(fid, '%d\n', 1);
s = fgets(fid);  % skip line
ishp = 1;  % index into shapeLibrary
for ii = 1:nmodules
	modules{ii}.fname = fscanf(fid, '%s ', 1);
	modules{ii}.dur = fscanf(fid, '%d ', 1);
	modules{ii}.hasRF = fscanf(fid, '%d ', 1);
	modules{ii}.hasDAQ = fscanf(fid, '%d\n', 1);
	%[desc rho th gx gy gz]	= readmod(modules{ii}.fname,false);
end
fclose(fid);

% Loop through scanloop.txt and get waveforms
d = readloop('scanloop.txt');
nt = size(d,1);    % number of startseq calls
max_pg_iamp = 2^15-2;
ii = 1;
for is = 62:79
	% get waveforms and delay for one row (one startseq call)
	[~,~,~,~,~,rho{ii},th{ii},gxwav{ii},gywav{ii},gzwav{ii},delay] = scansim(is, is, d, TPARAMS, false);
	freqOffset{ii} = d(is,15);
	phaseOffset{ii} = d(is,12)/max_pg_iamp*pi;
	module = modules{d(is,1)}
	hasRF{ii} = module.hasRF;
	hasDAQ{ii} = module.hasDAQ;
	ii = ii+1;
end
rasterTime = 4e-6;
save testseq rho th gxwav gywav gzwav freqOffset phaseOffset hasRF hasDAQ rasterTime




