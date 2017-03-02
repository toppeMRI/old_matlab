function g = makebalanced(g,mxs)
% add gradient refocusing lobe at end 
% g = [gx gy], size [N 2] (G/cm)
% mxs = G/cm/sec [default = 15000]

if ~exist('mxs', 'var')
	mxs = 15000;
end


% add refocusing lobe to gx and gy
gx = sub_addlobe(g(:,1),mxs/sqrt(2)); % sqrt(2) since we want to keep combination of gx and gy within slew rate limit
gy = sub_addlobe(g(:,2),mxs/sqrt(2));

% make sure gx and gy have same length
ll = max(length(gx),length(gy));
if ll > length(gx)
	gx= [gx(:); zeros(ll-length(gx),1)];
else
	gy= [gy(:); zeros(ll-length(gy),1)];
end

g = [gx(:) gy(:)];

% g = makeevenlength(g);

return;



function g = sub_addlobe(g,mxs)

dt = 4e-6;           % sec
mxg = 5.0/sqrt(2);           % G/cm
%mxs = 10e3/sqrt(2);        % G/cm/sec

% make sure g ends at 0
if g(end) ~= 0
	s = sign(g(end))*mxs*dt; % max change in g per sample (G/cm)
	g = [g(:); g(end)*ones(2,1);[g(end):-s:0]'];
end

area = dt*sum(g);    % G/cm*sec

gtrap = trapwave(-area,dt,mxg,mxs);
areatrap = dt*sum(gtrap);

g = [g(:); gtrap(:)]; %/(-areatrap)*area];

sum(g);

return;


