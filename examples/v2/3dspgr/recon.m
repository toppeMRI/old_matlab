function [obj dat] = recon(pfile)
% %
% 3DFT recon

dat = loadpfile(pfile);               % [ndat ncoils nslices nechoes nviews]
d = dat(256:1063,:,2:end,:,:);        % skip dabslice=0

ncoils = size(d,2); % 32; $numel(1:16:32);

for ii=1:ncoils                           
	d1 = squeeze(d(:,ii,:,1,:));   
	d1 = permute(d1,[1 3 2]);         % [ndat nviews nslices]
	obj(:,:,:,ii) = ift3(d1,false);
	%figure; im(obj(:,:,:,ii)); title(num2str(ii)); 
end

obj = sqrt(sum(abs(obj).^2,ndims(obj)));
figure; im(obj);

