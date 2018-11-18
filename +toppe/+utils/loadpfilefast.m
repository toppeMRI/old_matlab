function dat = loadpfilefast(pfile)
% Loads in pfile data. Like loadpfile, but hopefully faster
% Doesn't depend on loaddat_ge
% Checks that pfile matches header size

import toppe.utils.*

%% Loadpfile code
fid = fopen(pfile,'r','l');
ver = fread(fid,1,'float32');
str = num2str(ver);
rdbm_rev = str2double(str);
fseek(fid,0,'bof');                 % NB!
rdb_hdr = read_rdb_hdr(fid,rdbm_rev);


%% Header parameters
ndat    = rdb_hdr.frame_size;
nslices = rdb_hdr.nslices;
ptsize  = rdb_hdr.point_size;                    % Either 2 (data stored in short int format, int16) or 4 (extended precision)
nechoes = rdb_hdr.nechoes;
nviews  = rdb_hdr.nframes;
ncoils  = rdb_hdr.dab(2)-rdb_hdr.dab(1)+1;

%% Calculate size of data chunks. See pfilestruct.jpg, and rhrawsize calculation in .e file.
echores  = ndat*(nviews+1);                 % number of data points per 'echo' loaddab slot. Includes baseline (0) view.
sliceres = nechoes*echores;                 % number of data points per 'slice'
coilres  = nslices*sliceres;                % number of data points per receive coil

pfilesize = rdb_hdr.off_data + 2*ptsize*ncoils*nslices*nechoes*(nviews+1)*ndat;   % this should match the Pfile size exactly
pfilename=dir(pfile);

if pfilesize ~= pfilename.bytes
    fprintf('Expected %0.1fMB file but read in %0.1fMB file.\n',pfilesize,pfilename.bytes/1e6)
    error('Pfile size/header mismatch');
end

fprintf(1,'\nndat = %d, nslices = %d, nechoes = %d, nviews = %d, ncoils = %d\n', ndat, nslices, nechoes, nviews, ncoils);

%% Read data from file
datr = int16(zeros(ndat,ncoils,nslices-1,nechoes,nviews));
dati = datr;
for icoil = 1:ncoils
    for islice = 2:nslices   % skip first slice (sometimes contains corrupted data)
        for iecho = 1:nechoes % Echo starts at 1
            for iview = 1:nviews
                offsetres = (icoil-1)*coilres + (islice-1)*sliceres + (iecho-1)*echores + iview*ndat;
                offsetbytes = 2*ptsize*offsetres;
                fseek(fid, rdb_hdr.off_data+offsetbytes, 'bof');
                dtmp = fread(fid, 2*ndat, 'int16=>int16');
                datr(:,icoil,islice-1,iecho,iview) = dtmp(1:2:end); %Real data
                dati(:,icoil,islice-1,iecho,iview) = dtmp(2:2:end); %Imag
            end
        end
    end
end
dat = double(complex(datr,dati)); %Combine real+imag
fclose(fid);

return;
