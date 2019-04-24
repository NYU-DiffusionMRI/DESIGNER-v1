function des2dke(inDir)
% DES2DKE  Convert designer output to DKE inputs
%   ADDME('input directory') creates a folder called 'DKE' within
%   Designer's output folder and prepares all files for DKE.
%
%   Author: Siddhartha Dhiman
%   Email:  dhiman@musc.edu

%% Load Paths
dwi_Path = fullfile(inDir,'dwi_designer.nii');
bval_Path = fullfile(inDir,'dwi_designer.bval');
bvec_Path = fullfile(inDir,'dwi_designer.bvec');
mask_Path = fullfile(inDir,'brain_mask.nii');
dke_Path = fullfile(inDir,'DKE');
mkdir(dke_Path);

%% Read Files
fprintf('1: Reading Files\n');
%   Image
hdr = niftiinfo(dwi_Path);
dwi = niftiread(hdr);
fprintf('\tA:...loaded image\n');
brainmask = logical(niftiread(mask_Path));
fprintf('\tB:...loaded brainmask\n');

%   BVEC and BVAL
bval = load(bval_Path);
fprintf('\tC:...loaded BVAL\n');
bvec = load(bvec_Path);
bval = ceil(bval);
fprintf('\tD:...loaded BVEC\n');

%% Form Indexes
b0_idx = find(bval == 0);
b1_idx = find(bval == 1000);
b2_idx = find(bval == 2000);

%% Average B0s
fprintf('2: Computing B0 Means\n');
b0_mean = dwi(b0_idx(1));
for i = 2:length(b0_idx)
    b0_mean = b0_mean + dwi(:,:,:,b0_idx(i));
end
b0_mean = b0_mean / length(b0_idx);

%% Replace Mean and Modify BVAL/BVEC
% Update path variables
dwi_Path = fullfile(dke_Path,'dke.nii');
bval_Path = fullfile(dke_Path,'dke.bval');
bvec_Path = fullfile(dke_Path,'dke.bvec');

fprintf('3: Writing Files\n');
dwi(:,:,:,1) = b0_mean;

dwi(:,:,:,b0_idx(2:end)) = [];
for i = 1:size(dwi,4)
    dwi(:,:,:,i) = dwi(:,:,:,i) .* brainmask;
end
hdr.ImageSize = size(dwi);
fprintf('\tA:...writing image\n');
niftiwrite(dwi,dwi_Path,hdr);

bval(b0_idx(2:end)) = [];
bvec(:,b0_idx(2:end)) = [];
fprintf('B:...writing BVAL\n');
dlmwrite(bval_Path,bval,',');
fprintf('C:...writing BVAL\n');
dlmwrite(bvec_Path,bvec,',');
b1_idx_new = find(bval == 1000); %indices have changed, update


[p,~,~] = fileparts(bvec_Path);
Gradient1 = bvec';
Gradient1 = Gradient1(b1_idx_new,:);
fprintf('\tD:...writing gradient\n');
save(fullfile(p,'gradient_dke.txt'),'Gradient1','-ASCII');

%% Create DKE Parameter File
fprintf('4: Creating parameter files\n');
fid=fopen('dke_parameters.txt'); %Original file
fout=fullfile(dke_Path,'dke_parameters.txt');% new file

fidout=fopen(fout,'w');

while(~feof(fid))
    s=fgetl(fid);
    s=strrep(s,'dir-sub-changeme',dke_Path); %s=strrep(s,'A201', subject_list{i}) replace subject
    s=strrep(s,'ndir = changeme',sprintf('ndir = %d',length(b1_idx_new)));
    s=strrep(s,'fn-gradients-changeme',fullfile(dke_Path,'gradient_dke.txt'));
    fprintf(fidout,'%s\n',s);
end
fclose(fid);
fclose(fidout);
fprintf('.....Completed.....\n');
