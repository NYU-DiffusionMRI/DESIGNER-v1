%% Clearing
clc; close all; clear all;

%% Load Paths
folder = '/Users/sid/Downloads/Median_Tests/HARDI5-Output/DKE';
dwi_Path = fullfile(folder,'dwi_designer.nii');
bval_Path = fullfile(folder,'dwi_designer.bval');
bvec_Path = fullfile(folder,'dwi_designer.bvec');
mask_Path = fullfile(folder,'brain_mask.nii');

%% Read Files
%   Image
hdr = niftiinfo(dwi_Path);
dwi = niftiread(hdr);
brainmask = logical(niftiread(mask_Path));

%   BVEC and BVAL
bval = load(bval_Path);
bvec = load(bvec_Path);
bval = ceil(bval);

%% Form Indexes
b0_idx = find(bval == 0);
b1_idx = find(bval == 1000);
b2_idx = find(bval == 2000);

%% Average B0s
b0_mean = dwi(b0_idx(1));
for i = 2:length(b0_idx)
b0_mean = b0_mean + dwi(:,:,:,b0_idx(i));
end
b0_mean = b0_mean / length(b0_idx);

%% Replace Mean and Modify BVAL/BVEC
dwi(:,:,:,1) = b0_mean;

dwi(:,:,:,b0_idx(2:end)) = [];
for i = 1:size(dwi,4)
    dwi(:,:,:,i) = dwi(:,:,:,i) .* brainmask;
end
hdr.ImageSize = size(dwi);
niftiwrite(dwi,dwi_Path,hdr);

bval(b0_idx(2:end)) = [];
bvec(:,b0_idx(2:end)) = [];
[p,~,~] = fileparts(bvec_Path);
Gradient1 = bvec';
Gradient1 = Gradient1(b1_idx,:);
save(fullfile(p,'gradient_dke.txt'),'Gradient1','-ASCII');

