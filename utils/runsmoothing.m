function runsmoothing(root,width,DKIroot)
addpath(genpath(DKIroot));

nii = load_untouch_nii(fullfile(root,'CSFmask.nii')); csf = logical(nii.img);
nii = load_untouch_nii(fullfile(root,'dwibc.nii')); dwi = double(nii.img);
kernelsize = [5 5];
width = str2double(width);

dwism = smoothing(dwi,kernelsize,width,csf);

nii.img = dwism;
nii.hdr.dime.glmax = max(dwism(:));
save_untouch_nii(nii,fullfile(root,'dwism.nii'));

end
