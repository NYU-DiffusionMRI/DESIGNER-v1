function runsmoothing(root,width,DKIroot)
addpath(genpath(DKIroot));

nii = niftiread(fullfile(root,'CSFmask.nii')); csf = logical(nii);
info = niftiinfo(fullfile(root,'dwibc.nii'));
dwi = double(niftiread(info));

kernelsize = [5 5];
width = str2double(width);

dwism = smoothing(dwi,kernelsize,width,csf);
info.Datatype = 'double';
niftiwrite(dwism, fullfile(root,'dwism.nii'), info);
end
