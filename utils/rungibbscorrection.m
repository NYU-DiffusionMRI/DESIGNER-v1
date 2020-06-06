function rungibbscorrection(gibbsdir,root,DKIroot)
addpath(genpath(DKIroot));
addpath(gibbsdir);

nii = niftiread(fullfile(root,'dwidn.nii')); dwi = double(nii);

params = [1 3 20];
dwigc = unring(dwi,params);

info = niftiinfo(fullfile(root,'dwidn.nii'));
info.DataType = 'double';
niftiwrite(dwigc, fullfile(root,'dwigc.nii'), info);
end
