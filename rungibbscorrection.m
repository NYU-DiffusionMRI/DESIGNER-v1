function rungibbscorrection(gibbsdir,root)
addpath(genpath(pwd));
addpath(gibbsdir);

nii = load_untouch_nii(fullfile(root,'dwidn.nii')); dwi = double(nii.img);

params = [1 3 20];
dwigc = unring(dwi,params);
nii.img = dwigc;
nii.hdr.dime.glmax = max(dwigc(:));
save_untouch_nii(nii,fullfile(root,'dwigc.nii'));

end