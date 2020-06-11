function runrpg(pf,dim,root,DKI_root)

    %% Removal of Partial-fourier induced Gibbs ringing (RPG)
    % Author: Hong-Hsi Lee, February 2020
    %
    % C++ codes are modified from the code based on local subvoxel-shifts by
    % Kellner et al. MRM 2016 (doi: 10.1002/mrm.26054)
    % (https://bitbucket.org/reisert/unring/src/master/)

%     clear
    restoredefaultpath

    % Directory to this script
    addpath(genpath(fullfile(DKI_root,'utils/lib')));
    addpath(genpath(fullfile(DKI_root,'utils')));

    %% Compile mex files
    %rpg = rpgdegibbs();
    %rpg.compilefiles(fullfile(DKI_root,'utils/lib'));

    %% degibbs

    % load img w gibbs
    img_path = fullfile(root,'working.nii');
    out = fullfile(root,'working_rpg.nii');
    vol = niftiread(img_path);
    vol_ = vol;

    % Remove Gibbs ringing
    % partial-fourier dimension, 2 (phase encoding), 1 (frequency encoding)
    % partial-fourier factor, pf
    c = regexp(pf,'[0-9]+','match');
    num = str2double(c{1});
    den = str2double(c{2});
    pf = num/den;
    
    rpg = rpgdegibbs(); % rpgdegibbs_s() removes more Gibbs ringing but introduces more smoothing.   
    P_pf_dg = rpg.degibbs(abs(vol_),dim,pf);

    info = niftiinfo(img_path);
    niftiwrite(P_pf_dg, out, info);

end
