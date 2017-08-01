function designer(root, varargin)

% DESIGNER pipeline for the processing of diffusion weighted MR data and
% the estimation of diffusion kurtosis tensor parameters and white matter
% integrity metrics. 
%
% please read: https://github.com/NYU-DiffusionMRI/Diffusion-Kurtosis-Imaging/wiki/Diffusion-Kurtosis-Imaging:-some-background
% for more information in the steps included in this pipeline.
%
% the pipeline required the installation of MRtrix3.0, FSL (latest verion),
% and "unring". please read
% https://github.com/NYU-DiffusionMRI/Diffusion-Kurtosis-Imaging/wiki/DESIGNER-installation
% for more info on those software packages.
%
%
% how to use this pipeline?
% designer(targetfolder, dwifiles, [PE, [reversed_encoded_files]])
% 
% dwifiles =  list of filenames that needs to be included in the DKI/WMTI
%             processing. For now, we assume all those images are acquired with the
%             same phase encoding.
% PE = phase encoding direction: 'AP', 'PA', 'LR', 'RL', ... of the
%      dwifiles
% reversed_encoded_files = filename of image that is acquired with reversed
%                           phase encoding. Important: only the first image is considerd in the processing. 
%
% jelle.veraart@nyumc.org



params.gibbscorrect = [1 3 20];
params.denoising.kernel = [5, 5, 5];

params.smoothing.enable = true;
params.smoothing.kernel = [5, 5];
params.smoothing.width = 1.2;

params.fit.maxbval = 3;
params.fit.constraints = [0 1 0];

params.b1correction = true;
params.ricianbiascorrection = true;



nargin  = numel(varargin);
if nargin == 1
    pefiles = varargin{1};
    pe = 'AP';
elseif nargin == 2
    pefiles = varargin{1};
    pe = varargin{2};
elseif nargin == 3
    pefiles = varargin{1};
    pe = varargin{2};   
    rpefile = varargin{3};
else
    return
end


pwdir = pwd;
cd(root)

%% concatinate data
mkdir('./Analysis/sorted/')
mkdir('./Analysis/processed/')
mkdir('./Analysis/parameters/')


cmd = 'mrcat -axis 3 ';

bvals = [];
bvecs = [];
index = [];

if iscell(pefiles);
    N = numel(pefiles);
else
    N = 1;
end
    
for i = 1:N
  if iscell(pefiles) 
    [rt, filename] = fileparts(pefiles{i});
  else
    [rt, filename] = fileparts(pefiles);
  end
      
  if isempty(rt)
      rt = './';
  end
  info = load_nii_hdr(fullfile(rt, [filename, '.nii']));   
  
  ndwis = info.dime.dim(5);
  
  if exist(fullfile(rt, [filename, '.bval']), 'file') & exist(fullfile(rt, [filename, '.bvec']), 'file')
     bval = textread(fullfile(rt, [filename, '.bval'])); bval = bval(1:ndwis)'; 
     bvec = textread(fullfile(rt, [filename, '.bvec'])); bvec = bvec(:, 1:ndwis)';
  end

  if params.b1correction
    system(['dwibiascorrect -fslgrad ', fullfile(rt, [filename, '.bvec']), ' ', fullfile(rt, [filename, '.bval']), ' -fsl ', fullfile(rt, [filename, '.nii']), ' ./Analysis/sorted/bfc', num2str(i), '.nii'])
    cmd = [cmd, ['./Analysis/sorted/bfc', num2str(i), '.nii ']];
  else
    cmd = [cmd, fullfile(rt, [filename, '.nii']), ' '];
  end
  bvals = cat(1, bvals, bval);
  bvecs = cat(1, bvecs, bvec);
  index = cat(1, index, i*ones(size(bval)));
end

bvals = bvals';
bvecs = bvecs';
save('./Analysis/sorted/dwi.bvec', 'bvecs', '-ascii');
save('./Analysis/sorted/dwi.bval', 'bvals', '-ascii');
save('./Analysis/sorted/dwi.index', 'index', '-ascii');

cmd = [cmd, './Analysis/sorted/dwi.nii'];
if N>1
    system(cmd);
else
    system(['cp ', fullfile(rt, [filename, '.nii']), ' ./Analysis/sorted/dwi.nii']);
end


system('dwiextract -bzero -fslgrad ./Analysis/sorted/dwi.bvec ./Analysis/sorted/dwi.bval ./Analysis/sorted/dwi.nii ./Analysis/sorted/PEb0s.nii')
system('mrconvert -coord 3 0 ./Analysis/sorted/PEb0s.nii ./Analysis/sorted/PEb0.nii')

if exist('rpefile', 'var')
     if iscell(rpefile)
        [rt, filename] = fileparts(rpefile{1});
     else
        [rt, filename] = fileparts(rpefile);
     end
     if isempty(rt)
      rt = './';
     end
     
     system(['cp ', fullfile(rt, [filename, '.nii']), ' ./Analysis/sorted/rPEb0s.nii']);
     
     info = load_nii_hdr('./Analysis/sorted/rPEb0s.nii');
     if info.dime.dim(5) > 1
         system('mrconvert -coord 3 0 ./Analysis/sorted/rPEb0s.nii ./Analysis/sorted/rPEb0.nii')
     else
         system(['cp ./Analysis/sorted/rPEb0s.nii ./Analysis/sorted/rPEb0.nii']);
     end

     system('mrcat ./Analysis/sorted/PEb0.nii ./Analysis/sorted/rPEb0.nii ./Analysis/sorted/EPrPEb0.nii')
end

%% denoising
system(['dwidenoise -extent ', num2str(params.denoising.kernel(1)), ',', num2str(params.denoising.kernel(2)), ',', num2str(params.denoising.kernel(3)), ' -noise ./Analysis/sorted/sigma.nii ./Analysis/sorted/dwi.nii ./Analysis/processed/dn.nii']);

%% Gibbs ringing

try 
    nii = load_untouch_nii('./Analysis/processed/dn.nii'); dwi = single(nii.img);
    gc = unring(dwi, params.gibbscorrect);
    nii.img = gc;
    nii.hdr.dime.glmax = max(gc(:));
    nii.hdr.dime.datatype = 16;
    nii.hdr.dime.bitpix =  32;
    save_untouch_nii(nii, './Analysis/processed/gc.nii');
catch
    system(['unring ./Analysis/processed/dn.nii ./Analysis/processed/gc.nii nii -minW ', num2str(params.gibbscorrect(1)) ,' -maxW ', num2str(params.gibbscorrect(2)),' -nsh ', num2str(params.gibbscorrect(3))]);
end

%% topup + eddy, through MRTRIX3.0 

if exist('rpefile', 'var')
    system(['dwipreproc ./Analysis/processed/gc.nii ./Analysis/processed/ec.nii -eddy_options "--repol --data_is_shelled " -rpe_pair -se_epi ./Analysis/sorted/EPrPEb0.nii -pe_dir ', pe, ' -fslgrad ./Analysis/sorted/dwi.bvec ./Analysis/sorted/dwi.bval'])
    %system('dwipreproc  -rpe_pair ./Analysis/sorted/PEb0.nii ./Analysis/sorted/rPEb0.nii -fslgrad ./Analysis/sorted/dwi.bvec ./Analysis/sorted/dwi.bval ', pe, ' ./Analysis/processed/gc.nii ./Analysis/processed/ec.nii')
else
    system(['dwipreproc ./Analysis/processed/gc.nii ./Analysis/processed/ec.nii -eddy_options "--repol --data_is_shelled " -rpe_none -pe_dir ', pe, ' -fslgrad ./Analysis/sorted/dwi.bvec ./Analysis/sorted/dwi.bval'])
end

%% smooth data (after selecting/excluding CSF) + fitting
system('dwiextract -bzero -fslgrad ./Analysis/sorted/dwi.bvec ./Analysis/sorted/dwi.bval ./Analysis/processed/ec.nii ./Analysis/processed/b0s.nii')
system('mrmath -axis 3 ./Analysis/processed/b0s.nii mean ./Analysis/processed/b0.nii');

system('bet ./Analysis/processed/b0.nii ./Analysis/processed/b0.nii -m -f 0.25');
if exist('./Analysis/processed/b0.nii.gz', 'file')
    system('gunzip -f ./Analysis/processed/b0.nii.gz')
    system('gunzip -f ./Analysis/processed/b0_mask.nii.gz')

end
system('fast -n 4 ./Analysis/processed/b0.nii')
if exist('./Analysis/processed/b0_pve_0.nii.gz', 'file')
    system('gunzip -f ./Analysis/processed/b0_pve_0.nii.gz')
    system('gunzip -f ./Analysis/processed/b0_pve_1.nii.gz')
    system('gunzip -f ./Analysis/processed/b0_pve_2.nii.gz')
    system('gunzip -f ./Analysis/processed/b0_pve_3.nii.gz')
end

nii = load_untouch_nii('./Analysis/processed/b0.nii'); b0 = nii.img;

for i = 0:3
    nii = load_untouch_nii(['./Analysis/processed/b0_pve_', num2str(i),'.nii']);     
    avgS(i+1) = prctile(b0(nii.img(:)>0.95), 95);
end
[~,pos] = max(avgS);


nii = load_untouch_nii('./Analysis/sorted/sigma.nii'); sigma= single(nii.img);
nii = load_untouch_nii('./Analysis/processed/b0_mask.nii'); mask= single(nii.img)>0;
nii = load_untouch_nii(['./Analysis/processed/b0_pve_', num2str(pos-1),'.nii']); csfmask= single(nii.img)>0.95;
nii = load_untouch_nii('./Analysis/processed/ec.nii'); img = single(nii.img);



if params.smoothing.enable
    img = smoothing(img, params.smoothing.kernel, params.smoothing.width, csfmask);
    nii.img = img;
    save_untouch_nii(nii, './Analysis/processed/sm.nii');
end

%% Rician bias correction

if params.ricianbiascorrection
    for i = 1:size(img, 4)
        img(:,:,:,i) = sqrt(abs(img(:,:,:,i).^2 - sigma.^2));
    end
    nii.img = img;
    save_untouch_nii(nii, './Analysis/processed/rc.nii');
end
%% fitting
bval = textread('./Analysis/sorted/dwi.bval')';  
bvec = textread('./Analysis/sorted/dwi.bvec')';

options.Excludeb0 = 0;
outliers = irlls(img,mask,bvec,bval,[],options);          

[~, dt] = dki_fit(img, [bvec, bval], mask, params.fit.constraints, outliers, params.fit.maxbval);
[fa, md, rd, ad, fe, mk,  rk, ak] = dki_parameters(dt, mask);
[awf, eas, ias] = wmti_parameters(dt, mask);

%% save data

nii = load_untouch_nii('./Analysis/processed/ec.nii');

nii.hdr.dime.dim(5) = 1;
nii.hdr.dime.datatype = 16;
nii.hdr.dime.bitpix =  32;
nii.hdr.dime.scl_slope = 1;
nii.hdr.dime.scl_inter = 0;

nii.img = feval('single',fa); save_untouch_nii(nii, './Analysis/parameters/fa.nii');
nii.img = feval('single',md); save_untouch_nii(nii, './Analysis/parameters/md.nii');
nii.img = feval('single',rd); save_untouch_nii(nii, './Analysis/parameters/rd.nii');
nii.img = feval('single',ad); save_untouch_nii(nii, './Analysis/parameters/ad.nii');
nii.img = feval('single',mk); save_untouch_nii(nii, './Analysis/parameters/mk.nii');
nii.img = feval('single',rk); save_untouch_nii(nii, './Analysis/parameters/rk.nii');
nii.img = feval('single',ak); save_untouch_nii(nii, './Analysis/parameters/ak.nii');
nii.img = feval('single',awf); save_untouch_nii(nii, './Analysis/parameters/awf.nii');

nii.img = feval('single', ias.da1); save_untouch_nii(nii, './Analysis/parameters/Da_1.nii');
nii.img = feval('single', ias.da2); save_untouch_nii(nii, './Analysis/parameters/Da_2.nii');
nii.img = feval('single', ias.da3); save_untouch_nii(nii, './Analysis/parameters/Da_3.nii');
nii.img = feval('single', ias.da_perp); save_untouch_nii(nii, './Analysis/parameters/Daperp.nii');

nii.img = feval('single', eas.de1); save_untouch_nii(nii, './Analysis/parameters/De_1.nii');
nii.img = feval('single', eas.de2); save_untouch_nii(nii, './Analysis/parameters/De_2.nii');
nii.img = feval('single', eas.de3); save_untouch_nii(nii, './Analysis/parameters/De_3.nii');
nii.img = feval('single', eas.de_perp); save_untouch_nii(nii, './Analysis/parameters/Deperp.nii');
nii.img = feval('single', eas.tort); save_untouch_nii(nii, './Analysis/parameters/tort.nii');

cd(pwdir)
end
