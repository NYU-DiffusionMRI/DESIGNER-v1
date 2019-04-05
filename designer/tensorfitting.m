function tensorfitting(root,outdir,detectoutliers,cumulants,dti,dki,wmti,fitconstraints,medianfilter,akc,DKIroot)

addpath(genpath(DKIroot));

maskex = exist(fullfile(root,'brain_mask.nii'),'file');
if maskex == 2
    nii = load_untouch_nii(fullfile(root,'brain_mask.nii')); mask = logical(nii.img);
end
nii = load_untouch_nii(fullfile(root,'dwi_designer.nii')); dwi = double(nii.img);
ndwis = size(dwi,4); pixdim5 = nii.hdr.dime.pixdim(5);
if maskex == 0
    mask = logical(ones(size(dwi,1),size(dwi,2),size(dwi,3)));
end

bvallist = dir(fullfile(root,'dwi_designer.bval')); bvaldir = bvallist.name;
bveclist = dir(fullfile(root,'dwi_designer.bvec')); bvecdir = bveclist.name;
bval = textread(fullfile(root,bvaldir)); bval = bval(:, 1:ndwis)'; bval = bval./1000;
bvec = textread(fullfile(root,bvecdir)); bvec = bvec(:, 1:ndwis)';
maxbval = max(bval);

detectoutliers = logical(detectoutliers);
medianfilter = logical(medianfilter);
dti = logical(dti);
dki = logical(dki);
wmti = logical(wmti);
akc = logical(akc);
if ischar(fitconstraints)
    conts = split(fitconstraints,',');
    constraints(1) = str2double(conts(1));
    constraints(2) = str2double(conts(2));
    constraints(3) = str2double(conts(3));
else
    constraints = [0,1,0];
end
cumulants = logical(cumulants);

if detectoutliers
    disp('...running IRWLLS')
    options.Excludeb0 = 0;
    outliers = irlls(dwi,mask,bvec,bval,[],options);
    nii.hdr.dime.dim(1) = 4;
    nii.hdr.dime.dim(5) = ndwis;
    nii.hdr.dime.pixdim(5) = pixdim5;
    nii.img = outliers;  save_untouch_nii(nii,fullfile(root,'irwlls_out.nii'));
    disp(['...fitting with constraints ',num2str(constraints)])
    [b0,dt,violMask] = dki_fit(dwi,[bvec,bval],mask,constraints,outliers,maxbval);
else
    disp(['...fitting with constraints ',num2str(constraints)])
    [b0,dt,violMask] = dki_fit(dwi,[bvec,bval],mask,constraints,[],maxbval);
end

if akc
    disp('...running AKC correction')
    akc_out = outlierdetection(dt,mask);
    akc_out(isnan(akc_out)) = 0;
    for v = 1:size(dt,4)
        dt_v = dt(:,:,:,v);
        dt_v(logical(akc_out)) = NaN;
        if verLessThan('matlab','9.2')
            dt_f = fillmissing(dt_v,'linear');
        else
            dt_f = fillmissing(dt_v,'movmedian',5);
        end
        dt(:,:,:,v) = dt_f;
    end
    
    nii.hdr.dime.dim(1) = 3;
    nii.hdr.dime.dim(5) = 1;
    nii.hdr.dime.pixdim(5) = 0;
    nii.img = akc_out;
    save_untouch_nii(nii,fullfile(root,'akc_out.nii'));
end

if cumulants
    dt = vectorize(dt, mask)
    D_apprSq = (sum(dt([1 4 6],:),1)/3).^2;
    dt(7:21,:) = dt(7:21,:) .* D_apprSq(ones(15,1),:);
    dt = vectorize(dt, mask)
end

disp('...getting DTI and DKI params')
[fa, md, rd, ad, fe, mk, rk, ak, kfa, mkt] = dki_parameters(dt,mask,violMask,medianfilter);
fe = cat(4,fa.*fe(:,:,:,1),fa.*fe(:,:,:,2),fa.*fe(:,:,:,3));

if dti
    disp('...saving DTI params')
    nii.hdr.dime.dim(1) = 3;
    nii.hdr.dime.dim(5) = 1;
    nii.hdr.dime.pixdim(5) = 0;
    nii.img = fa; nii.hdr.dime.glmax = max(fa(:)); save_untouch_nii(nii,fullfile(outdir,'fa.nii'));
    nii.img = md; nii.hdr.dime.glmax = max(md(:)); save_untouch_nii(nii,fullfile(outdir,'md.nii'));
    nii.img = rd; nii.hdr.dime.glmax = max(rd(:)); save_untouch_nii(nii,fullfile(outdir,'rd.nii'));
    nii.img = ad; nii.hdr.dime.glmax = max(ad(:)); save_untouch_nii(nii,fullfile(outdir,'ad.nii'));
    nii.img = b0; nii.hdr.dime.glmax = max(b0(:)); save_untouch_nii(nii,fullfile(outdir,'b0.nii'));
    nii.hdr.dime.dim(1) = 4;
    nii.hdr.dime.dim(5) = 3;
    nii.hdr.dime.pixdim(5) = pixdim5;
    nii.img = fe; save_untouch_nii(nii,fullfile(outdir,'fe.nii'));
    nii.hdr.dime.dim(5) = 21;
    nii.img = dt; save_untouch_nii(nii,fullfile(outdir,'dkt.nii'));
end
if dki
    disp('...saving DKI params')
    nii.hdr.dime.dim(1) = 3;
    nii.hdr.dime.dim(5) = 1;
    nii.hdr.dime.pixdim(5) = 0;
    nii.img = mk; nii.hdr.dime.glmax = max(mk(:)); save_untouch_nii(nii,fullfile(outdir,'mk.nii'));
    nii.img = rk; nii.hdr.dime.glmax = max(rk(:)); save_untouch_nii(nii,fullfile(outdir,'rk.nii'));
    nii.img = ak; nii.hdr.dime.glmax = max(ak(:)); save_untouch_nii(nii,fullfile(outdir,'ak.nii'));
    nii.img = kfa; nii.hdr.dime.glmax = max(ak(:)); save_untouch_nii(nii,fullfile(outdir,'kfa.nii'));
    nii.img = mkt; nii.hdr.dime.glmax = max(ak(:)); save_untouch_nii(nii,fullfile(outdir,'mkt.nii'));
end

if wmti
    disp('...getting WMTI params')
    nii.hdr.dime.dim(1) = 3;
    nii.hdr.dime.dim(5) = 1;
    nii.hdr.dime.pixdim(5) = 0;
    [awf, eas, ias] = wmti_parameters(dt, mask);
    disp('...saving WMTI params')
    nii.img = awf; nii.hdr.dime.glmax = max(awf(:)); save_untouch_nii(nii,fullfile(outdir,'awf.nii'));
    fields = fieldnames(ias);
    for ii=1:numel(fields)
        paramsii = getfield(ias, fields{ii});
        savename = fullfile(outdir, ['ias_', fields{ii}, '.nii']);
        nii.img = paramsii; nii.hdr.dime.glmax = max(paramsii(:)); save_untouch_nii(nii,savename);
    end
    fields = fieldnames(eas);
    for ii=1:numel(fields)
        paramsii = getfield(eas, fields{ii});
        savename = fullfile(outdir, ['eas_', fields{ii}, '.nii']);
        nii.img = paramsii; nii.hdr.dime.glmax = max(paramsii(:)); save_untouch_nii(nii,savename);
    end
end
    addpath(genpath(DKIroot));

    maskex = exist(fullfile(root,'brain_mask.nii'),'file');
    if maskex == 2
        nii = load_untouch_nii(fullfile(root,'brain_mask.nii')); mask = logical(nii.img);
    end
    nii = load_untouch_nii(fullfile(root,'dwi_designer.nii')); dwi = double(nii.img);
    ndwis = size(dwi,4); pixdim5 = nii.hdr.dime.pixdim(5);
    if maskex == 0
        mask = logical(ones(size(dwi,1),size(dwi,2),size(dwi,3)));
    end

    bvallist = dir(fullfile(root,'dwi_designer.bval')); bvaldir = bvallist.name;
    bveclist = dir(fullfile(root,'dwi_designer.bvec')); bvecdir = bveclist.name;
    bval = textread(fullfile(root,bvaldir)); bval = bval(:, 1:ndwis)'; bval = bval./1000;
    bvec = textread(fullfile(root,bvecdir)); bvec = bvec(:, 1:ndwis)';
    maxbval = max(bval);

    detectoutliers = logical(detectoutliers);
    dti = logical(dti);
    dki = logical(dki);
    wmti = logical(wmti);
    akc = logical(akc);
    if ischar(fitconstraints)
        conts = split(fitconstraints,',');
        constraints(1) = str2double(conts(1));
        constraints(2) = str2double(conts(2));
        constraints(3) = str2double(conts(3));
    else
        constraints = [0,1,0];
    end
    cumulants = logical(cumulants);

    if detectoutliers
        disp('...running IRWLLS')
        options.Excludeb0 = 0;
        outliers = irlls(dwi,mask,bvec,bval,[],options);
        nii.hdr.dime.dim(1) = 4;
        nii.hdr.dime.dim(5) = ndwis;
        nii.hdr.dime.pixdim(5) = pixdim5;
        nii.img = outliers;  save_untouch_nii(nii,fullfile(root,'irwlls_out.nii'));
        disp(['... fitting with constraints ',num2str(constraints)])
        [b0,dt,violMask] = dki_fit(dwi,[bvec,bval],mask,constraints,outliers,maxbval);
    else
        disp(['... fitting with constraints ',num2str(constraints)])
        [b0,dt,violMask] = dki_fit(dwi,[bvec,bval],mask,constraints,[],maxbval);
    end

    if akc
        disp('...running AKC correction')
        akc_out = outlierdetection(dt,mask);
        akc_out(isnan(akc_out)) = 0;
        for v = 1:size(dt,4)
            dt_v = dt(:,:,:,v);
            dt_v(logical(akc_out)) = NaN;
            if verLessThan('matlab','9.2')
                dt_f = fillmissing(dt_v,'linear');
            else
                dt_f = fillmissing(dt_v,'movmedian',5);
            end
            dt(:,:,:,v) = dt_f;
        end

        nii.hdr.dime.dim(1) = 3;
        nii.hdr.dime.dim(5) = 1;
        nii.hdr.dime.pixdim(5) = 0;
        nii.img = akc_out;
        save_untouch_nii(nii,fullfile(root,'akc_out.nii'));
    end

    if cumulants
        dt = vectorize(dt, mask)
        D_apprSq = (sum(dt([1 4 6],:),1)/3).^2;
        dt(7:21,:) = dt(7:21,:) .* D_apprSq(ones(15,1),:);
        dt = vectorize(dt, mask)
    end

    disp('...getting DTI and DKI params')
    [fa, md, rd, ad, fe, mk, rk, ak] = dki_parameters(dt,mask,violMask,medianfilter);
    fe = cat(4,fa.*fe(:,:,:,1),fa.*fe(:,:,:,2),fa.*fe(:,:,:,3));

    if dti
        disp('...saving DTI params')
        nii.hdr.dime.dim(1) = 3;
        nii.hdr.dime.dim(5) = 1;
        nii.hdr.dime.pixdim(5) = 0;
        nii.img = fa; nii.hdr.dime.glmax = max(fa(:)); save_untouch_nii(nii,fullfile(outdir,'fa.nii'));
        nii.img = md; nii.hdr.dime.glmax = max(md(:)); save_untouch_nii(nii,fullfile(outdir,'md.nii'));
        nii.img = rd; nii.hdr.dime.glmax = max(rd(:)); save_untouch_nii(nii,fullfile(outdir,'rd.nii'));
        nii.img = ad; nii.hdr.dime.glmax = max(ad(:)); save_untouch_nii(nii,fullfile(outdir,'ad.nii'));
        nii.img = b0; nii.hdr.dime.glmax = max(b0(:)); save_untouch_nii(nii,fullfile(outdir,'b0.nii'));
        nii.hdr.dime.dim(1) = 4;
        nii.hdr.dime.dim(5) = 3;
        nii.hdr.dime.pixdim(5) = pixdim5;
        nii.img = fe; save_untouch_nii(nii,fullfile(outdir,'fe.nii'));
        nii.hdr.dime.dim(5) = 21;
        nii.img = dt; save_untouch_nii(nii,fullfile(outdir,'dkt.nii'));
    end
    if dki
        disp('...saving DKI params')
        nii.hdr.dime.dim(1) = 3;
        nii.hdr.dime.dim(5) = 1;
        nii.hdr.dime.pixdim(5) = 0;
        nii.img = mk; nii.hdr.dime.glmax = max(mk(:)); save_untouch_nii(nii,fullfile(outdir,'mk.nii'));
        nii.img = rk; nii.hdr.dime.glmax = max(rk(:)); save_untouch_nii(nii,fullfile(outdir,'rk.nii'));
        nii.img = ak; nii.hdr.dime.glmax = max(ak(:)); save_untouch_nii(nii,fullfile(outdir,'ak.nii'));
    end

    if wmti
        disp('...getting WMTI params')
        nii.hdr.dime.dim(1) = 3;
        nii.hdr.dime.dim(5) = 1;
        nii.hdr.dime.pixdim(5) = 0;
        [awf, eas, ias] = wmti_parameters(dt, mask);
        disp('...saving WMTI params')
        nii.img = awf; nii.hdr.dime.glmax = max(awf(:)); save_untouch_nii(nii,fullfile(outdir,'awf.nii'));
        fields = fieldnames(ias);
        for ii=1:numel(fields)
            paramsii = getfield(ias, fields{ii});
            savename = fullfile(outdir, ['ias_', fields{ii}, '.nii']);
            nii.img = paramsii; nii.hdr.dime.glmax = max(paramsii(:)); save_untouch_nii(nii,savename);
        end
        fields = fieldnames(eas);
        for ii=1:numel(fields)
            paramsii = getfield(eas, fields{ii});
            savename = fullfile(outdir, ['eas_', fields{ii}, '.nii']);
            nii.img = paramsii; nii.hdr.dime.glmax = max(paramsii(:)); save_untouch_nii(nii,savename);
        end
    end

    function [s, mask] = vectorize(S, mask)
        if nargin == 1
            mask = ~isnan(S(:,:,:,1));
        end
        if ismatrix(S)
            n = size(S, 1);
            [x, y, z] = size(mask);
            s = NaN([x, y, z, n], 'like', S);
            for i = 1:n
                tmp = NaN(x, y, z, 'like', S);
                tmp(mask(:)) = S(i, :);
                s(:,:,:,i) = tmp;
            end
        else
            for i = 1:size(S, 4)
                Si = S(:,:,:,i);
                s(i, :) = Si(mask(:));
            end
        end
    end

%   Split and save dt into DT and KT for MUSC DKE_FT
disp('...saving diffusion and kurtosis tensors')
DT = dt(:,:,:,1:6);
KT = dt(:,:,:,7:21);

%   Reshape 4D tensors into [num_voxels x tensors]
DT = reshape(DT,[],6);
KT = reshape(KT,[],15);

%   Initialize post-transformation matrices
DT_U = zeros(size(DT));
KT_U = zeros(size(KT));

%   The two pipelines take in tensors in the following layout:
%       =============================
%           MUSC            Designer
%       --------------D--------------
%    1  |   D11               D11
%    2  |   D22               D12
%    3  |   D33               D13
%    4  |   D12               D22
%    5  |   D13               D23
%    6  |   D23               D33
%       --------------K--------------
%    1  |  W1111             W1111
%    2  |  W2222             W1112
%    3  |  W3333             W1113
%    4  |  W1112             W1122
%    5  |  W1113             W1123
%    6  |  W1222             W1133
%    7  |  W1333             W1222
%    8  |  W2223             W1223
%    9  |  W2333             W1233
%   10  |  W1122             W1333
%   11  |  W1133             W2222
%   12  |  W2233             W2223
%   13  |  W1123             W2233
%   14  |  W1223             W2333
%   15  |  W1233             W3000

%   Using the table above, the following transformation matrices tranform
%   designer tensors onto MUSC's tensors.
DT_trans = [1 4 6 2 3 5];
KT_trans = [1 11 15 2 3 7 10 12 14 4 6 13 5 8 9];

%   Apply transformation
for i = 1:length(DT_trans)
    DT_U(:,i) = DT(:,DT_trans(i));
end

for i = 1:length(KT_trans)
    KT_U(:,i) = KT(:,KT_trans(i));
end

clear DT KT

%   Transpose to obtain [tensors x num_voxels]
DT = DT_U'; KT = KT_U';

%   Save tensors
save(fullfile(outdir,'DT.mat'),'DT');
save(fullfile(outdir,'KT.mat'),'KT');

%   Save violation mask
if ~exist(fullfile(outdir,'QC'))
    mkdir(fullfile(outdir,'QC'));
else
    ;
end
nii.img = violMask; nii.hdr.dime.glmax = max(b0(:)); save_untouch_nii(nii,fullfile(outdir,'QC','violation_map.nii'));
end

