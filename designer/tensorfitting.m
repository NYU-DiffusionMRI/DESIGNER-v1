function tensorfitting(root,outdir,detectoutliers,dti,dki,wmti,fitconstraints,akc,DKIroot)
addpath(genpath(DKIroot));

nii = load_untouch_nii(fullfile(root,'brain_mask.nii')); mask = logical(nii.img);
nii = load_untouch_nii(fullfile(root,'dwi_designer.nii')); dwi = double(nii.img);
ndwis = size(dwi,4); pixdim5 = nii.hdr.dime.pixdim(5);

bvallist = dir(fullfile(root,'dwi_designer.bval')); bvaldir = bvallist.name;
bveclist = dir(fullfile(root,'dwi_designer.bvec')); bvecdir = bveclist.name;
bval = textread(fullfile(root,bvaldir)); bval = bval(:, 1:ndwis)'; bval = bval./1000;
bvec = textread(fullfile(root,bvecdir)); bvec = bvec(:, 1:ndwis)';
maxbval = max(bval);

%f = find(bval~=.25)
%dwi = dwi(:,:,:,f);
%bvec = bvec(f,:);
%bval = bval(f);

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
    constraints = 0;
end

if detectoutliers
    disp('...running IRWLLS')
    options.Excludeb0 = 0;
    outliers = irlls(dwi,mask,bvec,bval,[],options);
    nii.hdr.dime.dim(1) = 4;
    nii.hdr.dime.dim(5) = ndwis;
    nii.hdr.dime.pixdim(5) = pixdim5;  
    nii.img = outliers;  save_untouch_nii(nii,fullfile(root,'irwlls_out.nii'));
    if sum(constraints) == 0
        [b0,dt] = dki_fit(dwi,[bvec,bval],mask,[0,1,0],outliers,maxbval);
    else
	disp('...running constrained fit')
        [b0,dt] = dki_fit(dwi,[bvec,bval],mask,constraints,outliers,maxbval);
    end
else
    if sum(constraints) == 0
        [b0,dt] = dki_fit(dwi,[bvec,bval],mask,[0,1,0],[],maxbval);
    else
	disp('...running constrained fit')
	[b0,dt] = dki_fit(dwi,[bvec,bval],mask,constraints,[],maxbval);
    end
end

if akc
    akc_out = outlierdetection(dt);
    akc_out(isnan(akc_out)) = 0;
    for v = 1:size(dt,4)
        dt_v = dt(:,:,:,v);
        dt_v(logical(akc_out)) = NaN;
        dt_f = fillmissing(dt_v,'movmedian',5);
        dt(:,:,:,v) = dt_f;
    end

    nii.hdr.dime.dim(1) = 3;
    nii.hdr.dime.dim(5) = 1;
    nii.hdr.dime.pixdim(5) = 0;
    nii.img = akc_out;
    save_untouch_nii(nii,fullfile(root,'akc_out.nii'));
end

[fa, md, rd, ad, fe, mk, rk, ak] = dki_parameters(dt,mask);

if dti
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
    nii.hdr.dime.dim(1) = 3;
    nii.hdr.dime.dim(5) = 1;
    nii.hdr.dime.pixdim(5) = 0;
    nii.img = mk; nii.hdr.dime.glmax = max(mk(:)); save_untouch_nii(nii,fullfile(outdir,'mk.nii'));
    nii.img = rk; nii.hdr.dime.glmax = max(rk(:)); save_untouch_nii(nii,fullfile(outdir,'rk.nii'));
    nii.img = ak; nii.hdr.dime.glmax = max(ak(:)); save_untouch_nii(nii,fullfile(outdir,'ak.nii'));
end

if wmti
    nii.hdr.dime.dim(1) = 3;
    nii.hdr.dime.dim(5) = 1;
    nii.hdr.dime.pixdim(5) = 0;
    [awf, eas, ias] = wmti_parameters(dt, mask);
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

end
