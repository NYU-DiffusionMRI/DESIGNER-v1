function tensorfitting(root,outdir,detectoutliers,cumulants,dti,dki,wmti,fitconstraints,akc,DKIroot)

    addpath(genpath(DKIroot));

    maskex = exist(fullfile(root,'brain_mask.nii'),'file');
    if maskex == 2
        nii = niftiread(fullfile(root,'brain_mask.nii')); mask = logical(nii);
    end
    nii = niftiread(fullfile(root,'dwi_designer.nii')); dwi = double(nii);
    info = niftiinfo(fullfile(root,'dwi_designer.nii'));
    ndwis = size(dwi,4); %pixdim5 = nii.hdr.dime.pixdim(5);
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
        
        info.Datatype = 'double';
        niftiwrite(outliers, fullfile(root,'irwlls_out.nii'), info);
        disp(['... fitting with constraints ',num2str(constraints)])
        [b0,dt] = dki_fit(dwi,[bvec,bval],mask,constraints,outliers,maxbval);
    else
        disp(['... fitting with constraints ',num2str(constraints)])
        [b0,dt] = dki_fit(dwi,[bvec,bval],mask,constraints,[],maxbval);
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
        
        info.ImageSize = info.ImageSize(1:3);
        info.PixelDimensions = info.PixelDimensions(1:3);
        niftiwrite(akc_out, fullfile(root,'akc_out.nii'), info);
    end

    if cumulants
        dt = vectorize(dt, mask);
        D_apprSq = (sum(dt([1 4 6],:),1)/3).^2;
        dt(7:21,:) = dt(7:21,:) .* D_apprSq(ones(15,1),:);
        dt = vectorize(dt, mask);
    end

    disp('...getting DTI and DKI params')
    [fa, md, rd, ad, fe, mk, rk, ak] = dki_parameters(dt,mask);
    fe = cat(4,fa.*fe(:,:,:,1),fa.*fe(:,:,:,2),fa.*fe(:,:,:,3));
    
    if dti
        disp('...saving DTI params')
        info.ImageSize = info.ImageSize(1:3);
        info.PixelDimensions = info.PixelDimensions(1:3);
        info.Datatype = 'double';
        
        info.DisplayIntensityRange = [0 1];
        niftiwrite(fa, fullfile(outdir,'fa.nii'), info);
        info.DisplayIntensityRange = [0 3];
        niftiwrite(md, fullfile(outdir,'md.nii'), info);
        niftiwrite(rd, fullfile(outdir,'rd.nii'), info);
        niftiwrite(ad, fullfile(outdir,'ad.nii'), info);
        info.DisplayIntensityRange = [0 0];
        niftiwrite(b0, fullfile(outdir,'b0.nii'), info);
        info.ImageSize = cat(2, info.ImageSize, 3);
        info.PixelDimensions = cat(2, info.PixelDimensions, NaN);
        niftiwrite(fe, fullfile(outdir,'fe.nii'), info);
        info.ImageSize(4) = 21;
        niftiwrite(dt, fullfile(outdir,'dkt.nii'), info);
    end
    if dki
        disp('...saving DKI params')
        info.ImageSize = info.ImageSize(1:3);
        info.PixelDimensions = info.PixelDimensions(1:3);
        info.DisplayIntensityRange = [0 3];
        info.Datatype = 'double';
        niftiwrite(mk, fullfile(outdir,'mk.nii'), info);
        niftiwrite(rk, fullfile(outdir,'rk.nii'), info);
        niftiwrite(ak, fullfile(outdir,'ak.nii'), info);
    end

    if wmti
        disp('...getting WMTI params')
        info.ImageSize = info.ImageSize(1:3);
        info.PixelDimensions = info.PixelDimensions(1:3);
        info.Datatype = 'double';
        info.DisplayIntensityRange = [0 1];
        [awf, eas, ias] = wmti_parameters(dt, mask);
        disp('...saving WMTI params')
        niftiwrite(awf, fullfile(outdir,'awf.nii'), info);
        fields = fieldnames(ias);
        for ii=1:numel(fields)
            paramsii = getfield(ias, fields{ii});
            savename = fullfile(outdir, ['ias_', fields{ii}, '.nii']);
            info.DisplayIntensityRange = [0 max(paramsii(:))];
            niftiwrite(paramsii, savename, info);
        end
        fields = fieldnames(eas);
        for ii=1:numel(fields)
            paramsii = getfield(eas, fields{ii});
            savename = fullfile(outdir, ['eas_', fields{ii}, '.nii']);
            info.DisplayIntensityRange = [0 max(paramsii(:))];
            niftiwrite(paramsii, savename, info);        
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

end

