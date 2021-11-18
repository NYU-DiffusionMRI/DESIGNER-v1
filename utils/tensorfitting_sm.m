function tensorfitting(dwi,root,outdir,detectoutliers,cumulants,dti,dki,wmti,fitconstraints,akc,DKIroot,fitWDKI,smooth)

    addpath(genpath(DKIroot));
    
%     if isempty(gcp('nocreate'))
%         pc = parcluster('local');
%         pc.JobStorageLocation = '/cbi05data/data1/Hamster/scratch';
%         nw = 12;
%         parpool(pc,nw);
%     end

    maskex = exist(fullfile(root,'mask.nii.gz'),'file');
    if maskex == 2
        nii = niftiread(fullfile(root,'mask.nii.gz')); mask = logical(nii);
    end
    %nii = niftiread(fullfile(root,'dwi_designer.nii')); dwi = double(nii);
    info = niftiinfo(fullfile(root,'data.nii'));
    ndwis = size(dwi,4); %pixdim5 = nii.hdr.dime.pixdim(5);
    if maskex == 0
        mask = logical(ones(size(dwi,1),size(dwi,2),size(dwi,3)));
    end

    bvallist = dir(fullfile(root,'data.bval')); bvaldir = bvallist.name;
    bveclist = dir(fullfile(root,'data.bvec')); bvecdir = bveclist.name;
    bval = textread(fullfile(root,bvaldir)); bval = bval(:, 1:ndwis)'; bval = bval./1000;
    bvec = textread(fullfile(root,bvecdir)); bvec = bvec(:, 1:ndwis)';
    maxbval = 3;
    

    detectoutliers = logical(detectoutliers);
    dti = logical(dti);
    dki = logical(dki);
    wmti = logical(wmti);
    akc = logical(akc);
    fitWDKI = logical(fitWDKI);
    if ischar(fitconstraints)
        conts = split(fitconstraints,',');
        constraints(1) = str2double(conts(1));
        constraints(2) = str2double(conts(2));
        constraints(3) = str2double(conts(3));
    else
        constraints = [0,1,0];
    end
    cumulants = logical(cumulants);
    
%     mesoset1 = dwi(:,:,:,1:23);
%     mesoset2 = dwi(:,:,:,24:end);
%     b01 = mesoset1(:,:,:,1);
%     b02 = mesoset2(:,:,:,1);
%     ratio = b01./b02;
%     %mesoset1 = mesoset1./b01;%.*ratio;
%     mesoset2 = mesoset2.*ratio;
%     dwi = cat(4,mesoset1,mesoset2);
    

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
         [dt,dwi,akc_out] = correctDt(dt,dwi,mask,bval,bvec);
         disp(['N outliers = ',num2str(sum(akc_out(:)))]);
        
        info.ImageSize = info.ImageSize(1:3);
        info.PixelDimensions = info.PixelDimensions(1:3);
        info.Datatype = 'single';
        niftiwrite(akc_out, fullfile(outdir,'akc_out1.nii'), info);
        
        
         akc_out = outlierdetection(dt,mask);
         akc_out(isnan(akc_out)) = 0;
         disp(['N outliers = ',num2str(sum(akc_out(:)))]);
        info.ImageSize = info.ImageSize(1:3);
        info.PixelDimensions = info.PixelDimensions(1:3);
        niftiwrite(akc_out, fullfile(outdir,'akc_out2.nii'), info);
    end
    
    if smooth
        disp('...smoothing');
        dwi = nlmsmooth(dwi,mask,akc_out);
        [b0,dt] = dki_fit(dwi,[bvec,bval],mask,constraints,[],maxbval);
        se = strel('cube',3);
        [dt,dwi,akc_out] = correctDt(dt,dwi,imerode(mask,se),bval,bvec);
        disp(['N outliers = ',num2str(sum(akc_out(:)))]);
    end
        

    if cumulants
        dt = vectorize(dt, mask);
        D_apprSq = (sum(dt([1 4 6],:),1)/3).^2;
        dt(7:21,:) = dt(7:21,:) .* D_apprSq(ones(15,1),:);
        dt = vectorize(dt, mask);
    end

    %disp('...getting DTI and DKI params')
    
    if dti
        disp('...getting DTI params')
        list = round(bval) == 0 | round(bval) == 1;
        dwi_dti = dwi(:,:,:,list);
        bvec_dti = bvec(list,:);
        bval_dti = bval(list);
        [b0_dti, dt_dti, md_dti, rd_dti, ad_dti, fa_dti] = dti_fit(dwi_dti,[bvec_dti,bval_dti],mask);
        
        disp('...saving DTI params')
        info.ImageSize = info.ImageSize(1:3);
        info.PixelDimensions = info.PixelDimensions(1:3);
        info.Datatype = 'double';
        info.DisplayIntensityRange = [0 1];
        niftiwrite(fa_dti, fullfile(outdir,'fa_dti.nii'), info);
        info.DisplayIntensityRange = [0 3];
        niftiwrite(md_dti, fullfile(outdir,'md_dti.nii'), info);
        niftiwrite(rd_dti, fullfile(outdir,'rd_dti.nii'), info);
        niftiwrite(ad_dti, fullfile(outdir,'ad_dti.nii'), info);
        info.DisplayIntensityRange = [0 0];
        niftiwrite(b0_dti, fullfile(outdir,'b0_dti.nii'), info);
        info.ImageSize = cat(2, info.ImageSize, 3);
        info.PixelDimensions = cat(2, info.PixelDimensions, NaN);
        info.ImageSize(4) = 6;
        niftiwrite(dt_dti, fullfile(outdir,'dt.nii'), info);
    end
    if dki
        disp('...getting DKI params')
        [fa, md, rd, ad, fe, mk, rk, ak, l1, l2, l3] = dki_parameters(dt,mask);
        fe = cat(4,fa.*fe(:,:,:,1),fa.*fe(:,:,:,2),fa.*fe(:,:,:,3));
        disp('...saving DKI params')
        info.ImageSize = info.ImageSize(1:3);
        info.PixelDimensions = info.PixelDimensions(1:3);
        info.DisplayIntensityRange = [0 3];
        info.Datatype = 'double';
        niftiwrite(mk, fullfile(outdir,'mk.nii'), info);
        niftiwrite(rk, fullfile(outdir,'rk.nii'), info);
        niftiwrite(ak, fullfile(outdir,'ak.nii'), info);
        info.ImageSize = info.ImageSize(1:3);
        info.PixelDimensions = info.PixelDimensions(1:3);
        info.DisplayIntensityRange = [0 1];
        niftiwrite(fa, fullfile(outdir,'fa.nii'), info);
        info.DisplayIntensityRange = [0 3];
        niftiwrite(md, fullfile(outdir,'md.nii'), info);
        niftiwrite(rd, fullfile(outdir,'rd.nii'), info);
        niftiwrite(ad, fullfile(outdir,'ad.nii'), info);
        niftiwrite(l1, fullfile(outdir,'l1.nii'), info);
        niftiwrite(l2, fullfile(outdir,'l2.nii'), info);
        niftiwrite(l3, fullfile(outdir,'l3.nii'), info);
        info.DisplayIntensityRange = [0 0];
        niftiwrite(b0, fullfile(outdir,'b0.nii'), info);
        info.ImageSize = cat(2, info.ImageSize, 3);
        info.PixelDimensions = cat(2, info.PixelDimensions, NaN);
        niftiwrite(fe, fullfile(outdir,'fe.nii'), info);
        info.ImageSize(4) = 21;
        niftiwrite(dt, fullfile(outdir,'dkt.nii'), info);
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
    
    if fitWDKI
        disp('...getting W params')
        [fa, md, rd, ad, fe, mw,  rw, aw] = dwi_parameters(dt,mask);
        fe = cat(4,fa.*fe(:,:,:,1),fa.*fe(:,:,:,2),fa.*fe(:,:,:,3));
        disp('...saving DWI params')
        info.ImageSize = info.ImageSize(1:3);
        info.PixelDimensions = info.PixelDimensions(1:3);
        info.DisplayIntensityRange = [0 3];
        info.Datatype = 'double';
        niftiwrite(mw, fullfile(outdir,'mw.nii'), info);
        niftiwrite(rw, fullfile(outdir,'rw.nii'), info);
        niftiwrite(aw, fullfile(outdir,'aw.nii'), info);
        info.ImageSize = info.ImageSize(1:3);
        info.PixelDimensions = info.PixelDimensions(1:3);
        info.DisplayIntensityRange = [0 1];
        niftiwrite(fa, fullfile(outdir,'fa_dwi.nii'), info);
        info.DisplayIntensityRange = [0 3];
        niftiwrite(md, fullfile(outdir,'md_dwi.nii'), info);
        niftiwrite(rd, fullfile(outdir,'rd_dwi.nii'), info);
        niftiwrite(ad, fullfile(outdir,'ad_dwi.nii'), info);

        info.DisplayIntensityRange = [0 0];
        info.ImageSize = cat(2, info.ImageSize, 3);
        info.PixelDimensions = cat(2, info.PixelDimensions, NaN);
        niftiwrite(fe, fullfile(outdir,'fe_dwi.nii'), info);     
    end
    
    function dwi = nlmsmooth(dwi,mask, akc)
        kernel = 5;
        k = floor(kernel/2);
        k_ = ceil(kernel/2);
        maskinds = find(mask);
        dwi_ = zeros(length(maskinds),size(dwi,4));
        dwi_norm = zeros(size(dwi));
        for i = 1:size(dwi,4)
            dwii = abs(dwi(:,:,:,i));
            dwi_norm(:,:,:,i) = dwii./max(dwii(:));
        end
        parfor index = 1 : length(maskinds)
            thisLinearIndex = maskinds(index);
            % Get the x,y,z location
            [x,y,z] = ind2sub(size(mask), thisLinearIndex);
            if (z < k_ || z > size(dwi,3)-k) || (x < k_ || x > size(dwi,1)-k) || (y < k_ || y > size(dwi,1)-k)
                continue
            end
    
            akcpatch = reshape(logical(akc(x-k:x+k,y-k:y+k,z-k:z+k)),[kernel^3,1]);
            ref = repmat(reshape(dwi_norm(x,y,z,:),[1,size(dwi,4)]),[kernel^3,1]);
            patch = reshape(dwi_norm(x-k:x+k,y-k:y+k,z-k:z+k,:),[kernel^3,size(dwi,4)]);
            patchorig = reshape(dwi(x-k:x+k,y-k:y+k,z-k:z+k,:),[kernel^3,size(dwi,4)]);
            intensities = sqrt(sum((patch-ref).^2,2))./size(dwi,4);
            [min_wgs,min_idx] = sort(intensities, 'ascend');
            
            wgs_max = min_wgs(end);
            min_wgs(akcpatch) = wgs_max;
            
            goodidx = min_wgs < prctile(min_wgs,10);
            min_idx = min_idx(goodidx);
            min_wgs = min_wgs(goodidx);
            wgs_max = max(min_wgs);
            
            wgs_inv = wgs_max - min_wgs;
            wgs_nrm = wgs_inv/sum(wgs_inv);
            wval = sum(patchorig(min_idx,:).*(wgs_nrm*ones(1,size(dwi,4))));
            dwi_(index,:)= wval;
         end
         for idx = 1:length(maskinds)
            thisLinearIndex = maskinds(idx);
            [x,y,z] = ind2sub(size(mask),thisLinearIndex);
            dwi(x,y,z,:) = dwi_(idx,:);
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

    function[dt,dwi,akc_out] = correctDt(dt,dwi,mask,bval,bvec)
        akc_out = outlierdetection(dt,mask);
        akc_out(isnan(akc_out)) = 0;
        
        [dwi_,nanlinearinds] = naninds(akc_out,dwi);
        for idx = 1:length(nanlinearinds)
            thisLinearIndex = nanlinearinds(idx);
            [x,y,z] = ind2sub(size(akc_out),thisLinearIndex);
            dwi(x,y,z,:) = dwi_(idx,:);
        end
        
        [~,dt_] = dki_fit(dwi,[bvec,bval],logical(akc_out),constraints,[],3);

        for d = 1:size(dt,4)
            dtf_ = dt_(:,:,:,d); 
            dto_ = dt(:,:,:,d);
            dto_(logical(akc_out)) = dtf_(logical(akc_out));
            dt(:,:,:,d) = dto_;
        end
    end
  
    
    function [dwi_,nanLinearIndexes] = naninds(nanlocations,dwi)
        nanLinearIndexes = find(nanlocations);
   
        % Get the x,y,z of all other locations that are non nan and non csf.
        kernel = 5;
        k = floor(kernel/2);
        k_ = ceil(kernel/2);
    
        dwi_ = zeros(length(nanLinearIndexes),size(dwi,4));
        dwi_norm = zeros(size(dwi));
        for i = 1:size(dwi,4)
            dwii = abs(dwi(:,:,:,i));
            dwi_norm(:,:,:,i) = dwii./max(dwii(:));
        end
        parfor index = 1 : length(nanLinearIndexes)
            thisLinearIndex = nanLinearIndexes(index);
            % Get the x,y,z location
            [x,y,z] = ind2sub(size(nanlocations), thisLinearIndex);
            if (z < k_ || z > size(dwi,3)-k) || (x < k_ || x > size(dwi,1)-k) || (y < k_ || y > size(dwi,1)-k)
                continue
            end
            %keyboard
            akcpatch = reshape(logical(nanlocations(x-k:x+k,y-k:y+k,z-k:z+k)),[kernel^3,1]);
            
            ref = repmat(reshape(dwi_norm(x,y,z,:),[1,size(dwi,4)]),[kernel^3,1]);
            patch = reshape(dwi_norm(x-k:x+k,y-k:y+k,z-k:z+k,:),[kernel^3,size(dwi,4)]);
            patchorig = reshape(dwi(x-k:x+k,y-k:y+k,z-k:z+k,:),[kernel^3,size(dwi,4)]);
            intensities = sqrt(sum((patch-ref).^2,2)./size(dwi,4));
            [min_wgs,min_idx] = sort(intensities, 'ascend');
            
            wgs_max = min_wgs(end);
            min_wgs(akcpatch) = wgs_max;
            
            goodidx = min_wgs < median(min_wgs);
            min_idx = min_idx(goodidx);
            min_wgs = min_wgs(goodidx);
            wgs_max = max(min_wgs);
            
            wgs_inv = wgs_max - min_wgs;
            wgs_nrm = wgs_inv/sum(wgs_inv);
            wval = sum(patchorig(min_idx,:).*(wgs_nrm*ones(1,size(dwi,4))));
            dwi_(index,:)= wval;
        end
                   
    end



end

