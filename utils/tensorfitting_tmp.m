function tensorfitting_tmp(root,outdir,detectoutliers,cumulants,dti,dki,wmti,fitconstraints,akc,DKIroot,fitWDKI)

    addpath(genpath(DKIroot));
    
%     if isempty(gcp('nocreate'))
%         pc = parcluster('local');
%         pc.JobStorageLocation = '/cbi05data/data1/Hamster/scratch';
%         nw = 12;
%         parpool(pc,nw);
%     end

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
        csfmask = niftiread(fullfile(root,'tissue_seg.nii.gz'));
         disp('...running AKC correction')
        [dt,akc_out] = shCorrectDt(dt,dwi,mask,bval,bvec,5,csfmask);
        
         info.ImageSize = info.ImageSize(1:3);
        info.PixelDimensions = info.PixelDimensions(1:3);
        niftiwrite(akc_out, fullfile(root,'akc_out1.nii'), info);
        
        disp(['N outliers = ',num2str(sum(akc_out(:)))]);
         akc_out = outlierdetection(dt,mask);
         akc_out(isnan(akc_out)) = 0;
         disp(['N outliers = ',num2str(sum(akc_out(:)))]);
         
%          [dt,akc_out] = shCorrectDt(dt,dwi,mask,bval,bvec,5);
%          disp(['N outliers = ',num2str(sum(akc_out(:)))]);
%         [dt,akc_out] = shCorrectDt(dt,dwi,mask,bval,bvec,5);
%         disp(['N outliers = ',num2str(sum(akc_out(:)))]);
%         akc_out = outlierdetection(dt,mask);
%         akc_out(isnan(akc_out)) = 0;
%         disp(['N outliers = ',num2str(sum(akc_out(:)))]);

        
%         for v = 1:size(dt,4)
%             dt_v = dt(:,:,:,v);
%             dt_v(logical(akc_out)) = NaN;
%             if verLessThan('matlab','9.2')
%                 dt_f = fillmissing(dt_v,'linear');
%             else
%                 dt_f = fillmissing(dt_v,'movmedian',5);
%             end
%             dt(:,:,:,v) = dt_f;
%         end
        
        info.ImageSize = info.ImageSize(1:3);
        info.PixelDimensions = info.PixelDimensions(1:3);
        niftiwrite(akc_out, fullfile(root,'akc_out2.nii'), info);
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
        disp('...getting axially symmetric W params')
        woutdir = fullfile(outdir,'Wparams');
        mkdir(woutdir)
        info.ImageSize = info.ImageSize(1:3);
        info.PixelDimensions = info.PixelDimensions(1:3);
        info.Datatype = 'double';
        [b0,Dlm,Wlm,Dl,Wl,DTI_scalars,DKI_scalars,ExtraScalars] = Fit_DKI_LTE(dwi,bval,bvec, mask, 1);
        info.DisplayIntensityRange = [0 1];
        niftiwrite(DTI_scalars.fa, fullfile(woutdir,'fa.nii'), info);
        info.DisplayIntensityRange = [0 3];
        niftiwrite(DTI_scalars.md, fullfile(woutdir,'md.nii'), info);
        niftiwrite(DTI_scalars.ad_axsym, fullfile(woutdir,'ad_axsm.nii'), info);
        niftiwrite(DTI_scalars.rd_axsym, fullfile(woutdir,'rd_axsm.nii'), info);
        niftiwrite(DKI_scalars.mw, fullfile(woutdir,'mw.nii'), info);
        niftiwrite(DKI_scalars.aw_axsym, fullfile(woutdir,'aw_axsm.nii'), info);
        niftiwrite(DKI_scalars.rw_axsym, fullfile(woutdir,'rw_axsm.nii'), info);
        niftiwrite(ExtraScalars.rd, fullfile(woutdir,'rd_exact.nii'), info);
        niftiwrite(ExtraScalars.ad, fullfile(woutdir,'ad_exact.nii'), info);
        niftiwrite(ExtraScalars.rd, fullfile(woutdir,'rd_exact.nii'), info);
        niftiwrite(ExtraScalars.mk, fullfile(woutdir,'mk_exact.nii'), info);
        niftiwrite(ExtraScalars.aw, fullfile(woutdir,'aw_exact.nii'), info);
        niftiwrite(ExtraScalars.rw, fullfile(woutdir,'rw_exact.nii'), info);
        niftiwrite(Wl(:,:,:,1), fullfile(woutdir,'W0.nii'), info);
        niftiwrite(Wl(:,:,:,2), fullfile(woutdir,'W1.nii'), info);
        niftiwrite(Wl(:,:,:,3), fullfile(woutdir,'W2.nii'), info);
        niftiwrite(Dl(:,:,:,1), fullfile(woutdir,'D0.nii'), info);
        niftiwrite(Dl(:,:,:,2), fullfile(woutdir,'D2.nii'), info);        
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

    function[dt,akc_out] = shCorrectDt(dt,dwi,mask,bval,bvec,lmax,csfmask)
        akc_out = outlierdetection(dt,mask);
        akc_out(isnan(akc_out)) = 0;
        [closeInds, ginds,nani] = naninds(akc_out,lmax,csfmask,dwi(:,:,:,1));
        
        dwitmp = dwi;
        for i = 1:size(dwi,4)
            dwitmp_ = dwitmp(:,:,:,i);
            dwitmp_(logical(akc_out)) = NaN;
            
            dwifilled_ = fillnans(dwitmp_,closeInds,ginds,nani,csfmask);
            %dwifilled_ = inpaint_nans3(dwitimp_);
            
%             dwifilled_ = fillmissing(dwitmp_,'movmean',[lmax,lmax],1);
%             dwifilled_ = fillmissing(dwifilled_,'movmean',[lmax,lmax],3);
%             dwifilled_ = fillmissing(dwifilled_,'movmean',[lmax,lmax],2);
%             
%             nanremain = isnan(dwifilled_);
%             if sum(nanremain(:)) > 0
%             	dwifilled_ = fillmissing(dwitmp_,'movmean',[lmax+2,lmax+2],1);
%                 dwifilled_ = fillmissing(dwifilled_,'movmean',[lmax+2,lmax+2],3);
%                 dwifilled_ = fillmissing(dwifilled_,'movmean',[lmax+2,lmax+2],2);
%             end
                
%             lmaxi = lmax;
%             while sum(nanremain(:)) > 0
%                 sum(nanremain(:))
%                 lmaxi = lmaxi + 1
%                 dwifilled_ = fillmissing(dwifilled_,'movmean',lmaxi);
%             end
            %dwifilled_ = fillnans(dwitmp_,akc_out,lmax);
            dwi(:,:,:,i) = dwifilled_;
        end
        [~,dt_] = dki_fit(dwi,[bvec,bval],logical(akc_out),[0,0,0],[],3);
        
%         bu = unique(bval);
%         dwitmp = [];
%         bvaltmp = [];
%         bvectmp = [];
%         for bi = 1:numel(bu)
%             if sum(bval==bu(bi))<=4
%                 lmax = 0;
%             end
%             inds = find(bval==bu(bi));
%             dwib_ = vectorize(dwi(:,:,:,inds), logical(akc_out));
%             shb = SHfit(dwib_,bvec(inds,:),lmax);
%             dwip_ = SHproj(shb,bvec(inds,:),lmax);
% %            dwip = vectorize(dwip_, logical(akc_out));
%             
%             pcthr = 0.2;
%             pc = mean(abs(dwip_-dwib_)./dwib_,2);
%             lowpcinds = inds(pc < pcthr);
%             dwitmp = cat(4,dwitmp,dwi(:,:,:,lowpcinds));
%             bvaltmp = cat(1,bvaltmp,bval(lowpcinds));
%             bvectmp = cat(1,bvectmp,bvec(lowpcinds,:));
%        
% 
% %             figure;
% %             subplot(1,2,1)
% %             plot(1:length(inds),mean(dwib_,2),'bo',1:length(inds),mean(dwip_,2),'ro')
% %             subplot(1,2,2)
% %             plot(1:length(inds),mean(abs(dwip_-dwib_)./dwib_,2),'ko')
% %             keyboard
%             
% %             for bj = 1:numel(inds)
% %                 dwip_ = dwip(:,:,:,bj);
% %                 dwi_ = dwi(:,:,:,inds(bj));
% %                 dwi_(logical(akc_out)) = dwip_(logical(akc_out));
% %                 dwi(:,:,:,inds(bj)) = dwi_;
% %             end
%         end
%         keyboard
%         [~,dt_] = dki_fit(dwitmp,[bvectmp,bvaltmp],logical(akc_out),[0,0,0],[],3);
        for d = 1:size(dt,4)
            dtf_ = dt_(:,:,:,d); 
            dto_ = dt(:,:,:,d);
            dto_(logical(akc_out)) = dtf_(logical(akc_out));
            dt(:,:,:,d) = dto_;
        end
    end
   
    function u = fillnans(u, nanlocations,ginds,nanLinearIndexes,seg)
        seg = seg+1;
        for index = 1 : size(nanlocations,1)
            thisLinearIndex = nanLinearIndexes(index);
            [x,y,z] = ind2sub(size(u), thisLinearIndex);
            ttype = seg(x,y,z);
            
            goodValue = 0;
            for ind = 1:size(nanlocations,2)
                idx = nanlocations(index,ind);
                goodValue = goodValue + u(ginds{ttype}(1,idx), ginds{ttype}(2,idx), ginds{ttype}(3,idx));
            end
            % Replace the bad nan value in u with the good value.
            gv = goodValue/size(nanlocations,2);
            u(x,y,z) = gv;
        end
    end

    function [indexOfClosest,ginds,nanLinearIndexes] = naninds(nanlocations,nvals,seg,b0)
        nanLinearIndexes = find(nanlocations);
        seg = seg + 1;
%         csfLinearIndices = find(csfmask);
        ginds = cell(1,max(seg(:)));
        for m = 1:max(seg(:))
            ttypeLinearIndexes = find(seg~=m);
            nonNanLinearIndexes = setdiff(1:numel(nanlocations), union(nanLinearIndexes, ttypeLinearIndexes));
            [xGood, yGood, zGood] = ind2sub(size(nanlocations), nonNanLinearIndexes);
            ginds{m} = [xGood; yGood; zGood];
        end
%         nonNanLinearIndexes = setdiff(1:numel(nanlocations), union(nanLinearIndexes,csfLinearIndices));
%         % Get the x,y,z of all other locations that are non nan and non csf.
%         [xGood, yGood, zGood] = ind2sub(size(nanlocations), nonNanLinearIndexes);
%         ginds = [xGood; yGood; zGood];

        indexOfClosest = zeros(length(nanLinearIndexes),nvals);
        for index = 1 : length(nanLinearIndexes)
            thisLinearIndex = nanLinearIndexes(index);
            % Get the x,y,z location
            [x,y,z] = ind2sub(size(nanlocations), thisLinearIndex);
            ttype = seg(x,y,z);
            % Get distances of this location to all the other locations
            distances = sqrt((x-ginds{ttype}(1,:)).^2 + (y - ginds{ttype}(2,:)).^ 2 + (z - ginds{ttype}(3,:)).^ 2);
            intensities = sqrt((b0(x,y,z) - b0(seg==ttype)).^2);
            [~,sortedIndexes] = mink(intensities.*distances', nvals);
            keyboard
            %[~, sortedIndexes] = sort(intensities.*distances, 'ascend');
            % The closest non-nan value will be located at index sortedIndexes(1)
            indexOfClosest(index,:) = sortedIndexes;
        end
    end

end

