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
[fa, md, rd, ad, fe, mk, rk, ak, kfa, mkt, medianFilter] = dki_parameters(dt,mask,violMask,medianfilter);
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
    [awf, eas, ias] = wmti_parameters(dt, mask, violMask, medianFilter);
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

%%   Split and save dt into DT and KT for MUSC DKE_FT
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

%%   Save violation mask
if ~exist(fullfile(outdir,'QC'))
    mkdir(fullfile(outdir,'QC'));
else
    ;
end
nii.img = violMask.Proportional; nii.hdr.dime.glmax = max(b0(:)); save_untouch_nii(nii,fullfile(outdir,'QC','Propotional_Violations.nii'));
nii.img = violMask.Directional; nii.hdr.dime.glmax = max(b0(:)); save_untouch_nii(nii,fullfile(outdir,'QC','Good_Directions.nii'));

%% Create SNR Plots
mkdir(fullfile(outdir,'QC'));

%   Load all files
%   Run SNR plotting only if user specific MP-PCA option
noisePath = dir(fullfile(outdir,'fullnoisemap.nii'));
if size(noisePath,1) > 0
    disp('...generating SNR plot');
    noisePath = fullfile(noisePath.folder,noisePath.name);
    Inoise = niftiread(noisePath);
    Inoise(find(isnan(Inoise))) = 0;
    
    bvalPath = dir(fullfile(outdir,'*bval'));
    bvalPath = fullfile(bvalPath.folder,bvalPath.name);
    bval = round(load(bvalPath));
    
    procPath = dir(fullfile(outdir,'dwi_designer.nii'));
    procPath = fullfile(procPath.folder,procPath.name);
    Iproc = niftiread(procPath);
    
    rawPath = dir(fullfile(outdir,'dwi_raw.nii'));
    rawPath = fullfile(rawPath.folder,rawPath.name);
    Iraw = niftiread(rawPath);
    
    brainPath = dir(fullfile(outdir,'brain_mask.nii'));
    brainPath = fullfile(brainPath.folder,brainPath.name);
    brainMask = logical(niftiread(brainPath));
    
    %   Create logical indexes of bvalues
    listBval = unique(bval);
    bIdx = zeros(length(listBval),length(bval));
    for j = 1:length(listBval)
        bIdx(j,:) = bval == listBval(j);
    end
    bIdx = logical(bIdx);
    
    %   Form index vector for mean 4D array
    nB0 = find(bIdx(1,:));
    meanIdx = horzcat(repmat(listBval(1),[1 numel(nB0)]),listBval(2:end));
    
    for j = 1:length(listBval)
        for k = 1:length(meanIdx)
            if k <= length(nB0);
                tmpProc = Iproc(:,:,:,nB0(k)) ./ Inoise;
                snrProc(:,k) = tmpProc(brainMask);
                
                tmpRaw = Iraw(:,:,:,nB0(k)) ./ Inoise;
                snrRaw(:,k) = tmpRaw(brainMask);
                
            elseif meanIdx(k) == listBval(j)
                tmpProc = mean(Iproc(:,:,:,bIdx(j,:)),4)./Inoise;
                snrProc(:,k) = tmpProc(brainMask);
                
                tmpRaw = mean(Iraw(:,:,:,bIdx(j,:)),4)./Inoise;
                snrRaw(:,k) = tmpRaw(brainMask);
            else
                continue;
            end
        end
    end
    
    %   Histogram counts --------------------------------------------------
    %   Generate subplot layout such that plots are in nPlots x 2 layout.
    %   The firsy column contains plots and the second contains a single
    %   legend. Plot layout sets the size of subplot and plotLocations
    %   determines where to plot them on a nPlots x 2 grid.
    
    nbins = 100;    % No. of bins for SNR plot
    plotLayout = [numel(listBval) 3];
    plotLocation = linspace(1,(numel(listBval)+4),numel(listBval));
    
    figure; fig = gcf;
    
    set(fig,'PaperUnits','inches','PaperPosition',[.25 .25 10 12],...
        'InvertHardcopy','off','Color','white','Visible','off');
    
    for j = 1:length(listBval)
        [Nproc Eproc] = histcounts(snrProc(:,meanIdx == listBval(j)),...
            nbins,'Normalization','Probability');
        [Nraw Eraw] = histcounts(snrRaw(:,meanIdx == listBval(j)),...
            nbins,'Normalization','Probability');
        
        Nproc = smooth(Nproc);  % Smooth vector for easy visualization
        Nraw = smooth(Nraw);
        
        %   Compute median value of each bin
        for k = 1:length(Nproc)
            Mproc(k) = median([Eproc(k) Eproc(k+1)]);
        end
        for k = 1:length(Nraw)
            Mraw(k) = median([Eraw(k) Eraw(k+1)]);
        end
        
        %   Plot graphs ---------------------------------------------------
        
        %   Plotting colors
        c1 = [86,187,131]/255;   % processed color
        c2 = [78,173,241]/255;   % raw color
        c3 = [235,235,235]/255;  % background color
        
        subplot(numel(listBval),3,[plotLocation(j) plotLocation(j)+1])
        hold on;
        pArea = area(Mproc,Nproc,...
            'EdgeColor',c1,'LineWidth',3,...
            'FaceColor',c1,'FaceAlpha',0.55);
        rArea = area(Mraw,Nraw,...
            'EdgeColor',c2,'LineWidth',3,...
            'FaceColor',c2,'FaceAlpha',0.55);
        hold off;
        xlabel('SNR');
        ylabel('%Voxels');
        
        if listBval(j) == 0
            xlim([0 500]);
        elseif listBval(j) == 1000
            xlim([0 80]);
        elseif listBval(j) == 2000
            xlim([1 40]);
        else
            xlim([1 20]);
        end
        
        title(sprintf('B%d',listBval(j)));
        grid on; box on; legend off;
        %   Set axial colors
        ax = gca;
        set(ax,'Color',c3,...
            'GridColor','white','GridAlpha',1,'MinorGridAlpha',0.15,...
            'fontname','helvetica','FontWeight','bold','fontsize',14);
        clear Nproc Eproc Nraw Eraw Mproc Mraw
    end
    legPlot = subplot(numel(listBval),3,numel(listBval)*3);
    posLeg = get(legPlot,'Position');
    leg1 = legend([pArea,rArea],...
        {'After Preprocessing','Before Preprocessing'});
    set(leg1,'Position',posLeg);
    title(leg1,'Legend');
    axis(legPlot,'off');
    print(fullfile(outdir,'QC','SNR_Plots'),'-dpng','-r800');
    clear tmpRaw tmpProc snrProc snrRaw
else
    disp('...MP-PCA disabled, skipping SNR plot');
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

