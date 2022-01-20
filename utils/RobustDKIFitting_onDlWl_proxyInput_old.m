function [dt, b0, list] = RobustDKIFitting_onDlWl_proxyInput(dwi, grad, mask, DlWlproxy,alpha_tuning)
%
% THIS CAN BE MADE FASTER, JUST NEED TO WRITE THE DERIVATIVES OF THE OBJ
% FUNCTION
%
%
%
    % [dt, b0, DlWlpred, list] = RobustDKIFitting_onWl(dwi, grad, mask,alpha_tuning)
    %
    % This is a modification of Jelle's RobustDKIFitting
    % Here the W0-W2-W4 of W tensor are regularized
    % For more details, contact: Santiago.Coelho@nyulangone.org 
    %
    %
    % RobustDKIFitting provides a novel strategy "Voxel Quality Transfer" for the robust estimation of the diffusion kurtosis tensor.
    % The technique is build upon the fact that (1) the powder kurtosis is a good proxy of the mean kurtosis and 
    % (2) the powder kurtosis is more robust metric than the mean kurtosis. Note that the powder kurtosis is the AKC, 
    % estimated from powder-averaged diffusion-weighted MRI data 
    % 
    % -----------------------------------------------------------------------------------
    % please cite: https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28730
    % ------------------------------------------------------------------------------------
    %     
    % Although we demonstrated a significant improvement in the precision
    % and robustness in the estimation of the DKI model parameters, the
    % proposed technique might benefit from furhter development, including
    % optimization in terms of computation speed.
    %
    % Note that this first code release does not allow for parallel
    % processing and the rather slow NLS estimator is used. 
    %
    % Feel free to contact Rafael Henriques or Jelle Veraart if you'd like
    % to contribute to this project. 
    %
    % Usage:
    % ------
    % [b0, dt] = RobustDKIFitting(dwi, grad, mask)
    %
    % Required input: 
    % ---------------
    %     1. dwi: diffusion-weighted images.
    %           [x, y, z, ndwis]
    %       
    %       Important: We recommend that you apply denoising, gibbs correction, motion-
    %       and eddy current correction to the diffusion-weighted image
    %       prior to tensor fitting. Thes steps are not includede in this
    %       tools, but we are happy to assist (Jelle.Veraart@nyumc.org).
    %
    %     2. grad: diffusion encoding information (gradient direction 'g = [gx, gy, gz]' and b-values 'b')
    %           [ndwis, 4] 
    %           format: [gx, gy, gx, b]
    %    
    %       IMPORTANT: "shell" your b-values before running this
    %       code.
    %
    %     3. mask (boolean; [x, y, x]), providing a mask limits the
    %       calculation to a user-defined region-of-interest.
    %    
    % Output: 
    % -------
    %    1. dt, diffusion kurtosis tensor. The tensor elements are stored
    %    in the following way 
    %         [Dxx Dxy Dxz Dyy Dyx Dzz Wxxxx Wxxxy Wxxxz Wxxyy  Wxxyz Wxxzz Wxyyy Wxyyz Wxyzz Wxzzz Wyyyy  Wyyyz Wyyzz Wyzzz Wyzzz]
    %    2. b0, estimated non-diffusion weighted image
    %    3. mkpred, predicted mean kurtosis using Voxel Quality Transfer.
    %               This is not the estimated mean kurtosis.
    %
    %
    % Copyright (c) 2021 New York University 
    %
    % This Source Code Form is subject to the terms of the Mozilla Public
    % License, v. 2.0. If a copy of the MPL was not distributed with this file,
    % You can obtain one at http://mozilla.org/MPL/2.0/
    % 
    % This code is distributed  WITHOUT ANY WARRANTY; without even the 
    % implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    % 
    % For more details, contact: Jelle.Veraart@nyulangone.org 
    
    sigma_sm=0.65;
    sm_pz=3;
    extraRobustness = true;

    
    %% Step 1: prepare data
    % We (a) do some quick input checks, (b) vectorize the data to avoid handling with 4D data
    % structures and nested loops, and (c) build the DKI-specific b-matrix.  
    
    dwi = double(dwi);
    dwi(dwi<=0)=eps;
    scaleFact = 1000/max(dwi(:));
    dwi = dwi*scaleFact;
    
    if ~exist('alpha_tuning','var') || isempty(alpha_tuning)
        alpha_tuning = [1 1/3 5 2 1/2];
    end
    
    [x, y, z, ndwis] = size(dwi);
    if ~exist('grad','var') || size(grad,1) ~= ndwis || size(grad,2) ~= 4
        error('');
    end
    grad = double(grad);
    grad(:,1:3) = bsxfun(@rdivide,grad(:,1:3),sqrt(sum(grad(:,1:3).^2,2))); grad(isnan(grad)) = 0;
    bval = grad(:, 4);
    if ~exist('mask','var') || isempty(mask)
        mask = true(x, y, z);
    end
   
    y = vectorize(dwi, mask);
    bmat = bmatrix(grad);
    [~, nvox] = size(y);  
   
    %% Step 2: Regular NLS fit
    % We perform a regular  (unconstrained) NLS fit. The fit is initiated
    % with a LLS fit. This estimate will be the final result for all "good"
    % voxels.
    
    start = bmat\log(y);         %
    options = optimoptions('lsqnonlin','Jacobian','on','TolFun',1e-8,'TolX', ...
                        1e-8,'MaxIter',10000,'Display','off');
    
     % this is still very slow. It would be possible to switch everything
     % to WLLS to reduce the computation time. 
     start(~isfinite(start))=eps;
     nls=0*start;
     parfor ii = 1:nvox
         data=double(y(:, ii));
         if any(~isfinite(data))
             continue
         end
         nls(:,ii) = lsqnonlin(@(x)ObjF(x,double(bmat),data),double(start(:,ii)),[],[],options);
     end
     clear start

     y_hat = exp(bmat*nls);
     dt.nls = W2K(nls(2:22, :));     
     
    %% Step 2.5:
    % Repeat nls fit on smooth signals
    dwi_sm=0*dwi;
    for ii=1:ndwis
        dwi_sm(:,:,:,ii)=smooth3(dwi(:,:,:,ii),'gaussian',sm_pz*[1 1 1],sigma_sm);
    end
%     IMGUI(cat(4,dwi(:,:,:,24),dwi_sm(:,:,:,24)))
    y_sm = vectorize(dwi_sm, mask);
    start_sm = bmat\log(y_sm);
    start_sm(~isfinite(start_sm))=eps;
    nls_sm=0*start_sm;
    parfor ii = 1:nvox
     data=double(y_sm(:, ii));
     if any(~isfinite(data))
         continue
     end
     nls_sm(:,ii) = lsqnonlin(@(x)ObjF(x,double(bmat),data),double(start_sm(:,ii)),[],[],options);
    end
    clear start_sm
    dt.nls_sm = W2K(nls_sm(2:22, :));     
   
    
    
    
    
     
    %% Step 3: Classify voxels 
    % Let's define the bad voxels, we expect AKC to be possitive in all
    % directions.
    
    [akc, Td, Tk] = AKC(dt.nls);
    mk.nls = mean(akc, 1);
    
    list.accept = all(akc>-eps);
    list.remove = any(akc<0);
    list.MK.negative = mk.nls<0;
    list.MK.outlier = mk.nls<-2;
    list.N = nvox;
  
    
    Dsum = dt.nls(1,:).^2 + dt.nls(4,:).^2 + dt.nls(6, :).^2 + 2*dt.nls(2, :).^2 + 2*dt.nls(3, :).^2 + 2*dt.nls(5, :).^2;
    md.nls = (dt.nls(1,:) + dt.nls(4,:) + dt.nls(6,:))/3;
    md.nls_sm = (dt.nls_sm(1,:) + dt.nls_sm(4,:) + dt.nls_sm(6,:))/3;
    
    
    
    %% Step 4: Powder Kurtosis
    % Computation of the powder kurtosis.
    % We must include strategy to correct for this automatically. TODO.
    
    bvals = unique(bval);  
    for j = 1:numel(bvals)
        idx = bval==bvals(j);
        s(j, :) = mean(y(idx, :), 1);
    end
   
    s(s<1) = 1;
    iB = [ones(size(bvals)), -bvals, 1/6*bvals.^2];
    x = pinv(iB)*log(s);

    iD = x(2, :); 
    iW = x(3, :);
    iK = iW ./ iD.^2;
    
    list.IK.negative = iK<0; 
    list.IK.outlier = iK<-2; 

    
    
    if extraRobustness
        
        % we've experienced some data sets in which the powder kurtosis shows a
        % few black voxels. We aim to improve the robustness by
        % replacing those few black voxels by the estimate of smoothed data. 
    
    
        S = vectorize(s, mask); 
        S(isnan(S)) = 0;
        h = fspecial('gaussian', [5, 5], 1.25/(2*sqrt(2*log(2)))); 
        for k = 1:size(S, 3)
            for l = 1:size(S, 4)
                S(:,:,k,l) = filter2(h, S(:,:,k,l), 'same');  
            end
        end
        tmpmask = vectorize(iK, mask)<0;
        s = vectorize(S, tmpmask);
        s(s<1) = 1;
        x = pinv(iB)*log(s);
        iK(list.IK.negative) = x(3, :) ./ x(2, :).^2;

        list.IK.negative = iK<0; 
        list.IK.outlier = iK<-2;
    end
    
    
% %     %% Step 5: Voxel Quality Transfer: Predict W0,W2,W4 from pca on the DWI
% %     idx = list.accept & ~list.IK.negative;
% %     trainingSize=min([5e4 sum(idx)]);
% %     Dij_nls=nls(2:7,:);
% %     Wijkl_nls=nls(8:22,:)./md.nls.^2;
% %     Dlm_nls=wrapper_STF_decomposition_backandforth(Dij_nls,'cart2STF',[],'JV');
% %     Wlm_nls=wrapper_STF_decomposition_backandforth(Wijkl_nls,'cart2STF',[],'JV');
% %     nls_D0=vectorize(Dlm_nls(1,:)/sqrt(4*pi),mask);
% %     nls_D2=vectorize(sqrt(sum(Dlm_nls(2:6,:).^2,1))/sqrt(5*4*pi),mask);
% %     nls_W0=vectorize(Wlm_nls(1,:)/sqrt(4*pi),mask);
% %     nls_W2=vectorize(sqrt(sum(Wlm_nls(2:6,:).^2,1))/sqrt(5*4*pi),mask);
% %     nls_W4=vectorize(sqrt(sum(Wlm_nls(7:15,:).^2,1))/sqrt(9*4*pi),mask);
% % 
% %     
% %     
% %     Dij_nls_sm=nls_sm(2:7,:);
% %     Wijkl_nls_sm=nls_sm(8:22,:)./md.nls_sm.^2;
% %     Dlm_nls=wrapper_STF_decomposition_backandforth(Dij_nls_sm,'cart2STF',[],'JV');
% %     Wlm_nls=wrapper_STF_decomposition_backandforth(Wijkl_nls_sm,'cart2STF',[],'JV');
% %     nls_sm_D0=vectorize(Dlm_nls(1,:)/sqrt(4*pi),mask);
% %     nls_sm_D2=vectorize(sqrt(sum(Dlm_nls(2:6,:).^2,1))/sqrt(5*4*pi),mask);
% %     nls_sm_W0=vectorize(Wlm_nls(1,:)/sqrt(4*pi),mask);
% %     nls_sm_W2=vectorize(sqrt(sum(Wlm_nls(2:6,:).^2,1))/sqrt(5*4*pi),mask);
% %     nls_sm_W4=vectorize(sqrt(sum(Wlm_nls(7:15,:).^2,1))/sqrt(9*4*pi),mask);
% % 
% %     
% %     % Save predictions
% %     DlWlpred = cat(4,nls_sm_D0,nls_sm_D2,nls_sm_W0,nls_sm_W2,nls_sm_W4);
% %     
% % %     IMGUI(cat(4,nls_sm_D0,nls_D0))
% % %     IMGUI(cat(4,nls_sm_D2,nls_D2))
% % %     IMGUI(cat(4,nls_sm_W0,nls_W0))
% % %     IMGUI(cat(4,nls_sm_W2,nls_W2))
% % %     IMGUI(cat(4,nls_sm_W4,nls_W4))
    
    %% Step 6 Compute Regularization term
    % get proxy from inputs
    D0_proxy=vectorize(DlWlproxy(:,:,:,1),mask);
    D2_proxy=vectorize(DlWlproxy(:,:,:,2),mask);
    W0_proxy=vectorize(DlWlproxy(:,:,:,3),mask);
    W2_proxy=vectorize(DlWlproxy(:,:,:,4),mask);
    W4_proxy=vectorize(DlWlproxy(:,:,:,5),mask);

    
    Dij_nls=nls(2:7,:);
    Wijkl_nls=nls(8:22,:)./md.nls.^2;
    Dlm_nls=wrapper_STF_decomposition_backandforth(Dij_nls,'cart2STF',[],'JV');
    Wlm_nls=wrapper_STF_decomposition_backandforth(Wijkl_nls,'cart2STF',[],'JV');
    nls_D0=vectorize(Dlm_nls(1,:)/sqrt(4*pi),mask);
    nls_D2=vectorize(sqrt(sum(Dlm_nls(2:6,:).^2,1))/sqrt(5*4*pi),mask);
    nls_W0=vectorize(Wlm_nls(1,:)/sqrt(4*pi),mask);
    nls_W2=vectorize(sqrt(sum(Wlm_nls(2:6,:).^2,1))/sqrt(5*4*pi),mask);
    nls_W4=vectorize(sqrt(sum(Wlm_nls(7:15,:).^2,1))/sqrt(9*4*pi),mask);
    
    D0=vectorize(nls_D0,mask);
    D2=vectorize(nls_D2,mask);
    W0=vectorize(nls_W0,mask);
    W2=vectorize(nls_W2,mask);
    W4=vectorize(nls_W4,mask);
    
    % Tuning the "alpha" parameter is a challenge for most regularized estimators. That's not different here. 
    mse.nls = nanmedian(sum((y_hat-y).^2)/size(y, 1));
%     mse.mk = nanmedian((mk.predicted(idx) - mk.nls(idx)).^2);
    mse.d0 = nanmedian((D0_proxy(idx) - D0(idx)).^2);
    mse.d2 = nanmedian((D2_proxy(idx) - D2(idx)).^2);
    mse.w0 = nanmedian((W0_proxy(idx) - W0(idx)).^2);
    mse.w2 = nanmedian((W2_proxy(idx) - W2(idx)).^2);
    mse.w4 = nanmedian((W4_proxy(idx) - W4(idx)).^2);
%     alpha = 0.1 *mse.nls ./ mse.mk; % can we increase?
    alpha_d0 = 0.2 *mse.nls ./ mse.d0; % can we increase?
    alpha_d2 = 0.1 *mse.nls ./ mse.d2; % can we increase?
    alpha_w0 = 0.5 *mse.nls ./ mse.w0; % can we increase?
    alpha_w2 = 0.2 *mse.nls ./ mse.w2; % can we increase?
    alpha_w4 = 0.01 *mse.nls ./ mse.w4; % can we increase?
    DlWl_alpha=[alpha_d0 alpha_d2 alpha_w0 alpha_w2 alpha_w4];
    Wl_alpha=[alpha_w0 alpha_w2 alpha_w4];
    %% Step 7: Regularized fit with constrained fit as a starting point.
    %  regularized 
    [C, d] = createConstraints();
    opts = optimset('Display', 'off', 'Algorithm', 'interior-point', 'MaxIter', 10000, 'TolCon', 1e-8, 'TolFun', 1e-8, 'TolX', 1e-8);
    options = optimset('fminunc'); 
    options = optimset(options,'Jacobian','on','TolFun',1e-8,'TolX', 1e-8,'MaxIter',10000,'Display','off');
    start = nls;
    idx_refit = find(list.remove);
    for ii = 1:numel(idx_refit)   
        try
            start(:, idx_refit(ii)) = lsqlin(double(bmat),double((log(y(:, idx_refit(ii))))),-C, d, [],[],[],[],[],opts);
        catch
        end
    end
    
    CSphase=0;
    if length(alpha_tuning)==6
        dirs256 = get256dirs();
        bmat_stf_256=get_even_SH(dirs256,4,CSphase);
        alpha_tuning=alpha_tuning(1:5);
    end
    
    DlWl_alpha=[alpha_d0 alpha_d2 alpha_w0 alpha_w2 alpha_w4].*alpha_tuning;
    D_ij=start(2:7,:);
    W_ijkl_hat=start(8:22,:);
    [Dlm] = IrreducibleDecomposition_rank2_tensor(D_ij([1 4 6 2 3 5],:));
    [Wlm_hat] = IrreducibleDecomposition_rank4_tensor(W_ijkl_hat([1 11 15 4 6 2 3 5 13 7 8 12 9 10 14],:));
    start_lm=[start(1,:);Dlm;Wlm_hat];
    
    Y_LM_matrix = get_even_SH(grad(:,1:3),4,CSphase);
    Y00=Y_LM_matrix(:,1);
    Y2m=Y_LM_matrix(:,2:6);
    Y4m=Y_LM_matrix(:,7:15);
    reg = start_lm;
    % Do not re-fit b0
    b0_nls=exp(min(15,nls(1, :)))/scaleFact;
    bmat_stf=[-grad(:,4).*Y00 -grad(:,4).*Y2m 1/6*grad(:,4).^2.*Y00 1/6*grad(:,4).^2.*Y2m 1/6*grad(:,4).^2.*Y4m];
        
    Niter=zeros(numel(idx_refit),1);
    fval=zeros(numel(idx_refit),2);
    for ii = 1:numel(idx_refit)   
        if any(~isfinite(double(start_lm(:,idx_refit(ii)))))||max(abs(double(y(:, idx_refit(ii)))))<1e-10
            continue
        end
        DlWl_proxies=[D0_proxy(idx_refit(ii)) D2_proxy(idx_refit(ii)) W0_proxy(idx_refit(ii)) W2_proxy(idx_refit(ii)) W4_proxy(idx_refit(ii))]';  
%         Wl_proxies=[W0_proxy(idx_refit(ii)) W2_proxy(idx_refit(ii)) W4_proxy(idx_refit(ii))]';  
        x0=double(start_lm(2:22,idx_refit(ii)));
        
        if exist('bmat_stf_256','var')
            objFun=@(x)regObjF_DlWl_nob0_Wpositive(x,double(bmat_stf),double(y(:, idx_refit(ii))), DlWl_proxies, DlWl_alpha,scaleFact*b0_nls(idx_refit(ii)),double(bmat_stf_256));
        else
            objFun=@(x)regObjF_DlWl_nob0(x,double(bmat_stf),double(y(:, idx_refit(ii))), DlWl_proxies, DlWl_alpha,scaleFact*b0_nls(idx_refit(ii)));
        end
        [sol,~,exitflag(ii),output] = fminunc(objFun,x0,options);
        fval(ii,:)=[objFun(x0) objFun(sol)];
        Niter(ii)=output.iterations;
        reg(2:22,idx_refit(ii)) = sol;
    end    
    [Dij] = STF2cart_rank2_tensor(reg(2:7, :));
    [Wijkl] = STF2cart_rank4_tensor(reg(8:22, :));

    %% Step 9: write output
    dt = vectorize(W2K([Dij([1 4 5 2 6 3],:);Wijkl([1 6 7 4 8 5 10 11 13 14 2 12 9 15 3],:)]), mask);
    b0 = vectorize(b0_nls, mask);    

%     %% Step 8: Quality Control
%     [akc, Td, Tk] = AKC(dt.reg);    
%     list.notSolvedAfterFirstTry = mean(akc, 1)<0;    
% %     %% Step 8: Update the remaining black voxels with stronger regularization and constrained fit
% %     if extraRobustness 
% %         
% %         options = optimset('fmincon'); 
% %         options = optimset(options,'Jacobian','on','TolFun',1e-8,'TolX', ...
% %                         1e-8,'MaxIter',10000,'Display','off');
% %           
% %         idx = find(list.notSolvedAfterFirstTry);
% %         for  ii = 1:numel(idx)   
% %             try
% % %                 reg(:,idx(ii)) = fmincon(@(x)regObjF(x,double(bmat),double(y(:, idx(ii))), mk.predicted(idx(ii)), Td, Tk, 10*alpha),double(start(:,idx(ii))), -C, d, [],[], [], [], [], options);
% %                 reg(:,idx(ii)) = fmincon(@(x)regObjF_W(x,double(bmat),double(y(:, idx(ii))), W0_proxy(idx(ii)), 10*alpha),double(start(:,idx(ii))), -C, d, [],[], [], [], [], options);
% %             catch
% %                 reg(:,idx(ii)) = NaN;
% %             end
% %         end
% %         dt.reg = W2K(reg(2:22, :));
% %         b0.reg = exp(min(15,reg(1, :)))/scaleFact;
% % 
% %     end
%     
%     %% Step 9: write output
%     dt = vectorize(dt.reg, mask);
%     b0 = vectorize(b0.reg, mask);    
end

function [C, d] = createConstraints()

    dir = [0.382517725304416 -0.748614094922528 0.541532202838631;-0.266039846728327 0.963894740823927 0.0113898448076587;-0.128563443377023 0.800867622029815 0.584878186472493;0.691696803043553 0.485345502199397 0.534785261721136;0.776929615593511 -0.627201085171846 0.0547646891069225;-0.314229418625565 0.891550800503996 0.326185595314880;-0.984699447847175 0.0338717154803320 0.170937720529697;0.729869942283584 0.134539815263771 0.670215566411097;0.0491118066650937 0.613801560175467 0.787931262974286;0.615167937666214 0.786762996419759 0.0507187926916626;-0.504930428375015 -0.548805916531712 0.666226184175323;0.514775318788445 0.353967263592948 0.780841563616317;-0.306616550550256 0.577152309970495 0.756889359169743;-0.644455563348338 0.445243323148325 0.621639292565402;0.888177219438464 0.244852048242751 0.388829913126405;-0.115867623474531 0.331617270421714 0.936271691224516;0.312724982544119 -0.262437525100548 0.912868901163732;-0.348318356641730 -0.328727572647744 0.877845376707953;0.622993255900061 -0.170127400464004 0.763502502100944;-0.870285082134136 -0.0832402149084147 0.485463636575162;0.879693901263504 -0.0847887289384472 0.467920411528283;0.375735817168569 0.624963320743740 0.684283160264539;-0.763508679267313 0.569075961898599 0.305298290648126;0.895786299773340 -0.371201461149296 0.244492086536589;0.431182280410191 0.0594580709470589 0.900303603713504;-0.927083085686508 -0.288337655567580 0.239537781186549;0.208899044398678 0.833216349905585 0.511968459477078;-0.671275756453876 -0.252498824452251 0.696873878436771;-0.385511621254227 -0.908766073079027 0.159765497834991;-0.501120479467597 0.703268192077924 0.504273849281930;-0.578440272143465 0.801933906922628 0.149361509400528;0.986601726072896 -0.0507533113495985 0.155052041254001;0.0472262384668294 -0.790665651327184 0.610424041311968;0.957038035873056 0.279601625450131 0.0768188058868917;-0.497573767044291 -0.0342706449545790 0.866744408256408;0.537095370960702 0.746985750118871 0.391842891490881;0.174500355902118 -0.559258086805823 0.810419655568845;-0.0648836431571087 -0.997212186937296 0.0368506048036694;-0.200896533381969 0.00230971954655800 0.979609742739793;-0.436037609875685 0.290696319170598 0.851684714430502;0.332217034685261 0.924381756972555 0.187483890618005;0.115538097684954 0.0265470728743124 0.992948236770250;-0.167448247267712 -0.594347070791611 0.786582890691377;-0.931940478288593 0.352679013976080 0.0842879471104161;0.749660835628331 -0.375215305423717 0.545180801295148;-0.112213298457421 -0.929475988744578 0.351400856567817;-0.596160909541517 -0.730179079768923 0.333812344592648;0.211077351410955 0.350067854669890 0.912632921194583;-0.325168748302559 -0.780267672863407 0.534273004943793;-0.717210613875971 0.128994202815781 0.684813427864535;-0.0218381924490005 -0.303713916706869 0.952512965869302;0.213275433291729 -0.924500457406975 0.315931153589701;-0.810453788321924 -0.574858547954973 0.112704511168551;0.665549405791414 -0.637998127347686 0.387301404530814;0.489520321770316 -0.495410872500818 0.717591751612200;0.514060443295042 -0.837385561080287 0.185815184346054;-0.757892441488769 -0.466692556842526 0.455847676885577;0.00471100435105065 0.958734616992657 0.284263505603424;0.800137357904460 0.555340864988139 0.226664360144898;-0.872328992553570 0.265326196661285 0.410663046932313;];
    ndir = size(dir, 1);
    [W_ind, W_cnt] = createTensorOrder(4);    
    C = [zeros(ndir, 7), W_cnt(ones(ndir, 1), :).*dir(:,W_ind(:, 1)).*dir(:,W_ind(:, 2)).*dir(:,W_ind(:, 3)).*dir(:,W_ind(:, 4))];
    d = zeros([1, size(C, 1)]);

end

function b = bmatrix(grad)

    bval = grad(:, 4);
    ndwis = size(grad, 1);
    [D_ind, D_cnt] = createTensorOrder(2);
    [W_ind, W_cnt] = createTensorOrder(4);
    
    bS = ones(ndwis, 1);
    bD = D_cnt(ones(ndwis, 1), :).*grad(:,D_ind(:, 1)).*grad(:,D_ind(:, 2));
    bW = W_cnt(ones(ndwis, 1), :).*grad(:,W_ind(:, 1)).*grad(:,W_ind(:, 2)).*grad(:,W_ind(:, 3)).*grad(:,W_ind(:, 4));
    
    b = [bS, -bval(:, ones(1, 6)).*bD, (bval(:, ones(1, 15)).^2/6).*bW];
                
                 
end

function [X, cnt] = createTensorOrder(order)
    X = nchoosek(kron([1, 2, 3], ones(1, order)), order);
    X = unique(X, 'rows');
    for i = 1:size(X, 1)
        cnt(i) = factorial(order) / factorial(nnz(X(i, :) ==1))/ factorial(nnz(X(i, :) ==2))/ factorial(nnz(X(i, :) ==3));
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


function [akc, Td, Tk] = AKC(dt)

        dir = get256dirs();
        [D_ind, D_cnt] = createTensorOrder(2);
        [W_ind, W_cnt] = createTensorOrder(4);
        
        Td = (dir(:,D_ind(1:6,1)).*dir(:,D_ind(1:6,2))) * diag(D_cnt);
        adc = Td * dt(1:6, :);            
        md = sum(dt([1 4 6],:),1)/3;
        
        Tk = (prod(reshape(dir(:,W_ind),[],15,4),3))*diag(W_cnt);
        akc =  Tk*dt(7:21, :);
        akc = (akc .* repmat(md.^2, [size(adc, 1), 1]))./(adc.^2);
            
            
end

function [E, J] = regObjF(x,bmat, y, iW, Td, Tk, alpha)
    
        S = exp(bmat*x);
        ny = numel(y);
        E1 = sum((y - S).^2)/ny;
        
        tmp = -2*(y-S).*S;
        J1 = sum(tmp(:, ones(1, 22)).*bmat)/ny;
       
        g = Td*x(2:7); dg = Td;             
        f = Tk*x(8:22); df = Tk;
        
        gg = 1 ./(g.^2); gg = gg(:, ones(1, 15));
        ff = f ./ (g.^3); ff = [ff, ff, ff, ff, ff, ff]; 
        
        mw = sum(f./g.^2)/256;

        E2 = (iW-mw).^2;  

        term1 = [zeros(1, 7), sum(gg.*df)/256];
        term2 = [0, sum(2*ff.*dg)/256, zeros(1, 15)];

        J2 =  -2*(iW-mw)* (term1 - term2); 
        
        E = E1 + alpha*E2;
        J = J1 + alpha*J2;
        
        E = E';
end

function [E, J] = regObjF_W0(x,bmat, y, MW_proxy, alpha)
    
        S = exp(bmat*x);
        ny = numel(y);
        E1 = sum((y - S).^2)/ny;
        
        tmp = -2*(y-S).*S;
        J1 = sum(tmp(:, ones(1, 22)).*bmat)/ny;
        
        MD=sum(x(1+[1 4 6]))/3;
        mw = 1/5*sum(x(7+[1 11 15 4 4 6 6 13 13]))/MD^2;

        E2 = (MW_proxy-mw).^2;  
        
        a=1/5*sum(x(7+[1 11 15 4 4 6 6 13 13]))*(-2)/MD^3*1/3;
        D_term=[a 0 0 a 0 a];
        W_term=[1 0 0 2 0 2 0 0 0 0 1 0 2 0 1]/(5*MD^2);
        J2 =  -2*(MW_proxy-mw) * [0, D_term, W_term]; 
        
        E = E1 + alpha*E2;
        J = J1 + alpha*J2;
        
        E = E';
end

function [E] = regObjF_Wl(x,bmat, y, MW_proxy, alpha)
    
        S = exp(bmat*x);
        ny = numel(y);
        E1 = sum((y - S).^2)/ny;
        
        tmp = -2*(y-S).*S;
        J1 = sum(tmp(:, ones(1, 22)).*bmat)/ny;
        
        MD=x(2)/sqrt(4*pi);
        wl_hat=[x(8)/sqrt(4*pi) sqrt(sum(x(9:13).^2)/(5*4*pi)) sqrt(sum(x(14:22).^2)/(9*4*pi))]';
        wl=wl_hat/MD^2;       
        E2 = (MW_proxy-wl).^2;  
        
        sqrt_sq_sum_w2m_hat=sqrt(sum(x(9:13).^2));
        sqrt_sq_sum_w4m_hat=sqrt(sum(x(14:22).^2));
        w2m_factor=sqrt(4*pi/5)/x(2)^2;
        w4m_factor=sqrt(4*pi/9)/x(2)^2;
        
        dW0_ddt=[0,(-2)*sqrt(4*pi)/x(2)^3*abs(x(8))            ,0,0,0,0,0, sqrt(4*pi)/x(2)^2,           0,0,0,0,0,    0,0,0,0,0,0,0,0,0];
        dW2_ddt=[0,(-2)*sqrt(4*pi/5)/x(2)^3*sqrt_sq_sum_w2m_hat,0,0,0,0,0,                 0, w2m_factor*x(9:13)',    0,0,0,0,0,0,0,0,0];
        dW4_ddt=[0,(-2)*sqrt(4*pi/9)/x(2)^3*sqrt_sq_sum_w4m_hat,0,0,0,0,0,                 0,           0,0,0,0,0, w4m_factor*x(14:22)'];
        dWl_ddt= [dW0_ddt;dW2_ddt;dW4_ddt]; 
        J2 =  -2*repmat((MW_proxy-wl),1,22) .* dWl_ddt; 
        
        E = E1 + alpha*E2;
        J = J1 + alpha*J2;
        
        E = E';
end

function [E] = regObjF_DlWl_nob0(x,bmat, y, DlWl_proxy, alpha,b0)
        S = b0*exp(bmat*x);
        ny = numel(y);
        E1 = sum((y - S).^2)/ny;
        dl=[x(1)/sqrt(4*pi) sqrt(sum(x(2:6).^2)/(5*4*pi))]';
        MD=x(1)/sqrt(4*pi);
        wl_hat=[x(7)/sqrt(4*pi) sqrt(sum(x(8:12).^2)/(5*4*pi)) sqrt(sum(x(13:21).^2)/(9*4*pi))]';
        wl=wl_hat/MD^2;  
        dlwl=[dl;wl];
        E2 = (DlWl_proxy-dlwl).^2;  
        E = E1 + alpha*E2;
        E = E';
%         % Getting derivatives
%         tmp = -2*(y-S).*S;
%         J1 = sum(tmp(:, ones(1, 21)).*bmat)/ny;
%         
%         sqrt_sq_sum_w2m_hat=sqrt(sum(x(9:13).^2));
%         sqrt_sq_sum_w4m_hat=sqrt(sum(x(14:22).^2));
%         w2m_factor=sqrt(4*pi/5)/x(2)^2;
%         w4m_factor=sqrt(4*pi/9)/x(2)^2;
%         
%         dD0_ddt=[1/sqrt(4*pi)                                ,0,0,0,0,0,                 0,           0,0,0,0,0,    0,0,0,0,0,0,0,0,0];
%         dD2_ddt=[0                                           ,0,0,0,0,0,                 0,           0,0,0,0,0,    0,0,0,0,0,0,0,0,0];
%         dW0_ddt=[(-2)*sqrt(4*pi)/x(2)^3*abs(x(8))            ,0,0,0,0,0, sqrt(4*pi)/x(1)^2,           0,0,0,0,0,    0,0,0,0,0,0,0,0,0];
%         dW2_ddt=[(-2)*sqrt(4*pi/5)/x(2)^3*sqrt_sq_sum_w2m_hat,0,0,0,0,0,                 0, w2m_factor*x(9:13)',    0,0,0,0,0,0,0,0,0];
%         dW4_ddt=[(-2)*sqrt(4*pi/9)/x(2)^3*sqrt_sq_sum_w4m_hat,0,0,0,0,0,                 0,           0,0,0,0,0, w4m_factor*x(14:22)'];
%         dWl_ddt= [dW0_ddt;dW2_ddt;dW4_ddt]; 
%         J2 =  -2*repmat((MW_proxy-wl),1,22) .* dWl_ddt; 
% 
%         J = J1 + alpha*J2;
end

function [E] = regObjF_DlWl_nob0_Wpositive(x,bmat, y, DlWl_proxy, alpha,b0,dirs256_bmat_W)
        S = b0*exp(bmat*x);
        ny = numel(y);
        E1 = sum((y - S).^2)/ny;
        dl=[x(1)/sqrt(4*pi) sqrt(sum(x(2:6).^2)/(5*4*pi))]';
        MD=x(1)/sqrt(4*pi);
        wl_hat=[x(7)/sqrt(4*pi) sqrt(sum(x(8:12).^2)/(5*4*pi)) sqrt(sum(x(13:21).^2)/(9*4*pi))]';
        wl=wl_hat/MD^2;  
        dlwl=[dl;wl];
        E2 = (DlWl_proxy-dlwl).^2;  
        
        Wproj=dirs256_bmat_W*x(7:21);
        Wproj(Wproj>0)=0;
        
        E = E1 + alpha*E2 - alpha(1)*sum(Wproj);
        E = E';
%         % Getting derivatives
%         tmp = -2*(y-S).*S;
%         J1 = sum(tmp(:, ones(1, 21)).*bmat)/ny;
%         
%         sqrt_sq_sum_w2m_hat=sqrt(sum(x(9:13).^2));
%         sqrt_sq_sum_w4m_hat=sqrt(sum(x(14:22).^2));
%         w2m_factor=sqrt(4*pi/5)/x(2)^2;
%         w4m_factor=sqrt(4*pi/9)/x(2)^2;
%         
%         dD0_ddt=[1/sqrt(4*pi)                                ,0,0,0,0,0,                 0,           0,0,0,0,0,    0,0,0,0,0,0,0,0,0];
%         dD2_ddt=[0                                           ,0,0,0,0,0,                 0,           0,0,0,0,0,    0,0,0,0,0,0,0,0,0];
%         dW0_ddt=[(-2)*sqrt(4*pi)/x(2)^3*abs(x(8))            ,0,0,0,0,0, sqrt(4*pi)/x(1)^2,           0,0,0,0,0,    0,0,0,0,0,0,0,0,0];
%         dW2_ddt=[(-2)*sqrt(4*pi/5)/x(2)^3*sqrt_sq_sum_w2m_hat,0,0,0,0,0,                 0, w2m_factor*x(9:13)',    0,0,0,0,0,0,0,0,0];
%         dW4_ddt=[(-2)*sqrt(4*pi/9)/x(2)^3*sqrt_sq_sum_w4m_hat,0,0,0,0,0,                 0,           0,0,0,0,0, w4m_factor*x(14:22)'];
%         dWl_ddt= [dW0_ddt;dW2_ddt;dW4_ddt]; 
%         J2 =  -2*repmat((MW_proxy-wl),1,22) .* dWl_ddt; 
% 
%         J = J1 + alpha*J2;
end


function [F,J] = ObjF(dt,b,dwi)
    dwi_hat = exp(b*dt);
    F = dwi_hat-dwi;
    J = dwi_hat(:,ones(1, size(b,2))).*b;
end


function dt = W2K(dt)
     D_apprSq = 1./(sum(dt([1 4 6],:),1)/3).^2;
     dt(7:21,:) = dt(7:21,:) .* D_apprSq(ones(15,1),:);
end
        

function [X, combs] = createDesignMatrix(params)

    [Ntraining, ~] = size(params);

    combs = [0 0 0; ...
                 1 0 0; ...
                 2 0 0; ...
                 3 0 0; ...
                 0 1 0; ...
                 0 2 0; ...
                 0 3 0; ...
                 0 0 1; ...
                 0 0 2; ...
                 0 0 3; ...
                 1 1 0; ...
                 0 1 1; ...
                 1 0 1];
             
    ncombs = size(combs, 1);
    for i = 1:Ntraining
        for j = 1:ncombs
            X(i, j) = prod(params(i, :).^combs(j, :));
        end
    end
end

function dir = get256dirs()
 dir = [0         0    1.0000
    0.5924         0    0.8056
   -0.7191   -0.1575   -0.6768
   -0.9151   -0.3479    0.2040
    0.5535    0.2437    0.7964
   -0.0844    0.9609   -0.2636
    0.9512   -0.3015    0.0651
   -0.4225    0.8984    0.1202
    0.5916   -0.6396    0.4909
    0.3172    0.8818   -0.3489
   -0.1988   -0.6687    0.7164
   -0.2735    0.3047   -0.9123
    0.9714   -0.1171    0.2066
   -0.5215   -0.4013    0.7530
   -0.3978   -0.9131   -0.0897
    0.2680    0.8196    0.5063
   -0.6824   -0.6532   -0.3281
    0.4748   -0.7261   -0.4973
    0.4504   -0.4036    0.7964
   -0.5551   -0.8034   -0.2153
    0.0455   -0.2169    0.9751
    0.0483    0.5845    0.8099
   -0.1909   -0.1544   -0.9694
    0.8383    0.5084    0.1969
   -0.2464    0.1148    0.9623
   -0.7458    0.6318    0.2114
   -0.0080   -0.9831   -0.1828
   -0.2630    0.5386   -0.8005
   -0.0507    0.6425   -0.7646
    0.4476   -0.8877    0.1081
   -0.5627    0.7710    0.2982
   -0.3790    0.7774   -0.5020
   -0.6217    0.4586   -0.6350
   -0.1506    0.8688   -0.4718
   -0.4579    0.2131    0.8631
   -0.8349   -0.2124    0.5077
    0.7682   -0.1732   -0.6163
    0.0997   -0.7168   -0.6901
    0.0386   -0.2146   -0.9759
    0.9312    0.1655   -0.3249
    0.9151    0.3053    0.2634
    0.8081    0.5289   -0.2593
   -0.3632   -0.9225    0.1305
    0.2709   -0.3327   -0.9033
   -0.1942   -0.9790   -0.0623
    0.6302   -0.7641    0.1377
   -0.6948   -0.3137    0.6471
   -0.6596   -0.6452    0.3854
   -0.9454    0.2713    0.1805
   -0.2586   -0.7957    0.5477
   -0.3576    0.6511    0.6695
   -0.8490   -0.5275    0.0328
    0.3830    0.2499   -0.8893
    0.8804   -0.2392   -0.4095
    0.4321   -0.4475   -0.7829
   -0.5821   -0.1656    0.7961
    0.3963    0.6637    0.6344
   -0.7222   -0.6855   -0.0929
    0.2130   -0.9650   -0.1527
    0.4737    0.7367   -0.4825
   -0.9956    0.0891    0.0278
   -0.5178    0.7899   -0.3287
   -0.8906    0.1431   -0.4317
    0.2431   -0.9670    0.0764
   -0.6812   -0.3807   -0.6254
   -0.1091   -0.5141    0.8507
   -0.2206    0.7274   -0.6498
    0.8359    0.2674    0.4794
    0.9873    0.1103    0.1147
    0.7471    0.0659   -0.6615
    0.6119   -0.2508    0.7502
   -0.6191    0.0776    0.7815
    0.7663   -0.4739    0.4339
   -0.5699    0.5369    0.6220
    0.0232   -0.9989    0.0401
    0.0671   -0.4207   -0.9047
   -0.2145    0.5538    0.8045
    0.8554   -0.4894    0.1698
   -0.7912   -0.4194    0.4450
   -0.2341    0.0754   -0.9693
   -0.7725    0.6346   -0.0216
    0.0228    0.7946   -0.6067
    0.7461   -0.3966   -0.5348
   -0.4045   -0.0837   -0.9107
   -0.4364    0.6084   -0.6629
    0.6177   -0.3175   -0.7195
   -0.4301   -0.0198    0.9026
   -0.1489   -0.9706    0.1892
    0.0879    0.9070   -0.4117
   -0.7764   -0.4707   -0.4190
    0.9850    0.1352   -0.1073
   -0.1581   -0.3154    0.9357
    0.8938   -0.3246    0.3096
    0.8358   -0.4464   -0.3197
    0.4943    0.4679    0.7327
   -0.3095    0.9015   -0.3024
   -0.3363   -0.8942   -0.2956
   -0.1271   -0.9274   -0.3519
    0.3523   -0.8717   -0.3407
    0.7188   -0.6321    0.2895
   -0.7447    0.0924   -0.6610
    0.1622    0.7186    0.6762
   -0.9406   -0.0829   -0.3293
   -0.1229    0.9204    0.3712
   -0.8802    0.4668    0.0856
   -0.2062   -0.1035    0.9730
   -0.4861   -0.7586   -0.4338
   -0.6138    0.7851    0.0827
    0.8476    0.0504    0.5282
    0.3236    0.4698   -0.8213
   -0.7053   -0.6935    0.1473
    0.1511    0.3778    0.9135
    0.6011    0.5847    0.5448
    0.3610    0.3183    0.8766
    0.9432    0.3304    0.0341
    0.2423   -0.8079   -0.5372
    0.4431   -0.1578    0.8825
    0.6204    0.5320   -0.5763
   -0.2806   -0.5376   -0.7952
   -0.5279   -0.8071    0.2646
   -0.4214   -0.6159    0.6656
    0.6759   -0.5995   -0.4288
    0.5670    0.8232   -0.0295
   -0.0874    0.4284   -0.8994
    0.8780   -0.0192   -0.4782
    0.0166    0.8421    0.5391
   -0.7741    0.2931   -0.5610
    0.9636   -0.0579   -0.2611
         0         0   -1.0000
   -0.5924         0   -0.8056
    0.7191    0.1575    0.6768
    0.9151    0.3479   -0.2040
   -0.5535   -0.2437   -0.7964
    0.0844   -0.9609    0.2636
   -0.9512    0.3015   -0.0651
    0.4225   -0.8984   -0.1202
   -0.5916    0.6396   -0.4909
   -0.3172   -0.8818    0.3489
    0.1988    0.6687   -0.7164
    0.2735   -0.3047    0.9123
   -0.9714    0.1171   -0.2066
    0.5215    0.4013   -0.7530
    0.3978    0.9131    0.0897
   -0.2680   -0.8196   -0.5063
    0.6824    0.6532    0.3281
   -0.4748    0.7261    0.4973
   -0.4504    0.4036   -0.7964
    0.5551    0.8034    0.2153
   -0.0455    0.2169   -0.9751
   -0.0483   -0.5845   -0.8099
    0.1909    0.1544    0.9694
   -0.8383   -0.5084   -0.1969
    0.2464   -0.1148   -0.9623
    0.7458   -0.6318   -0.2114
    0.0080    0.9831    0.1828
    0.2630   -0.5386    0.8005
    0.0507   -0.6425    0.7646
   -0.4476    0.8877   -0.1081
    0.5627   -0.7710   -0.2982
    0.3790   -0.7774    0.5020
    0.6217   -0.4586    0.6350
    0.1506   -0.8688    0.4718
    0.4579   -0.2131   -0.8631
    0.8349    0.2124   -0.5077
   -0.7682    0.1732    0.6163
   -0.0997    0.7168    0.6901
   -0.0386    0.2146    0.9759
   -0.9312   -0.1655    0.3249
   -0.9151   -0.3053   -0.2634
   -0.8081   -0.5289    0.2593
    0.3632    0.9225   -0.1305
   -0.2709    0.3327    0.9033
    0.1942    0.9790    0.0623
   -0.6302    0.7641   -0.1377
    0.6948    0.3137   -0.6471
    0.6596    0.6452   -0.3854
    0.9454   -0.2713   -0.1805
    0.2586    0.7957   -0.5477
    0.3576   -0.6511   -0.6695
    0.8490    0.5275   -0.0328
   -0.3830   -0.2499    0.8893
   -0.8804    0.2392    0.4095
   -0.4321    0.4475    0.7829
    0.5821    0.1656   -0.7961
   -0.3963   -0.6637   -0.6344
    0.7222    0.6855    0.0929
   -0.2130    0.9650    0.1527
   -0.4737   -0.7367    0.4825
    0.9956   -0.0891   -0.0278
    0.5178   -0.7899    0.3287
    0.8906   -0.1431    0.4317
   -0.2431    0.9670   -0.0764
    0.6812    0.3807    0.6254
    0.1091    0.5141   -0.8507
    0.2206   -0.7274    0.6498
   -0.8359   -0.2674   -0.4794
   -0.9873   -0.1103   -0.1147
   -0.7471   -0.0659    0.6615
   -0.6119    0.2508   -0.7502
    0.6191   -0.0776   -0.7815
   -0.7663    0.4739   -0.4339
    0.5699   -0.5369   -0.6220
   -0.0232    0.9989   -0.0401
   -0.0671    0.4207    0.9047
    0.2145   -0.5538   -0.8045
   -0.8554    0.4894   -0.1698
    0.7912    0.4194   -0.4450
    0.2341   -0.0754    0.9693
    0.7725   -0.6346    0.0216
   -0.0228   -0.7946    0.6067
   -0.7461    0.3966    0.5348
    0.4045    0.0837    0.9107
    0.4364   -0.6084    0.6629
   -0.6177    0.3175    0.7195
    0.4301    0.0198   -0.9026
    0.1489    0.9706   -0.1892
   -0.0879   -0.9070    0.4117
    0.7764    0.4707    0.4190
   -0.9850   -0.1352    0.1073
    0.1581    0.3154   -0.9357
   -0.8938    0.3246   -0.3096
   -0.8358    0.4464    0.3197
   -0.4943   -0.4679   -0.7327
    0.3095   -0.9015    0.3024
    0.3363    0.8942    0.2956
    0.1271    0.9274    0.3519
   -0.3523    0.8717    0.3407
   -0.7188    0.6321   -0.2895
    0.7447   -0.0924    0.6610
   -0.1622   -0.7186   -0.6762
    0.9406    0.0829    0.3293
    0.1229   -0.9204   -0.3712
    0.8802   -0.4668   -0.0856
    0.2062    0.1035   -0.9730
    0.4861    0.7586    0.4338
    0.6138   -0.7851   -0.0827
   -0.8476   -0.0504   -0.5282
   -0.3236   -0.4698    0.8213
    0.7053    0.6935   -0.1473
   -0.1511   -0.3778   -0.9135
   -0.6011   -0.5847   -0.5448
   -0.3610   -0.3183   -0.8766
   -0.9432   -0.3304   -0.0341
   -0.2423    0.8079    0.5372
   -0.4431    0.1578   -0.8825
   -0.6204   -0.5320    0.5763
    0.2806    0.5376    0.7952
    0.5279    0.8071   -0.2646
    0.4214    0.6159   -0.6656
   -0.6759    0.5995    0.4288
   -0.5670   -0.8232    0.0295
    0.0874   -0.4284    0.8994
   -0.8780    0.0192    0.4782
   -0.0166   -0.8421   -0.5391
    0.7741   -0.2931    0.5610
   -0.9636    0.0579    0.2611];
end

function [Slm] = IrreducibleDecomposition_rank2_tensor(S_ij)
% S: 2D array containing [S_11;S_22;S_33;S_12;S_13;S_23];   
%
C0=sqrt(1/(4*pi));
C2=sqrt(5/(4*pi));
Slm=0*S_ij;
Y_ij_00 =C0*eye(3);
Y_ij_2m2=[0,0.546274215296040,0;0.546274215296040,0,0;0,0,0];
Y_ij_2m1=[0,0,0;0,0,-0.546274215296040;0,-0.546274215296040,0];
Y_ij_20 =[-0.315391565252520,0,0;0,-0.315391565252520,0;0,0,0.630783130505040];
Y_ij_21 =[0,0,-0.546274215296040;0,0,0;-0.546274215296040,0,0];
Y_ij_22 =[0.546274215296040,0,0;0,-0.546274215296040,0;0,0,0];
Y_ij_00_vector=Y_ij_00([1 5 9 4 7 8]).*[1 1 1 2 2 2];
Y_ij_2m2_vector=Y_ij_2m2([1 5 9 4 7 8]).*[1 1 1 2 2 2];
Y_ij_2m1_vector=Y_ij_2m1([1 5 9 4 7 8]).*[1 1 1 2 2 2];
Y_ij_20_vector=Y_ij_20([1 5 9 4 7 8]).*[1 1 1 2 2 2];
Y_ij_21_vector=Y_ij_21([1 5 9 4 7 8]).*[1 1 1 2 2 2];
Y_ij_22_vector=Y_ij_22([1 5 9 4 7 8]).*[1 1 1 2 2 2];
Slm(1,:) =1/3*(1/C0^2)*Y_ij_00_vector*S_ij;
% Slm(1,:)=1/3*(1/C0)*sum(S_ij(1:3,:),1);
Slm(2,:)=2/3*(1/C2^2)*Y_ij_2m2_vector*S_ij;
Slm(3,:)=2/3*(1/C2^2)*Y_ij_2m1_vector*S_ij;
Slm(4,:)=2/3*(1/C2^2)*Y_ij_20_vector*S_ij;
Slm(5,:)=2/3*(1/C2^2)*Y_ij_21_vector*S_ij;
Slm(6,:)=2/3*(1/C2^2)*Y_ij_22_vector*S_ij;
end


function [Slm] = IrreducibleDecomposition_rank4_tensor(S_ijkl)
% Sijkl: 2D array containing [S_1111;S_2222;S_3333;S_1122;S_1133;...
%                             S_1112;S_1113;S_1123;S_2233;S_2212;...
%                             S_2213;S_2223;S_3312;S_3313;S_3323];   
%
%
Slm=0*S_ijkl;
C0=sqrt(1/(4*pi));
C2=sqrt(5/(4*pi));
C4=sqrt(9/(4*pi));

Y_ijkl_00_vector = [0.282094791773878,0.282094791773878,0.282094791773878,0.564189583547756,0.564189583547756,0,0,0,0.564189583547756,0,0,0,0,0,0];
Y_ijkl_2m2_vector= [0,0,0,0,0,1.09254843059208,0,0,0,1.09254843059208,0,0,1.09254843059208,0,0];
Y_ijkl_2m1_vector= [0,0,0,0,0,0,0,-1.09254843059208,0,0,0,-1.09254843059208,0,0,-1.09254843059208];
Y_ijkl_20_vector = [-0.315391565252520,-0.315391565252520,0.630783130505040,-0.630783130505040,0.315391565252520,0,0,0,0.315391565252520,0,0,0,0,0,0];
Y_ijkl_21_vector = [0,0,0,0,0,0,-1.09254843059208,0,0,0,-1.09254843059208,0,0,-1.09254843059208,0];
Y_ijkl_22_vector = [0.546274215296039,-0.546274215296039,0,0,0.546274215296040,0,0,0,-0.546274215296040,0,0,0,0,0,0];
Y_ijkl_4m4_vector= [0,0,0,0,0,2.50334294179671,0,0,0,-2.50334294179671,0,0,0,0,0];
Y_ijkl_4m3_vector= [0,0,0,0,0,0,0,-5.31039230933979,0,0,0,1.77013076977993,0,0,0];
Y_ijkl_4m2_vector= [0,0,0,0,0,-0.946174695757560,0,0,0,-0.946174695757560,0,0,5.67704817454536,0,0];
Y_ijkl_4m1_vector= [0,0,0,0,0,0,0,2.00713963067187,0,0,0,2.00713963067187,0,0,-2.67618617422916];
Y_ijkl_40_vector = [0.317356640745613,0.317356640745613,0.846284375321635,0.634713281491226,-2.53885312596490,0,0,0,-2.53885312596490,0,0,0,0,0,0];
Y_ijkl_41_vector = [0,0,0,0,0,0,2.00713963067187,0,0,0,2.00713963067187,0,0,-2.67618617422916,0];
Y_ijkl_42_vector = [-0.473087347878780,0.473087347878780,0,0,2.83852408727268,0,0,0,-2.83852408727268,0,0,0,0,0,0];
Y_ijkl_43_vector = [0,0,0,0,0,0,-1.77013076977993,0,0,0,5.31039230933979,0,0,0,0];
Y_ijkl_44_vector = [0.625835735449177,0.625835735449177,0,-3.75501441269506,0,0,0,0,0,0,0,0,0,0,0];

Slm(1,:) =1/5*(1/C0^2)*Y_ijkl_00_vector*S_ijkl;
Slm(2,:) =4/7*(1/C2^2)*Y_ijkl_2m2_vector*S_ijkl;
Slm(3,:) =4/7*(1/C2^2)*Y_ijkl_2m1_vector*S_ijkl;
Slm(4,:) =4/7*(1/C2^2)*Y_ijkl_20_vector*S_ijkl;
Slm(5,:) =4/7*(1/C2^2)*Y_ijkl_21_vector*S_ijkl;
Slm(6,:) =4/7*(1/C2^2)*Y_ijkl_22_vector*S_ijkl;
Slm(7,:) =8/35*(1/C4^2)*Y_ijkl_4m4_vector*S_ijkl;
Slm(8,:) =8/35*(1/C4^2)*Y_ijkl_4m3_vector*S_ijkl;
Slm(9,:) =8/35*(1/C4^2)*Y_ijkl_4m2_vector*S_ijkl;
Slm(10,:)=8/35*(1/C4^2)*Y_ijkl_4m1_vector*S_ijkl;
Slm(11,:)=8/35*(1/C4^2)*Y_ijkl_40_vector*S_ijkl;
Slm(12,:)=8/35*(1/C4^2)*Y_ijkl_41_vector*S_ijkl;
Slm(13,:)=8/35*(1/C4^2)*Y_ijkl_42_vector*S_ijkl;
Slm(14,:)=8/35*(1/C4^2)*Y_ijkl_43_vector*S_ijkl;
Slm(15,:)=8/35*(1/C4^2)*Y_ijkl_44_vector*S_ijkl;

end

function [Sij] = STF2cart_rank2_tensor(S2m)
% Sij: 2D array containing [S_11;S_22;S_33;S_12;S_13;S_23]
C0=sqrt(1/(4*pi));
Y_ij_00=C0*eye(3);
Y_ij_2m2=[0,0.546274215296040,0;0.546274215296040,0,0;0,0,0];
Y_ij_2m1=[0,0,0;0,0,-0.546274215296040;0,-0.546274215296040,0];
Y_ij_20 =[-0.315391565252520,0,0;0,-0.315391565252520,0;0,0,0.630783130505040];
Y_ij_21 =[0,0,-0.546274215296040;0,0,0;-0.546274215296040,0,0];
Y_ij_22 =[0.546274215296040,0,0;0,-0.546274215296040,0;0,0,0];
Y_ij_all=cat(3,Y_ij_00,Y_ij_2m2,Y_ij_2m1,Y_ij_20,Y_ij_21,Y_ij_22);
idx=[1 1;2 2;3 3;1 2;1 3;2 3];
Sij=0*S2m;
for ii=1:6
    for jj=1:6
        Sij(ii,:)=Sij(ii,:)+S2m(jj,:)*Y_ij_all(idx(ii,1),idx(ii,2),jj);
    end
end
end

function [Sijkl] = STF2cart_rank4_tensor(S4m)
% Sijkl: 2D array containing [S_1111;S_2222;S_3333;S_1122;S_1133;S_1112;S_1113;S_1123;S_2233;S_2212;S_2213;S_2223;S_3312;S_3313;S_3323];   
load('/mnt/labspace/Santiago/MyRobustDKI/stuff/Y_ijkl_all.mat')
idx=[ 1,1,1,1; 2,2,2,2; 3,3,3,3; 1,1,2,2; 1,1,3,3; 1,1,1,2; 1,1,1,3; 1,1,2,3; 2,2,3,3; 2,2,1,2; 2,2,1,3; 2,2,2,3; 3,3,1,2; 3,3,1,3; 3,3,2,3];
Sijkl=0*S4m;
for ii=1:15
    for jj=1:15
        Sijkl(ii,:)=Sijkl(ii,:)+S4m(jj,:)*Y_ijkl_all(idx(ii,1),idx(ii,2),idx(ii,3),idx(ii,4),jj);
    end
end
% clc,clear,close all
% y=getY(4);
% Y_ijkl_all(:,:,:,:,1)=y{2,1,1};
% % y{2,1,1}(:)'*y{2,1,1}(:)*1/5*(1/C0^2)
% 1/5*y{2,1,1}(1,1,1,1)+y{2,1,1}(2,2,2,2)+y{2,1,1}(3,3,3,3)+2*(y{2,1,1}(1,1,2,2)+y{2,1,1}(1,1,3,3)+y{2,1,1}(2,2,3,3))/sqrt(4*pi)
% 
% Y_ijkl_all(:,:,:,:,2)=y{2,2,1};
% Y_ijkl_all(:,:,:,:,3)=y{2,2,2};
% Y_ijkl_all(:,:,:,:,4) =y{2,2,3};
% Y_ijkl_all(:,:,:,:,5) =y{2,2,4};
% Y_ijkl_all(:,:,:,:,6) =y{2,2,5};
% Y_ijkl_all(:,:,:,:,7)=y{2,3,1};
% Y_ijkl_all(:,:,:,:,8)=y{2,3,2};
% Y_ijkl_all(:,:,:,:,9)=y{2,3,3};
% Y_ijkl_all(:,:,:,:,10)=y{2,3,4};
% Y_ijkl_all(:,:,:,:,11) =y{2,3,5};
% Y_ijkl_all(:,:,:,:,12) =y{2,3,6};
% Y_ijkl_all(:,:,:,:,13) =y{2,3,7};
% Y_ijkl_all(:,:,:,:,14) =y{2,3,8};
% Y_ijkl_all(:,:,:,:,15) =y{2,3,9};
% save('/Volumes/labspace/Santiago/MyRobustDKI/Y_ijkl_all.mat','Y_ijkl_all')
end
function Ylm_n = get_even_SH(dirs,Lmax,CS_phase)
% Ylm_n = get_even_SH(dirs,Lmax,CS_phase)
%
% if CS_phase=1, then the definition uses the Condon-Shortley phase factor
% of (-1)^m. Default is CS_phase=0 (so this factor is ommited)
%
% By: Santiago Coelho (10/06/2021)


if size(dirs,2)~=3
    dirs=dirs';
end
Nmeas=size(dirs,1);
[PHI,THETA,~]=cart2sph(dirs(:,1),dirs(:,2),dirs(:,3)); THETA=pi/2-THETA;
l=0:2:Lmax;
l_all=[];
m_all=[];
for ii=1:length(l)
    l_all=[l_all, l(ii)*ones(1,2*l(ii)+1)];
    m_all=[m_all -l(ii):l(ii)];
end
K_lm=sqrt((2*l_all+1)./(4*pi) .* factorial(l_all-abs(m_all))./factorial(l_all+abs(m_all)));
if nargin==2 || isempty(CS_phase) || ~exist('CS_phase','var') || ~CS_phase
    extra_factor=ones(size(K_lm));
    extra_factor(m_all~=0)=sqrt(2);
else
    extra_factor=ones(size(K_lm));
    extra_factor(m_all~=0)=sqrt(2);
    extra_factor=extra_factor.*(-1).^(m_all);
end
P_l_in_cos_theta=zeros(length(l_all),Nmeas);
phi_term=zeros(length(l_all),Nmeas);
id_which_pl=zeros(1,length(l_all));
for ii=1:length(l_all)
    all_Pls=legendre(l_all(ii),cos(THETA));
    P_l_in_cos_theta(ii,:)=all_Pls(abs(m_all(ii))+1,:);
    id_which_pl(ii)=abs(m_all(ii))+1;
    if m_all(ii)>0
        phi_term(ii,:)=cos(m_all(ii)*PHI);
    elseif m_all(ii)==0
        phi_term(ii,:)=1;
    elseif m_all(ii)<0
        phi_term(ii,:)=sin(-m_all(ii)*PHI);
    end
end
Y_lm=repmat(extra_factor',1,Nmeas).*repmat(K_lm',1,Nmeas).*phi_term.*P_l_in_cos_theta;
Ylm_n=Y_lm';
end
