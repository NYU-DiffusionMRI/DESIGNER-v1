function [dt, b0, mkpred, list] = RobustDKIFitting(dwi, grad, mask)

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
    
    extraRobustness = true;
    
    %% Step 1: prepare data
    % We (a) do some quick input checks, (b) vectorize the data to avoid handling with 4D data
    % structures and nested loops, and (c) build the DKI-specific b-matrix.  
    
    dwi = double(dwi);
    dwi(dwi<=0)=eps;
    scaleFact = 1000/max(dwi(:));
    dwi = dwi*scaleFact;
    
    
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
    options = optimset('lsqnonlin'); 
    options = optimset(options,'Jacobian','on','TolFun',1e-8,'TolX', ...
                        1e-8,'MaxIter',10000,'Display','off');
    
     % this is still very slow. It would be possible to switch everything
     % to WLLS to reduce the computation time. 
     for i = 1:nvox
         nls(:,i) = lsqnonlin(@(x)ObjF(x,double(bmat),double(y(:, i))),double(start(:,i)),[],[],options);
     end
     clear start

     y_hat = exp(bmat*nls);
     dt.nls = W2K(nls(2:22, :));
    

     

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
    
    
    %% Step 5: Voxel Quality Transfer
    % We train a mapping between the powder kurtosis and the mean kurtosis from the good voxels ... and apply them to the bad.  
     
    psi = 2/5 * Dsum./(md.nls.^2) - 6/5;
    idx = list.accept & ~list.IK.negative & psi<1;

    err = (mk.nls(idx)-iK(idx)) ./ mk.nls(idx);
    bnds = prctile(err, [1, 99]);
    err = (mk.nls-iK) ./ mk.nls;
    
    idx = list.accept & ~list.IK.negative & psi<1 & err<bnds(2) & err>bnds(1);
    
    [P.training, P.combs] = createDesignMatrix([iK(1, idx)', md.nls(1, idx)', Dsum(1, idx)']);
    P.coef = pinv(P.training)*mk.nls(1, idx)';
    
    %Predict mean kurtosis
    X = createDesignMatrix([iK(1, :)', md.nls(1, :)', Dsum(1, :)']);
    mk.predicted = (X*P.coef)';
    mk.predicted(isnan(mk.predicted)) = 0;

    
    
    
    %% Step 6 Compute Regularization term
    % Tuning the "alpha" parameter is a challenge for most regularized estimators. That's not different here. 
    mse.nls = nanmedian(sum((y_hat-y).^2)/size(y, 1));
    mse.mk = nanmedian((mk.predicted(idx) - mk.nls(idx)).^2);
    alpha = 0.1 *mse.nls ./ mse.mk; % can we increase?
    
    %% Step 7: Regularized fit with constrained fit as a starting point.
    %  regularized 
    [C, d] = createConstraints();
    opts = optimset('Display', 'off', 'Algorithm', 'interior-point', 'MaxIter', 10000, 'TolCon', 1e-8, 'TolFun', 1e-8, 'TolX', 1e-8);
  
    
    options = optimset('fminunc'); 
    options = optimset(options,'Jacobian','on','TolFun',1e-8,'TolX', ...
                        1e-8,'MaxIter',10000,'Display','off');
    
                      
    
    reg = nls;
    idx = find(list.remove);
    for i = 1:numel(idx)   
        try
        start(:, idx(i)) = lsqlin(double(bmat),double((log(y(:, idx(i))))),-C, d, [],[],[],[],[],opts);
        reg(:,idx(i)) = fminunc(@(x)regObjF(x,double(bmat),double(y(:, idx(i))), mk.predicted(idx(i)), Td, Tk, alpha),double(start(:,idx(i))),options);
        catch
            reg(:,idx(i)) = NaN;
        end
    end
    dt.reg = W2K(reg(2:22, :));
    b0.reg = exp(min(15,reg(1, :)))/scaleFact;
    %% Step 8: Quality Control
    [akc, Td, Tk] = AKC(dt.reg);    
    list.notSolvedAfterFirstTry = mean(akc, 1)<0;    
    %% Step 8: Update the remaining black voxels with stronger regularization and constrained fit
    if extraRobustness 
        
        options = optimset('fmincon'); 
        options = optimset(options,'Jacobian','on','TolFun',1e-8,'TolX', ...
                        1e-8,'MaxIter',10000,'Display','off');
          
        idx = find(list.notSolvedAfterFirstTry);
        for  i = 1:numel(idx)   
            try
                reg(:,idx(i)) = fmincon(@(x)regObjF(x,double(bmat),double(y(:, idx(i))), mk.predicted(idx(i)), Td, Tk, 10*alpha),double(start(:,idx(i))), -C, d, [],[], [], [], [], options);
            catch
                reg(:,idx(i)) = NaN;
            end
        end
        dt.reg = W2K(reg(2:22, :));
        b0.reg = exp(min(15,reg(1, :)))/scaleFact;

    end
    
    %% Step 9: write output
    dt = vectorize(dt.reg, mask);
    b0 = vectorize(b0.reg, mask);
    mkpred = vectorize(mk.predicted, mask);
    
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

function [F, J] = regObjF2(x, bmat, y, mkpred, Td, Tk, alpha)
    
        % data fidelity term
        S = exp(bmat*x);
        ny = numel(y);
        F1 = sum((y - S).^2)/ny;                %[1x1]
        
        tmp = -2*(y-S).*S;
        J1 = sum(tmp(:, ones(1, 22)).*bmat)/ny; %[1x22]
       
        
        % regularization term - compute "MK"
        adc = Td*x(2:7);   
        awc = Tk*x(8:22);  
        mk = sum(awc ./ adc.^2)/256;
        F2 = (mkpred-mk).^2;  

          
        % derivative of MK to x(1)
        dmk_dx(1) = 0;
        
        % derivative of MK to x(2:7)           
        dmk_dx(2:7) = 1/256*sum(-2*awc(:, [1 1 1 1 1 1]).*Td ./ (adc(:,[1 1 1 1 1 1]).^3));
        
        % derivative of MK to x(8:22)
        dmk_dx(8:22) = 1/256*sum(1./(adc(:, [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]).^2) .* Tk);
       
        J2 =  -2*(mkpred-mk)* dmk_dx; 
        
        F = F1 + alpha*F2;
        J = J1 + alpha*J2;
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
        