function [b0, dt] = dki_fit(dwi, grad, mask, constraints, outliers)
    % Diffusion Kurtosis Imaging tensor estimation using 
    % (constrained) weighted linear least squares estimation 
    % -----------------------------------------------------------------------------------
    % please cite:  Veraart, J.; Sijbers, J.; Sunaert, S.; Leemans, A. & Jeurissen, B.,
    %               Weighted linear least squares estimation of diffusion MRI parameters: 
    %               strengths, limitations, and pitfalls. NeuroImage, 2013, 81, 335-346 
    %------------------------------------------------------------------------------------
    % 
    % Usage:
    % ------
    % [b0, dt] = dki_fit(dwi, grad [, mask [, constraints]])
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
    % Optional input:
    % ---------------
    %    3. mask (boolean; [x, y, x]), providing a mask limits the
    %       calculation to a user-defined region-of-interest.
    %       default: mask = full FOV
    %
    %    4 . constraints (boolean; [1, 3] as in [c1, c2, c3]), imposes
    %       user-defined constraint to the weighted linear leasts squares
    %       estimation of the diffusion kurtosis tensor.
    %       Following constraints are available:
    %           c1: Dapp > 0
    %           c2: Kapp > 0
    %           c3: Kapp < b/(3*Dapp)
    %       default: [0 1 0]
    % 
    % Copyright: Vision Lab, University of Antwerp (2016)
    % 
    % This Source Code Form is subject to the terms of the Mozilla Public
    % License, v. 2.0. If a copy of the MPL was not distributed with this file,
    % You can obtain one at http://mozilla.org/MPL/2.0/
    % 
    % This code is distributed  WITHOUT ANY WARRANTY; without even the 
    % implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    % 
    % For more details, contact: Jelle.Veraart@nyumc.org 

    %% parameter checks
    dwi = double(dwi);
    dwi(dwi<=0)=eps;
    [x, y, z, ndwis] = size(dwi);
    if ~exist('grad','var') || size(grad,1) ~= ndwis || size(grad,2) ~= 4
        error('');
    end
    grad = double(grad);
    grad(:,1:3) = bsxfun(@rdivide,grad(:,1:3),sqrt(sum(grad(:,1:3).^2,2))); grad(isnan(grad)) = 0;

    if ~exist('mask','var') || isempty(mask)
        mask = true(x, y, z);
    end
    
    if ~exist('outliers', 'var') || isempty(outliers)
        outliers = false(size(dwi));
    end
    
    dwi = vectorize(dwi, mask);
    outliers = vectorize(outliers, mask);
   
    if exist('constraints', 'var') && ~isempty(constraints) && numel(constraints)==3
    else
        constraints = [0 1 0];
    end
    constraints = constraints > 0;
    %% tensor fit

    ind = [ 1     1;  1     2;   1     3;   2     2;   2     3;  3     3];
    cnt = [ 1     2     2     1     2     1 ];
    s = size(grad);
    b = [ones([s(1) 1],'like', grad) -(grad(:,ones(1,6)*4).*grad(:,ind(1:6,1)).*grad(:,ind(1:6,2)))*diag(cnt)];
    W_ind = [1 1 1 1; 1 1 1 2; 1 1 1 3; 1 1 2 2; 1 1 2 3;
        1 1 3 3; 1 2 2 2; 1 2 2 3; 1 2 3 3; 1 3 3 3;
        2 2 2 2; 2 2 2 3; 2 2 3 3; 2 3 3 3; 3 3 3 3];
    W_cnt = [1 4 4 6 12 6 4 12 12 4 1 4 6 4 1];
    bsqd6 = grad(:, 4).^2 * ones(1,15)/6;
    b = [b, (bsqd6 .* prod(reshape(grad(:,W_ind),[],15,4),3))*diag(W_cnt)];

    % unconstrained LLS fit
    dt = b\log(dwi);
    w = exp(b*dt);

    nvoxels = size(dwi,2);

    % WLLS fit initialized with LLS
    if any(constraints) 
        dir = [0.382517725304416 -0.748614094922528 0.541532202838631;-0.266039846728327 0.963894740823927 0.0113898448076587;-0.128563443377023 0.800867622029815 0.584878186472493;0.691696803043553 0.485345502199397 0.534785261721136;0.776929615593511 -0.627201085171846 0.0547646891069225;-0.314229418625565 0.891550800503996 0.326185595314880;-0.984699447847175 0.0338717154803320 0.170937720529697;0.729869942283584 0.134539815263771 0.670215566411097;0.0491118066650937 0.613801560175467 0.787931262974286;0.615167937666214 0.786762996419759 0.0507187926916626;-0.504930428375015 -0.548805916531712 0.666226184175323;0.514775318788445 0.353967263592948 0.780841563616317;-0.306616550550256 0.577152309970495 0.756889359169743;-0.644455563348338 0.445243323148325 0.621639292565402;0.888177219438464 0.244852048242751 0.388829913126405;-0.115867623474531 0.331617270421714 0.936271691224516;0.312724982544119 -0.262437525100548 0.912868901163732;-0.348318356641730 -0.328727572647744 0.877845376707953;0.622993255900061 -0.170127400464004 0.763502502100944;-0.870285082134136 -0.0832402149084147 0.485463636575162;0.879693901263504 -0.0847887289384472 0.467920411528283;0.375735817168569 0.624963320743740 0.684283160264539;-0.763508679267313 0.569075961898599 0.305298290648126;0.895786299773340 -0.371201461149296 0.244492086536589;0.431182280410191 0.0594580709470589 0.900303603713504;-0.927083085686508 -0.288337655567580 0.239537781186549;0.208899044398678 0.833216349905585 0.511968459477078;-0.671275756453876 -0.252498824452251 0.696873878436771;-0.385511621254227 -0.908766073079027 0.159765497834991;-0.501120479467597 0.703268192077924 0.504273849281930;-0.578440272143465 0.801933906922628 0.149361509400528;0.986601726072896 -0.0507533113495985 0.155052041254001;0.0472262384668294 -0.790665651327184 0.610424041311968;0.957038035873056 0.279601625450131 0.0768188058868917;-0.497573767044291 -0.0342706449545790 0.866744408256408;0.537095370960702 0.746985750118871 0.391842891490881;0.174500355902118 -0.559258086805823 0.810419655568845;-0.0648836431571087 -0.997212186937296 0.0368506048036694;-0.200896533381969 0.00230971954655800 0.979609742739793;-0.436037609875685 0.290696319170598 0.851684714430502;0.332217034685261 0.924381756972555 0.187483890618005;0.115538097684954 0.0265470728743124 0.992948236770250;-0.167448247267712 -0.594347070791611 0.786582890691377;-0.931940478288593 0.352679013976080 0.0842879471104161;0.749660835628331 -0.375215305423717 0.545180801295148;-0.112213298457421 -0.929475988744578 0.351400856567817;-0.596160909541517 -0.730179079768923 0.333812344592648;0.211077351410955 0.350067854669890 0.912632921194583;-0.325168748302559 -0.780267672863407 0.534273004943793;-0.717210613875971 0.128994202815781 0.684813427864535;-0.0218381924490005 -0.303713916706869 0.952512965869302;0.213275433291729 -0.924500457406975 0.315931153589701;-0.810453788321924 -0.574858547954973 0.112704511168551;0.665549405791414 -0.637998127347686 0.387301404530814;0.489520321770316 -0.495410872500818 0.717591751612200;0.514060443295042 -0.837385561080287 0.185815184346054;-0.757892441488769 -0.466692556842526 0.455847676885577;0.00471100435105065 0.958734616992657 0.284263505603424;0.800137357904460 0.555340864988139 0.226664360144898;-0.872328992553570 0.265326196661285 0.410663046932313;];
        C = [];
        if constraints(1)>0
            C = [C; [zeros(60, 1), (dir(:,ind(:,1)).*dir(:,ind(:,2)))*diag(cnt), zeros(60, 15)]];
        end
        if constraints(2)>0
            C = [C; [zeros(60, 7), ( prod(reshape(dir(:,W_ind),[],15,4),3))*diag(W_cnt)]];
        end
        if constraints(3)>0
            C = [C; [zeros(60, 1), 3/max(bval)*(dir(:,ind(1:6,1)).*dir(:,ind(1:6,2)))*diag(cnt), -(prod(reshape(dir(:,W_ind),[],15,4),3))*diag(W_cnt)]];
        end
        d = zeros([1, size(C, 1)]);
        options = optimset('LargeScale', 'off', 'Display', 'off', 'MaxIter', 22000, 'TolCon', 1e-12, 'TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 220000);
        parfor i = 1:nvoxels
            in_ = outliers(:, i) == 0;
            wi = w(:,i); Wi = diag(wi(in_));             
            dt(:, i) = lsqlin(Wi*b(in_, :),Wi*log(dwi(in_,i)),-C, d, [],[],[],[],[],options);
        end
    else
        parfor i = 1:nvoxels
            in_ = outliers(:, i) == 0;
            b_ = b(in_, :);
            if isempty(b_) || cond(b(in_, :))>1e15
                dt(:, i) = NaN
            else
                wi = w(:,i); Wi = diag(wi(in_)); 
                logdwii = log(dwi(in_,i));
                dt(:,i) = (Wi*b_)\(Wi*logdwii);
            end
        end
    end

    b0 = exp(dt(1,:));
    dt = dt(2:end, :);
    D_apprSq = 1./(sum(dt([1 4 6],:),1)/3).^2;
    dt(7:21,:) = dt(7:21,:) .* D_apprSq(ones(15,1),:);
    b0 = vectorize(b0, mask);
    dt = vectorize(dt, mask);
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
