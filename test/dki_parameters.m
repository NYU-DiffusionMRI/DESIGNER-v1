function [fa, md, rd, ad, fe, mk,  rk, ak, kfa, mkt] = dki_parameters(dt, mask)
    % diffusion and kurtosis tensor parameter calculation
    %
    % -----------------------------------------------------------------------------------
    % please cite: Veraart et al.  
    %              More Accurate Estimation of Diffusion Tensor Parameters Using Diffusion Kurtosis Imaging,
    %              MRM 65 (2011): 138-145.
    %------------------------------------------------------------------------------------
    % 
    % Usage:
    % ------
    % [fa, md, ad, rd, fe, mk, ak, rk] = dki_parameters(dt [, mask [, branch]])
    %
    % Required input: 
    % ---------------
    %     1. dt: diffusion kurtosis tensor (cf. order of tensor elements cf. dki_fit.m)
    %           [x, y, z, 21]
    %
    % Optional input:
    % ---------------
    %    2. mask (boolean; [x, y, x]), providing a mask limits the
    %       calculation to a user-defined region-of-interest.
    %       default: mask = full FOV
    %
    %    3. branch selection, 1 or 2 (default: 1)
    %              1. De_parallel > Da_parallel
    %              2. Da_parallel > De_parallel
    % 
    % output:
    % -------
    %  1. fa:                fractional anisitropy
    %  2. md:                mean diffusivity
    %  3. rd:                radial diffusivity
    %  4. ad:                axial diffusivity
    %  5. fe:                principal direction of diffusivity
    %  6. mk:                mean kurtosis
    %  7. rk:                radial kurtosis
    %  8. ak:                axial kurtosis
    %
    % Important: The presence of outliers "black voxels" in the kurtosis maps 
    %            are we well-known, but inherent problem to DKI. Smoothing the 
    %            data in addition to the typical data preprocessing steps
    %            might minimize the impact of those voxels on the visual
    %            and statistical interpretation. However, smoothing comes
    %            with the cost of partial voluming. 
    %       
    % Copyright (c) 2017 New York University and University of Antwerp
    % 
    % This Source Code Form is subject to the terms of the Mozilla Public
    % License, v. 2.0. If a copy of the MPL was not distributed with this file,
    % You can obtain one at http://mozilla.org/MPL/2.0/
    % 
    % This code is distributed  WITHOUT ANY WARRANTY; without even the 
    % implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    % 
    % For more details, contact: Jelle.Veraart@nyumc.org 
    
    
n = size(dt, 4);
if ndims(dt)~=4
    error('size of dt needs to be [x, y, z, 21]')
end
if n~=21
    error('dt needs to contain 21')
end
if ~exist('mask','var') || isempty(mask)     
        mask = ~isnan(dt(:,:,:,1));
end
    
dt = vectorize(dt, mask);
nvoxels = size(dt, 2);


%% DTI parameters
for i = 1:nvoxels
    DT = dt([1:3 2 4 5 3 5 6], i);
    DT = reshape(DT, [3 3]);
    try
        [eigvec, eigval] = eigs(DT);
        eigval = diag(eigval);
    catch
        eigvec = NaN(3, 3);
        eigval = NaN(3, 1);
    end
    [eigval, idx] = sort(eigval, 'descend');
    eigvec = eigvec(:, idx);
    l1(i) = eigval(1,:);
    l2(i) = eigval(2,:);
    l3(i) = eigval(3,:);
    
    e1(:, i) = eigvec(:, 1); 
end
md = (l1+l2+l3)/3;
rd = (l2+l3)/2;
ad = l1;
fa = sqrt(1/2).*sqrt((l1-l2).^2+(l2-l3).^2+(l3-l1).^2)./sqrt(l1.^2+l2.^2+l3.^2);

%% DKI parameters
dirs = get256dirs();
akc = AKC(dt, dirs);

mk = mean(akc);
ak = zeros([1, size(e1,2)]);
rk = zeros([1, size(e1,2)]);

parfor i = 1:nvoxels
    dirs = [e1(:,i), -e1(:,i)]';
    akc = AKC(dt(:,i), dirs);
    ak(i) = mean(akc);
    dirs = radialsampling(e1(:,i), 256)';
    akc = AKC(dt(:,i), dirs);
    rk(i) = mean(akc);
    [kfa(i),mkt(i)] = ComputeKFA(dt(:,i),0,3)
end
               
%% return maps
fa  = vectorize(fa, mask);
md  = vectorize(md, mask);
ad  = vectorize(ad, mask);
rd  = vectorize(rd, mask);
mk  = vectorize(mk, mask);
ak  = vectorize(ak, mask);
rk  = vectorize(rk, mask);
fe  = vectorize(e1, mask);
kfa = vectorize(kfa, mask);
mkt = vectorize(mkt, mask);
end
