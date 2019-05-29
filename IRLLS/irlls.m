function [ reject, dt, conv, fa, md ] = irlls( dwi, mask, g, b, sigma, options )
%IRLLS This functions performs outlier detection and robust parameter 
% estimation for diffusion MRI using the iterative reweigthed linear
% least squares (IRLLS) approach.
%
%   INPUT VARIABLES*
%       dwi:     [ndwi x nvox] or [nvoxX x nvoxY x nvoxZ x ndwi] diffusion 
%                 weighted images.
%       grad:    [ndwi x 3] matrix of gradient directions.
%       bvals:   [ndwi x 1] vector of b-values.
%       sigma:   [nvox x 1] estimation of the noise level.  --> adjusted to
%                incorporate heterogeneous noise maps
%       options: [structure] see OPTIONS
%
%   OUTPUT VARIABLES*
%       reject:  [ndwi x nvox] outlier map.
%       dt:      [nparam x nvox] diffusion tensor.
%       fa:      [nvox x 1] fractional anisotropy.
%       md:      [nvox x 1] mean diffusivity.
%       conv:    [nvox x 1] map of number of iterations needed to reach
%                 convergence.
%
%   OPTIONS
%       Excludeb0: Exlude the b0 images when removing outliers?
%           0: no
%           1: yes (default)
%       MaxIter: Maximum number of iterations in the iterative reweighting 
%        loop.
%           default: 25
%       ConvCrit: Fraction of L2-norm of estimated diffusion parameter 
%        vector that the L2-norm of the difference vector should get under 
%        in order to reach convergence in the iterative reweighted loop.
%           default: 1e-3
%       Kurtosis: Use DTI or DKI model
%           0: DTI (default)
%           1: DKI
%       Leverage: Measurements with a leverage above this threshold will 
%        not be excluded after outlier detection.
%           default: 0.85
%       Bounds: Set the threshold of the number of standard deviations that
%        are needed to exclude a measurement.
%           default: 3
%
%   *
%       ndwi = number of diffusion weighted images/volumes.
%       nvox = number of voxels in each masked volume.
%       nvoxX/nvoxY/nvoxZ = number of voxels in X/Y/Z direction in a single
%        volume.
%       nparam = number of parameters that need to be estimated.
%           DTI: nparam = 7  
%           DKI: nparam = 22
% 
%  ------------------------------------------------------------------------
%   Created: April 4th 2014
%   Last edit: July 2016
% 
%   Copyright 2014 University of Antwerp, Antwerp, Belgium
% 
%   Written by Quinten Collier (quinten.collier@uantwerpen.be), Jelle 
%   Veraart and Ben Jeurissen, 2014-2015
% 
%   This program is software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with IRLLS.m.  If not, see <http://www.gnu.org/licenses/>.
% 
%   When using this code, please refer to: 
%   Collier, Q., Veraart, J., Jeurissen, B., den Dekker, A. J., 
%   & Sijbers, J. (2015). Iterative reweighted linear least squares 
%   for accurate, fast, and robust estimation of diffusion magnetic 
%   resonance parameters. Magnetic resonance in medicine, 73(6), 
%   2174-2184. DOI: 10.1002/mrm.25351
%--------------------------------------------------------------------------


%% check options
    if ~isfield(options, 'Excludeb0')
        options.Excludeb0 = 1;
    elseif ~(options.Excludeb0==0 || options.Excludeb0==1)
        error('option: Excludeb0 should be set to 0 or 1');
    end
    
    if ~isfield(options, 'MaxIter')
        options.MaxIter = 25;
    elseif options.MaxIter<1 || options.MaxIter>200
        error('option: Maxiter should be set to a value between 1 and 200');
    end
    
    if ~isfield(options,'ConvCrit')
        options.ConvCrit = 1e-3;
    elseif options.ConvCrit<0 || options.ConvCrit>1
        error('option: ConvCrit should be set to a value between 0 and 1');
    end
    
    if ~isfield(options,'Kurtosis')
        options.Kurtosis = 0;
    elseif ~(options.Kurtosis==0 || options.Kurtosis==1)
        error('option: Kurtosis should be set to 0 or 1');
    end
    
    if ~isfield(options,'Leverage')
        options.Leverage = 0.85;
    elseif options.Leverage<0 || options.Leverage>1
        error('option: Leverage should be set to a value between 0 and 1');
    end

    if ~isfield(options,'Bounds')
        options.Bounds = 3;
    elseif options.Bounds<1
        error('option: Bounds should be set to a value >= 1');
    end
    
    names = fieldnames(options);
    tmp = ismember(names,{'Excludeb0','MaxIter','ConvCrit','Kurtosis','Leverage','Bounds'});
    for i = 1:numel(tmp)
        if tmp(i)==0
           error([names{i}, ' is not a valid option']); 
        end
    end
    clear tmp names
    
%% Vectorization of dwi
    dwi = double(dwi);
    v = false;
    if ~ismatrix(dwi)
        v = true;
        if ~exist('mask','var') || isempty(mask)
            mask = ~isnan(dwi(:,:,:,1));
        end
        dwi_ = zeros([size(dwi,4) sum(mask(:))], class(dwi));
        for k = 1:size(dwi,4);
            tmp = dwi(:,:,:,k);
            dwi_(k,:) = tmp(mask(:));
        end
        dwi = dwi_;
        clear dwi_ tmp
    end
    [ndwi, nvox] = size(dwi);
    
%% scaling
    scaling = false;
    if numel(dwi(dwi<1))/numel(dwi) < 0.001
        dwi(dwi<1)=1;
    else
        scaling = true;
        if max(b)<10
            tmp = dwi(b<0.05,:);
        else
            tmp = dwi(b<50,:);
        end
        sc = median(tmp(:));
        dwi(dwi<sc/1000) = sc/1000;
        dwi = dwi*1000/sc;
    end

%% Create B-matrix
    if options.Kurtosis
        D_ind = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
        D_cnt = [1 2 2 1 2 1];
        W_ind = [1 1 1 1; 1 1 1 2; 1 1 1 3; 1 1 2 2; 1 1 2 3;
                 1 1 3 3; 1 2 2 2; 1 2 2 3; 1 2 3 3; 1 3 3 3;
                 2 2 2 2; 2 2 2 3; 2 2 3 3; 2 3 3 3; 3 3 3 3];
        W_cnt = [1 4 4 6 12 6 4 12 12 4 1 4 6 4 1];
        bmat = [ones([ndwi,1]), -repmat(b,[1,6]).*g(:,D_ind(:,1)).*g(:,D_ind(:,2))*diag(D_cnt), (1/6)*repmat(b,[1,15]).^2.*g(:,W_ind(:,1)).*g(:,W_ind(:,2)).*g(:,W_ind(:,3)).*g(:,W_ind(:,4))*diag(W_cnt)];
    else
        ind = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
        cnt = [1 2 2 1 2 1];
        bmat = [ones([ndwi,1]), -repmat(b,[1,6]).*g(:,ind(:,1)).*g(:,ind(:,2))*diag(cnt)];
    end
    nparam = size(bmat,2);
    ndof = ndwi - nparam;

%% Initialization
    b0_pos = false(size(b));
    if options.Excludeb0
        if max(b) < 10
            b0_pos = b < 0.01;            
        else
            b0_pos = b < 10;
        end
    end

    reject = false(size(dwi));
    conv = zeros(1,nvox);
    dt = zeros(nparam,nvox);
    fa = zeros(nvox,1);
    md = zeros(nvox,1);

%% Basic noise estimation
    if ~exist('sigma', 'var') || isempty(sigma)
        sigma_ = zeros([nvox,1]);
        parfor i=1:nvox
            dt_ = bmat\log(dwi(:,i));
            w = exp(bmat*dt_);
            dt_ = (bmat.*repmat(w,[1,nparam]))\(log(dwi(:,i)).*w);
            e = log(dwi(:,i)) - bmat*dt_;
            m = median(abs(e.*w - median(e.*w)));
            sigma_(i) = sqrt(ndwi/(ndof))*1.4826*m;
        end
        sigma = median(sigma_, 1);
        sigma = repmat(sigma,nvox,1);        
    elseif scaling        
        sigma = sigma*1000/sc;
    end


    parfor i = 1:nvox

%% preliminary rough outlier check
        dwi_i = dwi(:,i);
        dwi0 = median(dwi_i(b < 0.01));
        out = dwi_i > dwi0 + 3 * sigma(i);
        if sum(~out(b > 0.01))<size(bmat,2)-1
            out = false(size(out));
        end   
        out(b0_pos) = false;     
        bmat_i = bmat(~out,:);
        dwi_i = dwi_i(~out);
        n_i = numel(dwi_i);
        ndof_i = n_i - size(bmat_i,2);
        
%% WLLS estimation
        dt_i = bmat_i\log(dwi_i);
        w = exp(bmat_i*dt_i);
        dt_i = (bmat_i.*repmat(w,[1,nparam]))\(log(dwi_i).*w);
        dwi_hat = exp(bmat_i*dt_i);

%% Goodness-of-fit
        residu = log(dwi_i) - log(dwi_hat);
        residu_ = dwi_i - dwi_hat;
        chi2 = sum(residu_.*residu_./sigma(i).^2)/(ndof_i) - 1;                
        gof = abs(chi2) < 3*sqrt(2/(ndof_i));
        gof2 = gof;

%% Iterative reweighting procedure
        iter = 0;
        while ~gof && iter < options.MaxIter
             C = sqrt(n_i/(n_i-nparam))*1.4826*median(abs(residu_(:) - median(residu_(:))))./dwi_hat;
             GMM = C.^2./(residu.^2 + C.^2).^2;
             w = sqrt(GMM).*dwi_hat;
             dt_imin1 = dt_i;
             dt_i = (bmat_i.*repmat(w,[1,nparam]))\(log(dwi_i).*w);
             dwi_hat = exp(bmat_i*dt_i);
             dwi_hat(dwi_hat<1)=1;
             residu = log(dwi_i) - log(dwi_hat);
             residu_ = dwi_i - dwi_hat;

             % convergence check
             iter = iter+1;
             gof = norm(dt_i - dt_imin1) < norm(dt_i)*options.ConvCrit;
        end
        conv(i) = iter;

%% outlier detection
        if ~gof2
            leverage = diag(bmat_i* ((bmat_i'*diag(w.^2)*bmat_i)\(bmat_i'*diag(w.^2))) );
            lowerbound_linear = -options.Bounds*sqrt(1-leverage)*sigma(i)./dwi_hat;
            upperbound_nonlinear = options.Bounds*sqrt(1-leverage)*sigma(i);

            tmp = false(size(residu));
            tmp(residu < lowerbound_linear) = true;
            tmp(residu_ > upperbound_nonlinear) = true;
            tmp(leverage > options.Leverage) = false;
            tmp2 = true(size(b));
            tmp2(~out) = tmp;
            tmp2(b0_pos) = false;
            reject(:,i) = tmp2;
        else
            tmp2 = false(size(b));
            tmp2(out) = true;
            reject(:, i) = tmp2;
        end

%% Robust parameter estimation
        keep = ~reject(:, i);
        bmat_i = bmat(keep, :);
        dwi_i = dwi(keep,i);
        dt_ = bmat_i\log(dwi_i);
        w = exp(bmat_i*dt_);
        dt(:,i) = (bmat_i.*repmat(w,[1,nparam]))\(log(dwi_i).*w);
        dt_tmp = dt(:,i);
        dt2 = [dt_tmp(2), dt_tmp(3)/2, dt_tmp(4); dt_tmp(3)/2, dt_tmp(5), dt_tmp(6)/2; dt_tmp(4)/2, dt_tmp(6)/2, dt_tmp(7)];
        eigv = eig(dt2);
        fa(i) = sqrt(1/2)*sqrt( (eigv(1)-eigv(2))^2 + (eigv(1)-eigv(3))^2 + (eigv(2)-eigv(3))^2 ) / sqrt( eigv(1)^2 + eigv(2)^2 + eigv(3)^2 );
        md(i) = sum(eigv)/3;
    end
    
%% unscaling
    if scaling
        dt(1,:) = dt(1,:)+log(sc/1000);
    end

%% Unvectorizing output variables
    if v
       dims = size(mask);
       
       % unvec fa, md and conv
       fa_ = NaN([dims 1], 'double'); md_ = fa_; conv_ = fa_;
       fa_(mask) = fa; fa = fa_; clear fa_;
       md_(mask) = md; md = md_; clear md_;
       conv_(mask) = conv; conv = conv_; clear conv_;
       
       % unvec dt
       dt_ = NaN([dims size(dt,1)], class(dt));
       for k = 1:size(dt,1)
           tmp = NaN(dims, class(dt));
           tmp(mask) = dt(k,:);
           dt_(:,:,:,k) = tmp;
       end
       dt = dt_; clear dt_;
       
       % unvec reject
       reject_ = false([dims size(reject,1)]);
       for k = 1:size(reject,1)
           tmp = false(dims);
           tmp(mask) = reject(k,:);
           reject_(:,:,:,k) = tmp;
       end
       reject = reject_; clear reject_;       
    end

end



