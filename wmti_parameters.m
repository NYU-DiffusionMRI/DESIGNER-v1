function [awf, eas, ias] = wmti_parameters(dt, mask, branch)
    % Calculation of the white matter tract integrity metrics
    %
    % -----------------------------------------------------------------------------------
    % please cite: Fieremans, Els, Jens H. Jensen, and Joseph A. Helpern. 
    %              White matter characterization with diffusional kurtosis imaging,
    %              Neuroimage 58.1 (2011): 177-188.
    %------------------------------------------------------------------------------------
    % 
    % Usage:
    % ------
    % [awf, eas, ias] = wmti_parameters(dt [, mask [, branch]])
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
    %  1. awf:               axonal water fraction
    %  2. eas:              extra-axonal diffusivities, incl. tortuosity
    %  3. ias:              intra-axonal diffusivities
    %
    % Important: The presence of outliers "black voxels" in the kurtosis maps 
    %            are we well-known, but inherent problem to DKI. Smoothing the 
    %            data in addition to the typical data preprocessing steps
    %            might minimize the impact of those voxels on the visual
    %            and statistical interpretation. However, smoothing comes
    %            with the cost of partial voluming. 
    %       
    % Copyright (c) 2016 New York University and University of Antwerp
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
if ~exist('branch','var') || isempty(branch)     
        branch = 1;
end
    
dt = vectorize(dt, mask);
nvoxels = size(dt, 2);


%% DTI parameters
for i = 1:nvoxels
    DT = dt([1:3 2 4 5 3 5 6], i);
    DT = reshape(DT, [3 3]);
    [eigvec, eigval] = eigs(DT);
    eigval = diag(eigval);
    
    l1(i) = eigval(1,:);
    l2(i) = eigval(2,:);
    l3(i) = eigval(3,:);
    
    e1(:, i) = eigvec(:, 1); 
    e2(:, i) = eigvec(:, 2); 
    e3(:, i) = eigvec(:, 3);
end


%% WMTI paramerers
% creation of 10000 directions, regularly distributed on half a sphere.
load('dirs10000.mat', 'dir')

% axonal water fraction
nvxls = size(dt, 2);
maxk = zeros([1 nvxls]);
N = 10000;
nblocks = 10;
for i=1:nblocks;
    akc = AKC(dt, dir((N/nblocks*(i-1))+1:N/nblocks*i, :)); %#ok<NODEF>
    maxk = nanmax([maxk; akc], [], 1);
end

clear dir
awf = maxk./(maxk+3); awf(isinf(awf))=0;


for i=1:nvxls
    akc = AKC(dt(:,i), [e1(:,i), e2(:,i), e3(:,i)]);

    if branch == 1
        eas.de1(i) = l1(i) * (1 + sqrt((akc(1)*awf(i))/(3*(1-awf(i)))));
        ias.da1(i) = l1(i) * (1 - sqrt((akc(1)*(1-awf(i)))./(3*awf(i))));
    else
        eas.de1(i) = l1(i) * (1 - sqrt((akc(1)*awf(i))/(3*(1-awf(i)))));
        ias.da1(i) = l1(i) * (1 + sqrt((akc(1)*(1-awf(i)))./(3*awf(i))));
    end
    eas.de2(i) = l2(i) * (1 + sqrt((akc(2)*awf(i))/(3*(1-awf(i)))));
    ias.da2(i) = l2(i) * (1 - sqrt((akc(2)*(1-awf(i)))./(3*awf(i))));

    eas.de3(i) = l3(i) * (1 + sqrt((akc(3)*awf(i))/(3*(1-awf(i)))));
    ias.da3(i) = l3(i) * (1 - sqrt((akc(3)*(1-awf(i)))./(3*awf(i))));

    eas.tort(i) = eas.de1(i)./(.5*(eas.de2(i) + eas.de3(i)));
    eas.de_perp(i) = (eas.de2(i) + eas.de3(i))/2;

    ias.da_perp(i) = (ias.da2(i) + ias.da3(i))/2;
end


            
            
%% return maps
awf = vectorize(awf, mask);

fields = fieldnames(ias);
for i=1:numel(fields)
    ias = setfield(ias, fields{i}, vectorize(getfield(ias, fields{i}), mask));  %#ok<*GFLD>
end

fields = fieldnames(eas);
for i=1:numel(fields)
    eas = setfield(eas, fields{i}, vectorize(getfield(eas, fields{i}), mask));  %#ok<*SFLD>
end
    
    
end


function [akc, adc] = AKC(dt, dir)

W_ind = [1 1 1 1; 1 1 1 2; 1 1 1 3; 1 1 2 2; 1 1 2 3;
    1 1 3 3; 1 2 2 2; 1 2 2 3; 1 2 3 3; 1 3 3 3;
    2 2 2 2; 2 2 2 3; 2 2 3 3; 2 3 3 3; 3 3 3 3];
W_cnt = [1 4 4 6 12 6 4 12 12 4 1 4 6 4 1];

adc = ADC(dt(1:6, :), dir);
md = sum(dt([1 4 6],:),1)/3;
T = (prod(reshape(dir(:,W_ind),[],15,4),3))*diag(W_cnt);
akc =  T*dt(7:21, :);
akc = (akc .* repmat(md.^2, [size(adc, 1), 1]))./(adc.^2);
end

function [adc] = ADC(dt, dir)
ind = [ 1     1;  1     2;   1     3;   2     2;   2     3;  3     3];
cnt = [ 1     2     2     1     2     1 ];
adc = (dir(:,ind(1:6,1)).*dir(:,ind(1:6,2))) * diag(cnt) * dt;
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