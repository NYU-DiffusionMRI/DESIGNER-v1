function [awf, eas, ias] = wmti_parameters_savedi(dt, mask)
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
    %   
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


dir =[-0.208098000000000,0.525514000000000,0.850005000000000;0.202387000000000,0.526131000000000,0.851002000000000;0.409956000000000,0.175267000000000,0.918257000000000;-0.412630000000000,0.742620000000000,0.565889000000000;-0.207127000000000,0.959492000000000,0.280092000000000;-0.872713000000000,0.525505000000000,0.0647640000000000;-0.746815000000000,0.526129000000000,0.455449000000000;-0.415238000000000,0.175473000000000,0.915841000000000;-0.746636000000000,0.175268000000000,0.673642000000000;0.665701000000000,0.742619000000000,-0.217574000000000;-0.330391000000000,0.959489000000000,-0.110458000000000;-0.331275000000000,0.525513000000000,-0.809983000000000;-0.663936000000000,0.526130000000000,-0.569521000000000;-0.999332000000000,0.175474000000000,-0.111904000000000;-0.871398000000000,0.175267000000000,-0.501922000000000;0.00121400000000000,0.742616000000000,-0.700356000000000;0.00294900000000000,0.959483000000000,-0.348370000000000;0.667975000000000,0.525509000000000,-0.565356000000000;0.336490000000000,0.526126000000000,-0.807431000000000;0.202383000000000,-0.175470000000000,0.985002000000000;0.208094000000000,0.175265000000000,-0.983848000000000;0.666452000000000,0.742619000000000,-0.215262000000000;0.332212000000000,0.959489000000000,-0.104850000000000;0.205064000000000,0.958364000000000,0.285421000000000;0.412630000000000,0.742620000000000,0.565889000000000;0.746093000000000,0.175315000000000,0.674232000000000;0.744110000000000,0.525505000000000,0.460568000000000;0.871894000000000,0.526125000000000,0.0705070000000000;0.874264000000000,0.175471000000000,-0.496841000000000;1,0.175267000000000,-0.106112000000000;];
[D_ind, D_cnt] = createTensorOrder(2);
adc2dt = pinv((dir(:,D_ind(1:6,1)).*dir(:,D_ind(1:6,2))) * diag(D_cnt));

for i=1:nvxls
       % try
            [akc, adc] = AKC(dt(:, i), dir); 
            akc(akc(:)<0) = 0; % avoiding complext output. However, negative AKC might be taken care of by applying constraints.
            %f = repmat(awf(i), [size(dir, 1), 1]);

            De=adc.*(1+sqrt(akc*awf(i)/(3*(1-awf(i)))));
            Di=adc.*(1-sqrt(akc*(1-awf(i))/(3*awf(i))));

            dt_e = adc2dt*De; 
            dt_i = adc2dt*Di; 
            
            % eigenvalue decomposition of De
            DTe = dt_e([1:3 2 4 5 3 5 6]);
            DTe = reshape(DTe, [3 3]);
            [~, eigval] = eigs(DTe); 
            eigval = sort(diag(eigval), 'descend');

            eas.de1(i) = eigval(1,:);
            eas.de2(i) = eigval(2,:);
            eas.de3(i) = eigval(3,:);

            eas.tort(i) = eas.de1(i)./(.5*(eas.de2(i) + eas.de3(i)));
            eas.de_perp(i) = (eas.de2(i) + eas.de3(i))/2;

            % eigenvalue decomposition of Da
            DTi = dt_i([1:3 2 4 5 3 5 6]);
            DTi = reshape(DTi, [3 3]);
            
            [~, eigval] = eigs(DTi); 
            eigval = sort(diag(eigval), 'descend');

            ias.da1(i) = eigval(1,:);
            ias.da2(i) = eigval(2,:);
            ias.da3(i) = eigval(3,:);

            ias.da_perp(i) = (ias.da2(i) + ias.da3(i))/2;
            ias.Da(i) = sum(eigval(:)); 
       % catch
%             ias.da1(i) = 0;
%             ias.da2(i) = 0;
%             ias.da3(i) = 0;
%             ias.da_perp(i) = 0;
%             ias.Da(i) = 0;
%             
%             eas.de1(i) = 0;
%             eas.de2(i) = 0;
%             eas.de3(i) = 0;
%             eas.de_perp(i) = 0;
%             eas.tort(i) = 0;

        %end
end
%keyboard



            
            
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
    [W_ind, W_cnt] = createTensorOrder(4);

    adc = ADC(dt(1:6, :), dir);
    md = sum(dt([1 4 6],:),1)/3;

    ndir  = size(dir, 1);
    T =  W_cnt(ones(ndir, 1), :).*dir(:,W_ind(:, 1)).*dir(:,W_ind(:, 2)).*dir(:,W_ind(:, 3)).*dir(:,W_ind(:, 4));

    akc =  T*dt(7:21, :);
    akc = (akc .* repmat(md.^2, [size(adc, 1), 1]))./(adc.^2);
end

function [adc] = ADC(dt, dir)
    [D_ind, D_cnt] = createTensorOrder(2);
    ndir  = size(dir, 1);
    T =  D_cnt(ones(ndir, 1), :).*dir(:,D_ind(:, 1)).*dir(:,D_ind(:, 2));
    adc = T * dt;
end

function [X, cnt] = createTensorOrder(order)
    
%     X = nchoosek(kron([1, 2, 3], ones(1, order)), order);
%     X = unique(X, 'rows');
%     for i = 1:size(X, 1)
%         cnt(i) = factorial(order) / factorial(nnz(X(i, :) ==1))/ factorial(nnz(X(i, :) ==2))/ factorial(nnz(X(i, :) ==3));
%     end

    if order == 2
        X = [1 1; ...
             1 2; ...
             1 3; ...
             2 2; ...
             2 3; ...
             3 3];
        cnt = [1 2 2 1 2 1]; 
    end
    
    if order == 4
        X = [1 1 1 1; ...
             1 1 1 2; ...
             1 1 1 3; ...
             1 1 2 2; ...
             1 1 2 3; ...
             1 1 3 3; ...
             1 2 2 2; ...
             1 2 2 3; ...
             1 2 3 3; ... 
             1 3 3 3; ...
             2 2 2 2; ...
             2 2 2 3; ...
             2 2 3 3; ...
             2 3 3 3; ...
             3 3 3 3];
        cnt = [1 4 4 6 12 6 4 12 12 4 1 4 6 4 1]; 
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