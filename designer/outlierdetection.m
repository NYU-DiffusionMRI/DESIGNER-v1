function akc_out = outlierdetection(dt,mask)
    if ~exist('mask','var') || isempty(mask)     
        mask = ~isnan(dt(:,:,:,1));
    end
    dt = vec(dt, mask);
    load('dirs10000.mat','dir');
    
    nvxls = size(dt, 2);
    akc_out = zeros([1 nvxls]);
    N = 10000;
    nblocks = 10;
    for i=1:nblocks
        akc = AKC(dt, dir((N/nblocks*(i-1))+1:N/nblocks*i, :)); %#ok<NODEF>
        akc_out = single(any(akc < -2 | akc > 10, 1));
    end
    %[akc, adc] = AKC(dt, dir);
    %akc_out = single(any(akc < -2 | akc > 10, 1));
    if exist('mask','var'), akc_out = vec(akc_out, mask); end
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

function [X, cnt] = createTensorOrder(order) 
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

function [adc] = ADC(dt, dir)
    [D_ind, D_cnt] = createTensorOrder(2);
    ndir  = size(dir, 1);
    T =  D_cnt(ones(ndir, 1), :).*dir(:,D_ind(:, 1)).*dir(:,D_ind(:, 2));
    adc = T * dt;
end

function [s, mask] = vec(S, mask)
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

