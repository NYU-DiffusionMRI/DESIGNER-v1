function [b0, dt, md, rd, ad, fa] = dti_fit(dwi, grad, mask)

    %% check grad order
    bval = grad(:, 4);
    order = floor(log(abs(max(bval)+1))./log(10));
    if order >= 2
        grad(:, 4) = grad(:, 4)/1000;
        bval = grad(:, 4);
    end
   
    %% parameter checks 
    dwi = double(dwi);
    dwi(dwi<=0)=eps;
    [x, y, z, ndwis] = size(dwi);
    if ~exist('grad','var') || size(grad,1) ~= ndwis || size(grad,2) ~= 4
        error('');
    end
    grad = double(grad);
    normgrad = sqrt(sum(grad(:, 1:3).^2, 2)); normgrad(normgrad == 0) = 1;
    grad(:, 1:3) = grad(:, 1:3)./repmat(normgrad, [1 3]);
    grad(isnan(grad)) = 0;

    bval = grad(:, 4);

    if ~exist('mask','var') || isempty(mask)
        mask = true(x, y, z);
    end
    
    dwi = vectorize(dwi, mask);
   
    %% tensor fit
    [D_ind, D_cnt] = createTensorOrder(2);
    
    bS = ones(ndwis, 1);
    bD = D_cnt(ones(ndwis, 1), :).*grad(:,D_ind(:, 1)).*grad(:,D_ind(:, 2));
    
    b = [bS, -bval(:, ones(1, 6)).*bD, (bval(:, ones(1, 15)).^2/6).*bW];

    % unconstrained LLS fit
    dt = b\log(dwi);
    w = exp(b*dt);

    nvoxels = size(dwi,2);

    % WLLS fit initialized with LLS   
    parfor i = 1:nvoxels
        wi = diag(w(:,i)); 
        logdwii = log(dwi(:,i));
        dt(:,i) = (wi*b)\(wi*logdwii);
    end

    b0 = exp(dt(1,:));
    dt = dt(2:7, :);
    b0 = vectorize(b0, mask);
    dt = vectorize(dt, mask);
    
    DT = reshape([dt(1,:); dt(2,:); dt(3,:);...
        dt(2,:); dt(4,:); dt(5,:);...
        dt(3,:); dt(5,:); dt(6,:);],[3,3,size(dt,2)]);

    for i=1:size(dt,2)
        [v,l] = eig(DT(:,:,i));
        l = diag(l); 
        [l, idx] = sort(l,'descend'); 
        v = v(:, idx);
        
        v1(i,:)=v(:,1);
        v2(i,:)=v(:,2);
        v3(i,:)=v(:,3);
        l1(i)=abs(l(1));
        l2(i)=abs(l(2));
        l3(i)=abs(l(3));
    end
    md = vectorize((l1+l2+l3)./3, mask);
    rd = vectorize((l2+l3)./2, mask);
    ad = vectorize(l1, mask);
    fa = vectorize(sqrt(1/2).*sqrt((l1-l2).^2+(l2-l3).^2+(l3-l1).^2)./sqrt(l1.^2+l2.^2+l3.^2), mask);
    

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

