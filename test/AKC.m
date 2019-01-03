function [akc, adc] = AKC(dt, dir)

[W_ind, W_cnt] = createTensorOrder(4);

adc = ADC(dt(1:6, :), dir);
md = sum(dt([1 4 6],:),1)/3;

ndir  = size(dir, 1);
T =  W_cnt(ones(ndir, 1), :).*dir(:,W_ind(:, 1)).*dir(:,W_ind(:, 2)).*dir(:,W_ind(:, 3)).*dir(:,W_ind(:, 4));
 
akc =  T*dt(7:21, :);
akc = (akc .* repmat(md.^2, [size(adc, 1), 1]))./(adc.^2);
end