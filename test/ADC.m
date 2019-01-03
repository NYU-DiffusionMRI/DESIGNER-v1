function [adc] = ADC(dt, dir)
[D_ind, D_cnt] = createTensorOrder(2);
ndir  = size(dir, 1);
T =  D_cnt(ones(ndir, 1), :).*dir(:,D_ind(:, 1)).*dir(:,D_ind(:, 2));
adc = T * dt;
end