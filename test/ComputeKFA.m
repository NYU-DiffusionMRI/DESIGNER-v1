function [kfa,mkt] = ComputeKFA(dt,Kmax_final,Kmin_final)
%
% computes the kfa given a 15-vector of kurtosis tensor values
%
% Author: Mark Van Horn, Emilie McKinnon, and Siddhartha Dhiman 
% Last modified: 01/03/19

offset = 6; %   Number of DT elements prior to KT elements in dt

W1111 = dt(1+offset);
W1112 = dt(2+offset);
W1113 = dt(3+offset);
W1122 = dt(4+offset);
W1123 = dt(5+offset);
W1133 = dt(6+offset);
W1222 = dt(7+offset);
W1223 = dt(8+offset);
W1233 = dt(9+offset);
W1333 = dt(10+offset);
W2222 = dt(11+offset);
W2223 = dt(12+offset);
W2233 = dt(13+offset);
W2333 = dt(14+offset);
W3333 = dt(15+offset);

W_F = sqrt(W1111 ^ 2 + W2222 ^ 2 + W3333 ^ 2 + 6 * W1122 ^ 2 + 6 * W1133 ^ 2 + 6 * W2233 ^ 2 + ...
    4 * W1112 ^ 2 + 4 * W1113 ^ 2 + 4 * W1222 ^ 2 + 4 * W2223 ^ 2 + 4 * W1333 ^ 2 + ...
    4 * W2333 ^ 2 + 12 * W1123 ^ 2 + 12 * W1223 ^ 2 + 12 * W1233 ^ 2);

Wbar = 1/5 * (W1111 + W2222+ W3333 + 2 * (W1122 + W1133 + W2233));

if W_F < 1e-3,
    kfa = 0;
else
    W_diff_F = sqrt((W1111 - Wbar) ^ 2 + (W2222 - Wbar) ^ 2 + (W3333 - Wbar) ^ 2 + ...
        6 * (W1122 - Wbar / 3) ^ 2 + 6 * (W1133 - Wbar / 3) ^ 2 + 6 * (W2233 - Wbar / 3) ^ 2 + ...
        4 * W1112 ^ 2 + 4 * W1113 ^ 2 + 4 * W1222 ^ 2 + 4 * W2223 ^ 2 + ...
        4 * W1333 ^ 2 + 4 * W2333 ^ 2 + 12 * W1123 ^ 2 + 12 * W1223 ^ 2 + 12 * W1233 ^ 2);
    
    kfa = W_diff_F / W_F;
end

mkt=Wbar;
mkt(mkt > Kmax_final) = Kmax_final;
mkt(mkt < Kmin_final) = Kmin_final;
end
