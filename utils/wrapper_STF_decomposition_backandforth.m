function dt_out = wrapper_STF_decomposition_backandforth(dt_in,direction,mask,inputOrder,CSphase)
% dt_out = wrapper_STF_decomposition_backandforth(dt_in,direction,mask,inputOrder,CSphase)
%
%
% S_ij is ordered:   [S_11;S_22;S_33;S_12;S_13;S_23]
% S_ijkl is ordered: [S_1111;S_2222;S_3333;S_1122;S_1133;S_1112;S_1113;S_1123;S_2233;S_2212;S_2213;S_2223;S_3312;S_3313;S_3323];
%
% Dlm is ordered:    [S_{0,0}, S_{2,-2}, S_{2,-1}, S_{2,0}, S_{2,1}, S_{2,2}]
% Slm is ordered:    [S_{0,0}, S_{2,-2}, S_{2,-1}, S_{2,0}, S_{2,1}, S_{2,2}, S_{4,-4}, S_{4,-1}, S_{4,-3}, S_{4,-1}, S_{4,0}, S_{4,1}, S_{4,2}, S_{4,3}, S_{4,4}]
%
% dt_in can be a rank 2 or rank 4 tensor (6 or 15 independent elements)
% direction can be 'STF2cart' or 'cart2STF'
%
% If CSphase flag is empty or 0, this function does not use the Condon-Shortley
% phase in the definition of spherical harmonics (if 1 then yes)
%
%
% For more details see Rotational Invariants of the Cumulant Expansion
% (RICE), ISMRM 2022, by Santiago Coelho et. al
%
% By: Santiago Coelho (03/09/2021) Santiago.Coelho@nyulangone.org

if ~exist('CSphase','var') || isempty(CSphase) || ~CSphase
    CSphase=0; % 0 means we DO NOT use it
    fprintf('Not using the Condon-Shortley phase on the spherical harmonics definition!\n')
else
    CSphase=1; % 1 means we use it
    fprintf('Using the Condon-Shortley phase on the spherical harmonics definition!\n')
end

is_4D_array=length(size(dt_in))==4;

if is_4D_array
    [x, y, z, nelem] = size(dt_in);
    if ~exist('mask','var') || isempty(mask)
        mask = true(x, y, z);
    end
    dt_in=vectorize(dt_in,mask);
else
    nelem=size(dt_in,1);
end

if exist('inputOrder','var') && ~isempty(inputOrder) && strcmp(inputOrder,'JV') && nelem==6
    dt_in=dt_in([1 4 6 2 3 5],:);
elseif exist('inputOrder','var') && ~isempty(inputOrder) && strcmp(inputOrder,'JV') && nelem==15
    dt_in=dt_in([1 11 15 4 6 2 3 5 13 7 8 12 9 10 14],:);
end

if strcmp(direction,'STF2cart')
    if nelem==6
        dt_out = STF2cart_rank2_tensor(dt_in,CSphase);
    elseif nelem==15
        dt_out = STF2cart_rank4_tensor(dt_in,CSphase);
    else
        error('Rank 2 or rank 4 tensors are only allowed in this function')
    end
elseif strcmp(direction,'cart2STF')
    if nelem==6
        dt_out = IrreducibleDecomposition_rank2_tensor(dt_in,CSphase);
    elseif nelem==15
        dt_out = IrreducibleDecomposition_rank4_tensor(dt_in,CSphase);
    else
        error('Rank 2 or rank 4 tensors are only allowed in this function')
    end
else
    error('direction must be ''STF2cart'' or ''cart2STF''')
end

if is_4D_array
    dt_out=vectorize(dt_out,mask);
end

end

function [Slm] = IrreducibleDecomposition_rank2_tensor(S_ij,CSphase)
% S: 2D array containing [S_11;S_22;S_33;S_12;S_13;S_23];   
%
C0=sqrt(1/(4*pi));
C2=sqrt(5/(4*pi));
Slm=0*S_ij;
Y_ij_00 =C0*eye(3);
Y_ij_2m2=[0,0.546274215296040,0;0.546274215296040,0,0;0,0,0];
Y_ij_2m1=[0,0,0;0,0,-0.546274215296040;0,-0.546274215296040,0];
Y_ij_20 =[-0.315391565252520,0,0;0,-0.315391565252520,0;0,0,0.630783130505040];
Y_ij_21 =[0,0,-0.546274215296040;0,0,0;-0.546274215296040,0,0];
Y_ij_22 =[0.546274215296040,0,0;0,-0.546274215296040,0;0,0,0];
Y_ij_00_vector=Y_ij_00([1 5 9 4 7 8]).*[1 1 1 2 2 2];
Y_ij_2m2_vector=Y_ij_2m2([1 5 9 4 7 8]).*[1 1 1 2 2 2];
Y_ij_2m1_vector=Y_ij_2m1([1 5 9 4 7 8]).*[1 1 1 2 2 2];
Y_ij_20_vector=Y_ij_20([1 5 9 4 7 8]).*[1 1 1 2 2 2];
Y_ij_21_vector=Y_ij_21([1 5 9 4 7 8]).*[1 1 1 2 2 2];
Y_ij_22_vector=Y_ij_22([1 5 9 4 7 8]).*[1 1 1 2 2 2];
if CSphase
    Y_ij_2m1_vector=-Y_ij_2m1_vector;
    Y_ij_21_vector=-Y_ij_21_vector;
end
Slm(1,:) =1/3*(1/C0^2)*Y_ij_00_vector*S_ij;
% Slm(1,:)=1/3*(1/C0)*sum(S_ij(1:3,:),1);
Slm(2,:)=2/3*(1/C2^2)*Y_ij_2m2_vector*S_ij;
Slm(3,:)=2/3*(1/C2^2)*Y_ij_2m1_vector*S_ij;
Slm(4,:)=2/3*(1/C2^2)*Y_ij_20_vector*S_ij;
Slm(5,:)=2/3*(1/C2^2)*Y_ij_21_vector*S_ij;
Slm(6,:)=2/3*(1/C2^2)*Y_ij_22_vector*S_ij;
end


function [Slm] = IrreducibleDecomposition_rank4_tensor(S_ijkl,CSphase)
% Sijkl: 2D array containing [S_1111;S_2222;S_3333;S_1122;S_1133;...
%                             S_1112;S_1113;S_1123;S_2233;S_2212;...
%                             S_2213;S_2223;S_3312;S_3313;S_3323];   
%
%
Slm=0*S_ijkl;
C0=sqrt(1/(4*pi));
C2=sqrt(5/(4*pi));
C4=sqrt(9/(4*pi));

Y_ijkl_00_vector = [0.282094791773878,0.282094791773878,0.282094791773878,0.564189583547756,0.564189583547756,0,0,0,0.564189583547756,0,0,0,0,0,0];
Y_ijkl_2m2_vector= [0,0,0,0,0,1.09254843059208,0,0,0,1.09254843059208,0,0,1.09254843059208,0,0];
Y_ijkl_2m1_vector= [0,0,0,0,0,0,0,-1.09254843059208,0,0,0,-1.09254843059208,0,0,-1.09254843059208];
Y_ijkl_20_vector = [-0.315391565252520,-0.315391565252520,0.630783130505040,-0.630783130505040,0.315391565252520,0,0,0,0.315391565252520,0,0,0,0,0,0];
Y_ijkl_21_vector = [0,0,0,0,0,0,-1.09254843059208,0,0,0,-1.09254843059208,0,0,-1.09254843059208,0];
Y_ijkl_22_vector = [0.546274215296039,-0.546274215296039,0,0,0.546274215296040,0,0,0,-0.546274215296040,0,0,0,0,0,0];
Y_ijkl_4m4_vector= [0,0,0,0,0,2.50334294179671,0,0,0,-2.50334294179671,0,0,0,0,0];
Y_ijkl_4m3_vector= [0,0,0,0,0,0,0,-5.31039230933979,0,0,0,1.77013076977993,0,0,0];
Y_ijkl_4m2_vector= [0,0,0,0,0,-0.946174695757560,0,0,0,-0.946174695757560,0,0,5.67704817454536,0,0];
Y_ijkl_4m1_vector= [0,0,0,0,0,0,0,2.00713963067187,0,0,0,2.00713963067187,0,0,-2.67618617422916];
Y_ijkl_40_vector = [0.317356640745613,0.317356640745613,0.846284375321635,0.634713281491226,-2.53885312596490,0,0,0,-2.53885312596490,0,0,0,0,0,0];
Y_ijkl_41_vector = [0,0,0,0,0,0,2.00713963067187,0,0,0,2.00713963067187,0,0,-2.67618617422916,0];
Y_ijkl_42_vector = [-0.473087347878780,0.473087347878780,0,0,2.83852408727268,0,0,0,-2.83852408727268,0,0,0,0,0,0];
Y_ijkl_43_vector = [0,0,0,0,0,0,-1.77013076977993,0,0,0,5.31039230933979,0,0,0,0];
Y_ijkl_44_vector = [0.625835735449177,0.625835735449177,0,-3.75501441269506,0,0,0,0,0,0,0,0,0,0,0];
if CSphase
    Y_ijkl_2m1_vector=-Y_ijkl_2m1_vector;
    Y_ijkl_21_vector=-Y_ijkl_21_vector;
    Y_ijkl_4m3_vector=-Y_ijkl_4m3_vector;
    Y_ijkl_4m1_vector=-Y_ijkl_4m1_vector;
    Y_ijkl_41_vector=-Y_ijkl_41_vector;
    Y_ijkl_43_vector=-Y_ijkl_43_vector;
end
Slm(1,:) =1/5*(1/C0^2)*Y_ijkl_00_vector*S_ijkl;
Slm(2,:) =4/7*(1/C2^2)*Y_ijkl_2m2_vector*S_ijkl;
Slm(3,:) =4/7*(1/C2^2)*Y_ijkl_2m1_vector*S_ijkl;
Slm(4,:) =4/7*(1/C2^2)*Y_ijkl_20_vector*S_ijkl;
Slm(5,:) =4/7*(1/C2^2)*Y_ijkl_21_vector*S_ijkl;
Slm(6,:) =4/7*(1/C2^2)*Y_ijkl_22_vector*S_ijkl;
Slm(7,:) =8/35*(1/C4^2)*Y_ijkl_4m4_vector*S_ijkl;
Slm(8,:) =8/35*(1/C4^2)*Y_ijkl_4m3_vector*S_ijkl;
Slm(9,:) =8/35*(1/C4^2)*Y_ijkl_4m2_vector*S_ijkl;
Slm(10,:)=8/35*(1/C4^2)*Y_ijkl_4m1_vector*S_ijkl;
Slm(11,:)=8/35*(1/C4^2)*Y_ijkl_40_vector*S_ijkl;
Slm(12,:)=8/35*(1/C4^2)*Y_ijkl_41_vector*S_ijkl;
Slm(13,:)=8/35*(1/C4^2)*Y_ijkl_42_vector*S_ijkl;
Slm(14,:)=8/35*(1/C4^2)*Y_ijkl_43_vector*S_ijkl;
Slm(15,:)=8/35*(1/C4^2)*Y_ijkl_44_vector*S_ijkl;
end

function [Sij] = STF2cart_rank2_tensor(S2m,CSphase)
% Sij: 2D array containing [S_11;S_22;S_33;S_12;S_13;S_23]
C0=sqrt(1/(4*pi));
Y_ij_00=C0*eye(3);
Y_ij_2m2=[0,0.546274215296040,0;0.546274215296040,0,0;0,0,0];
Y_ij_2m1=[0,0,0;0,0,-0.546274215296040;0,-0.546274215296040,0];
Y_ij_20 =[-0.315391565252520,0,0;0,-0.315391565252520,0;0,0,0.630783130505040];
Y_ij_21 =[0,0,-0.546274215296040;0,0,0;-0.546274215296040,0,0];
Y_ij_22 =[0.546274215296040,0,0;0,-0.546274215296040,0;0,0,0];
if CSphase
    Y_ij_2m1=-Y_ij_2m1;
    Y_ij_21=-Y_ij_21;
end
Y_ij_all=cat(3,Y_ij_00,Y_ij_2m2,Y_ij_2m1,Y_ij_20,Y_ij_21,Y_ij_22);
idx=[1 1;2 2;3 3;1 2;1 3;2 3];
Sij=0*S2m;
for ii=1:6
    for jj=1:6
        Sij(ii,:)=Sij(ii,:)+S2m(jj,:)*Y_ij_all(idx(ii,1),idx(ii,2),jj);
    end
end
end

% function [Sij] = STF2cart_rank2_tensor(S2m)
% %
% C0=sqrt(1/(4*pi));
% Sij=0*S2m;
% Y_ij_00=C0*eye(3);
% % Y_ij_00=1/3*eye(3);
% Y_ij_2m2=[0,0.546274215296040,0;0.546274215296040,0,0;0,0,0];
% Y_ij_2m1=[0,0,0;0,0,-0.546274215296040;0,-0.546274215296040,0];
% Y_ij_20 =[-0.315391565252520,0,0;0,-0.315391565252520,0;0,0,0.630783130505040];
% Y_ij_21 =[0,0,-0.546274215296040;0,0,0;-0.546274215296040,0,0];
% Y_ij_22 =[0.546274215296040,0,0;0,-0.546274215296040,0;0,0,0];
% Y_ij_all=cat(3,Y_ij_00,Y_ij_2m2,Y_ij_2m1,Y_ij_20,Y_ij_21,Y_ij_22);
% idx=[1 1;2 2;3 3;1 2;1 3;2 3];
% for ii=1:6
%     for jj=1:6
%         Sij(ii,:)=Sij(ii,:)+S2m(jj,:).*Y_ij_all(idx(ii,1),idx(ii,2),jj);
%     end
% end
% % Sij(1,:)=S2m(1,:).*delta(1,1)+S2m(2,:).*Y_ij_2m2(1,1)+S2m(3,:).*Y_ij_2m1(1,1)+S2m(4,:).*Y_ij_20(1,1)+S2m(5,:).*Y_ij_21(1,1)+S2m(6,:).*Y_ij_22(1,1);
% % Sij(2,:)=S2m(1,:).*delta(2,2)+S2m(2,:).*Y_ij_2m2(2,2)+S2m(3,:).*Y_ij_2m1(2,2)+S2m(4,:).*Y_ij_20(2,2)+S2m(5,:).*Y_ij_21(2,2)+S2m(6,:).*Y_ij_22(2,2);
% % Sij(3,:)=S2m(1,:).*delta(3,3)+S2m(2,:).*Y_ij_2m2(3,3)+S2m(3,:).*Y_ij_2m1(3,3)+S2m(4,:).*Y_ij_20(3,3)+S2m(5,:).*Y_ij_21(3,3)+S2m(6,:).*Y_ij_22(3,3);
% % Sij(4,:)=S2m(1,:).*delta(1,2)+S2m(2,:).*Y_ij_2m2(1,2)+S2m(3,:).*Y_ij_2m1(1,2)+S2m(4,:).*Y_ij_20(1,2)+S2m(5,:).*Y_ij_21(1,2)+S2m(6,:).*Y_ij_22(1,2);
% % Sij(5,:)=S2m(1,:).*delta(1,3)+S2m(2,:).*Y_ij_2m2(1,3)+S2m(3,:).*Y_ij_2m1(1,3)+S2m(4,:).*Y_ij_20(1,3)+S2m(5,:).*Y_ij_21(1,3)+S2m(6,:).*Y_ij_22(1,3);
% % Sij(6,:)=S2m(1,:).*delta(2,3)+S2m(2,:).*Y_ij_2m2(2,3)+S2m(3,:).*Y_ij_2m1(2,3)+S2m(4,:).*Y_ij_20(2,3)+S2m(5,:).*Y_ij_21(2,3)+S2m(6,:).*Y_ij_22(2,3);
% end


function [Sijkl] = STF2cart_rank4_tensor(S4m,CSphase)
% Sijkl: 2D array containing [S_1111;S_2222;S_3333;S_1122;S_1133;S_1112;S_1113;S_1123;S_2233;S_2212;S_2213;S_2223;S_3312;S_3313;S_3323];   
load('/Users/coelhs01/Documents/SantiagoCoelho/NYU_Postdoc_MyScience/MATLAB/General_Scripts/Y_ijkl_all.mat')
idx=[ 1,1,1,1; 2,2,2,2; 3,3,3,3; 1,1,2,2; 1,1,3,3; 1,1,1,2; 1,1,1,3; 1,1,2,3; 2,2,3,3; 2,2,1,2; 2,2,1,3; 2,2,2,3; 3,3,1,2; 3,3,1,3; 3,3,2,3];
Sijkl=0*S4m;
signs=ones(15,1);
if CSphase
    signs([3 5 8 10 12 14])=-1;
end
for ii=1:15
    for jj=1:15
        Sijkl(ii,:)=Sijkl(ii,:)+S4m(jj,:)*Y_ijkl_all(idx(ii,1),idx(ii,2),idx(ii,3),idx(ii,4),jj)*signs(jj);
    end
end
% clc,clear,close all
% y=getY(4);
% Y_ijkl_all(:,:,:,:,1)=y{2,1,1};
% % y{2,1,1}(:)'*y{2,1,1}(:)*1/5*(1/C0^2)
% 1/5*y{2,1,1}(1,1,1,1)+y{2,1,1}(2,2,2,2)+y{2,1,1}(3,3,3,3)+2*(y{2,1,1}(1,1,2,2)+y{2,1,1}(1,1,3,3)+y{2,1,1}(2,2,3,3))/sqrt(4*pi)
% 
% Y_ijkl_all(:,:,:,:,2)=y{2,2,1};
% Y_ijkl_all(:,:,:,:,3)=y{2,2,2};
% Y_ijkl_all(:,:,:,:,4) =y{2,2,3};
% Y_ijkl_all(:,:,:,:,5) =y{2,2,4};
% Y_ijkl_all(:,:,:,:,6) =y{2,2,5};
% Y_ijkl_all(:,:,:,:,7)=y{2,3,1};
% Y_ijkl_all(:,:,:,:,8)=y{2,3,2};
% Y_ijkl_all(:,:,:,:,9)=y{2,3,3};
% Y_ijkl_all(:,:,:,:,10)=y{2,3,4};
% Y_ijkl_all(:,:,:,:,11) =y{2,3,5};
% Y_ijkl_all(:,:,:,:,12) =y{2,3,6};
% Y_ijkl_all(:,:,:,:,13) =y{2,3,7};
% Y_ijkl_all(:,:,:,:,14) =y{2,3,8};
% Y_ijkl_all(:,:,:,:,15) =y{2,3,9};
% save('/Users/coelhs01/Documents/SantiagoCoelho/NYU_Postdoc_MyScience/MATLAB/General_Scripts/Y_ijkl_all.mat','Y_ijkl_all')
end






function [s, mask] = vectorize(S, mask)
if nargin == 1
   mask = ~isnan(S(:,:,:,1));
end
if ismatrix(S)
n = size(S, 1);
[x, y, z] = size(mask);
s = NaN([x, y, z, n], 'like', S);
for ii = 1:n
    tmp = NaN(x, y, z, 'like', S);
    tmp(mask(:)) = S(ii, :);
    s(:,:,:,ii) = tmp;
end
else
for ii = 1:size(S, 4)
   Si = S(:,:,:,ii);
   s(ii, :) = Si(mask(:));
end
end
end
