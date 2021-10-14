function [b0,Dlm,Wlm,Dl,Wl,DTI_scalars,DKI_scalars,ExtraScalars] = Fit_DKI_LTE(DWI,b,dirs, mask, computeExtra)
% [b0,Dlm,Wlm,Dl,Wl,DTI_scalars,DKI_scalars,ExtraScalars] = Fit_DKI_LTE(DWI,b,dirs, mask, computeExtra)
%
% This function fits DKI on a 4D diffusion dataset (x y z ndwis)
%
% Parameter estimation:
% WLLS estimation on the log(dwi) (initialized with LLS on the log)
%
%
% The following are outputs:
% b0: b0 estimate
% Dlm: 4D array containing D_{0,0}, D_{2,-2}, D_{2,-1}, D_{2,0}, D_{2,1}, D_{2,2}
% Wlm: 4D array containing W_{0,0}, W_{2,-2}, W_{2,-1}, W_{2,0}, W_{2,1}, W_{2,2}, W_{4,-4}, W_{4,-1}, W_{4,-3}, W_{4,-1}, W_{4,0}, W_{4,1}, W_{4,2}, W_{4,3}, W_{4,4}
% Dl: 4D array containing D_0 and D_2
% Wl: 4D array containing W_0, W_2 and W_4
% DTI_scalars.fa: Fractional anisotropy of the diffusion tensor
% DTI_scalars.md: Mean diffusivity of the diffusion tensor
% DTI_scalars.ad_axsym: Axial diffusivity of the diffusion tensor (Assuming
% axial symmetry so there is no projection to the fibre basis)
% DTI_scalars.rd_axsym: Radial diffusivity of the diffusion tensor (Assuming
% axial symmetry so there is no projection to the fibre basis)
% DKI_scalars.mw: Spherical mean of the W tensor (its trace)
% DKI_scalars.aw_axsym: Axial projection of the kurtosis tensor (Assuming
% axial symmetry so there is no projection to the fibre basis) 
% DKI_scalars.rw_axsym: Radial projection of the kurtosis tensor (Assuming
% axial symmetry so there is no projection to the fibre basis) 
%
% Whenever kurtosis tensor is mentioned it is talking about W_{ijkl} (K is
% not a 4th order tensor)
%
% If computeExtra==1, the following parameters are also computed
%     ExtraScalars.rd radial diffusivity
%     ExtraScalars.ad axial diffusivity
%     ExtraScalars.mk mean(k(n)), mean directional kurtosis
%     ExtraScalars.aw Projection of W parallel to the fibre axis
%     ExtraScalars.rw Projection of W perpendicular to the fibre axis
%     ExtraScalars.e1  Main fibre axis
%
%
% By: Santiago Coelho (10/06/2021) Santiago.Coelho@nyulangone.org
% =========================================================================
[x, y, z, ndwis] = size(DWI);
if ~exist('mask','var') || isempty(mask)
    mask = true(x, y, z);
end
DWI = vectorize(DWI, mask);
if size(dirs,2)~=3
    dirs=dirs';
end
Nvoxels=size(DWI,2);
b=b(:);
Nmeas=length(b);
if size(dirs,1)~=Nmeas
    error('Number of directions must much number of b-values')
end
% Evaluate SH on input directions
Y_LM_matrix = get_even_SH(dirs,4);
% Y00=Y_LM_matrix(:,1);
Y2m=Y_LM_matrix(:,2:6);
Y4m=Y_LM_matrix(:,7:15);
X=[ones(size(b)) -b -b.*Y2m 1/6*b.^2 1/6*b.^2.*Y2m 1/6*b.^2.*Y4m];
% unconstrained LLS fit
dv = X\log(DWI);
w = exp(X*dv);
% WLLS fit initialized with LLS
parfor ii = 1:Nvoxels
    wi = diag(w(:,ii));
    logdwii = log(DWI(:,ii));
    dv(:,ii) = (wi*X)\(wi*logdwii);
end

% Recover b0
b0 = exp(dv(1,:));
b0 = vectorize(b0, mask);
% Recover Diffusion and Kurtosis tensor elements
dt=dv(2:22,:);
MDSq = dt(1,:).^2;
Dlm=vectorize(dt(1:6,:), mask);
Wlm=vectorize(dt(7:21,:)./MDSq, mask);
% Compute Diffusion and Kurtosis rotational invariants
D0 = Dlm(:,:,:,1);
D2=sqrt(sum(Dlm(:,:,:,2:6).^2,4))/sqrt(5*4*pi);
W0 = Wlm(:,:,:,1);
W2=sqrt(sum(Wlm(:,:,:,2:6).^2,4))/sqrt(5*4*pi);
W4=sqrt(sum(Wlm(:,:,:,7:15).^2,4))/sqrt(9*4*pi);
Dl=cat(4,D0,D2);
Wl=cat(4,W0,W2,W4);

% Computing scalars from the diffusion tensor
D00_sq=(3*D0/sqrt(4*pi)).^2;
D2m_norm_sq=sum(Dlm(:,:,:,2:6).^2,4)*sqrt(2/(3*5));
%fa=sqrt(3*D2m_norm_sq./(5*D00_sq+2*D2m_norm_sq));
fa=sqrt(75*D2.^2./(4*D0.^2+50*D2.^2));

DTI_scalars.fa=fa;
DTI_scalars.md=D0;
DTI_scalars.ad_axsym=D0+5*D2;
DTI_scalars.rd_axsym=D0-5/2*D2;
DKI_scalars.mw=W0;
DKI_scalars.aw_axsym=(W0+5*W2+9*W4).*D0.^2./(D0+5*D2).^2;
DKI_scalars.rw_axsym=(W0-5/2*W2+27/8*W4).*D0.^2./(D0-5/2*D2).^2;

if ~exist('computeExtra','var') || isempty(computeExtra) || ~computeExtra
    ExtraScalars=[];
else
    % Computing AW and RW with their exact definitions
    Y2m2=[0,0.546274215296040,0;0.546274215296040,0,0;0,0,0];
    Y2m1=[0,0,0;0,0,-0.546274215296040;0,-0.546274215296040,0];
    Y20= [-0.315391565252520,0,0;0,-0.315391565252520,0;0,0,0.630783130505040];
    Y21= [0,0,-0.546274215296040;0,0,0;-0.546274215296040,0,0];
    Y22= [0.546274215296040,0,0;0,-0.546274215296040,0;0,0,0];
    wlm=dt(7:21,:)./MDSq;
    AW=0*dt(1, :);
    e1=0*dt(1:3, :);
    l1=0*dt(1, :);
    l2=0*dt(1, :);
    l3=0*dt(1, :);
    parfor ii = 1:Nvoxels
        DT = dt(1,ii)*eye(3)+(dt(2, ii)*Y2m2+dt(3, ii)*Y2m1+dt(4, ii)*Y20+dt(5, ii)*Y21+dt(6, ii)*Y22);
        DT = reshape(DT, [3 3]);
        [eigvec, eigval] = eigs(DT);
        e1(:, ii) = eigvec(:, 1);
        eigval = diag(eigval);
        l1(ii) = eigval(1,:);
        l2(ii) = eigval(2,:);
        l3(ii) = eigval(3,:);
        Ylm_e1 = get_even_SH(eigvec(:, 1),4);
        AW(ii)=wlm(1,ii)+Ylm_e1(1,2:15)*wlm(2:15,ii);
        dirs_radial = radialsampling(eigvec(:, 1), 100);
        Ylm_radial = get_even_SH(dirs_radial,4);
        RW(ii)=mean(wlm(1,ii)+Ylm_radial(:,2:15)*wlm(2:15,ii));
    end
    RD = vectorize((l2+l3)/2, mask);
    AD = vectorize(l1, mask);
    AW = vectorize(AW, mask);
    RW = vectorize(RW, mask);
    AW=AW.*D0.^2./AD.^2;
    RW=RW.*D0.^2./RD.^2;

    dirs2 = get256dirs();
    Y_LM_matrix2 = get_even_SH(dirs2,4);
    Y2m_2=Y_LM_matrix2(:,2:6);
    Y4m_2=Y_LM_matrix2(:,7:15);
    adc = [ones(256,1) Y2m_2]*dt(1:6,:);
    akc = [ones(256,1) Y2m_2 Y4m_2]*(dt(7:21,:))./(adc.^2);
    
    ExtraScalars.rd=RD;
    ExtraScalars.ad=AD;
    ExtraScalars.mk=vectorize(mean(akc), mask);
    ExtraScalars.aw=AW;
    ExtraScalars.rw=RW;
    ExtraScalars.e1=vectorize(e1, mask);
end



end

function Ylm_n = get_even_SH(dirs,Lmax)
% Ylm_n = get_even_SH(dirs,Lmax)
%
%
if size(dirs,2)~=3
    dirs=dirs';
end
Nmeas=size(dirs,1);
[PHI,THETA,~]=cart2sph(dirs(:,1),dirs(:,2),dirs(:,3)); THETA=pi/2-THETA;
l=0:2:Lmax;
l_all=[];
m_all=[];
for ii=1:length(l)
    l_all=[l_all, l(ii)*ones(1,2*l(ii)+1)];
    m_all=[m_all -l(ii):l(ii)];
end
K_lm=sqrt((2*l_all+1)./(4*pi) .* factorial(l_all-abs(m_all))./factorial(l_all+abs(m_all)));
extra_factor=ones(size(K_lm));
extra_factor(m_all~=0)=sqrt(2);
P_l_in_cos_theta=zeros(length(l_all),Nmeas);
phi_term=zeros(length(l_all),Nmeas);
id_which_pl=zeros(1,length(l_all));
for ii=1:length(l_all)
    all_Pls=legendre(l_all(ii),cos(THETA));
    P_l_in_cos_theta(ii,:)=all_Pls(abs(m_all(ii))+1,:);
    id_which_pl(ii)=abs(m_all(ii))+1;
    if m_all(ii)>0
        phi_term(ii,:)=cos(m_all(ii)*PHI);
    elseif m_all(ii)==0
        phi_term(ii,:)=1;
    elseif m_all(ii)<0
        phi_term(ii,:)=sin(-m_all(ii)*PHI);
    end
end
Y_lm=repmat(extra_factor',1,Nmeas).*repmat(K_lm',1,Nmeas).*phi_term.*P_l_in_cos_theta;
Ylm_n=Y_lm';
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

function dirs = get256dirs()
% get 256 isotropically distributed directions
dirs =  [0         0    1.0000;
    0.5924         0    0.8056;
    -0.7191   -0.1575   -0.6768;
    -0.9151   -0.3479    0.2040;
    0.5535    0.2437    0.7964;
    -0.0844    0.9609   -0.2636;
    0.9512   -0.3015    0.0651;
    -0.4225    0.8984    0.1202;
    0.5916   -0.6396    0.4909;
    0.3172    0.8818   -0.3489;
    -0.1988   -0.6687    0.7164;
    -0.2735    0.3047   -0.9123;
    0.9714   -0.1171    0.2066;
    -0.5215   -0.4013    0.7530;
    -0.3978   -0.9131   -0.0897;
    0.2680    0.8196    0.5063;
    -0.6824   -0.6532   -0.3281;
    0.4748   -0.7261   -0.4973;
    0.4504   -0.4036    0.7964;
    -0.5551   -0.8034   -0.2153;
    0.0455   -0.2169    0.9751;
    0.0483    0.5845    0.8099;
    -0.1909   -0.1544   -0.9694;
    0.8383    0.5084    0.1969;
    -0.2464    0.1148    0.9623;
    -0.7458    0.6318    0.2114;
    -0.0080   -0.9831   -0.1828;
    -0.2630    0.5386   -0.8005;
    -0.0507    0.6425   -0.7646;
    0.4476   -0.8877    0.1081;
    -0.5627    0.7710    0.2982;
    -0.3790    0.7774   -0.5020;
    -0.6217    0.4586   -0.6350;
    -0.1506    0.8688   -0.4718;
    -0.4579    0.2131    0.8631;
    -0.8349   -0.2124    0.5077;
    0.7682   -0.1732   -0.6163;
    0.0997   -0.7168   -0.6901;
    0.0386   -0.2146   -0.9759;
    0.9312    0.1655   -0.3249;
    0.9151    0.3053    0.2634;
    0.8081    0.5289   -0.2593;
    -0.3632   -0.9225    0.1305;
    0.2709   -0.3327   -0.9033;
    -0.1942   -0.9790   -0.0623;
    0.6302   -0.7641    0.1377;
    -0.6948   -0.3137    0.6471;
    -0.6596   -0.6452    0.3854;
    -0.9454    0.2713    0.1805;
    -0.2586   -0.7957    0.5477;
    -0.3576    0.6511    0.6695;
    -0.8490   -0.5275    0.0328;
    0.3830    0.2499   -0.8893;
    0.8804   -0.2392   -0.4095;
    0.4321   -0.4475   -0.7829;
    -0.5821   -0.1656    0.7961;
    0.3963    0.6637    0.6344;
    -0.7222   -0.6855   -0.0929;
    0.2130   -0.9650   -0.1527;
    0.4737    0.7367   -0.4825;
    -0.9956    0.0891    0.0278;
    -0.5178    0.7899   -0.3287;
    -0.8906    0.1431   -0.4317;
    0.2431   -0.9670    0.0764;
    -0.6812   -0.3807   -0.6254;
    -0.1091   -0.5141    0.8507;
    -0.2206    0.7274   -0.6498;
    0.8359    0.2674    0.4794;
    0.9873    0.1103    0.1147;
    0.7471    0.0659   -0.6615;
    0.6119   -0.2508    0.7502;
    -0.6191    0.0776    0.7815;
    0.7663   -0.4739    0.4339;
    -0.5699    0.5369    0.6220;
    0.0232   -0.9989    0.0401;
    0.0671   -0.4207   -0.9047;
    -0.2145    0.5538    0.8045;
    0.8554   -0.4894    0.1698;
    -0.7912   -0.4194    0.4450;
    -0.2341    0.0754   -0.9693;
    -0.7725    0.6346   -0.0216;
    0.0228    0.7946   -0.6067;
    0.7461   -0.3966   -0.5348;
    -0.4045   -0.0837   -0.9107;
    -0.4364    0.6084   -0.6629;
    0.6177   -0.3175   -0.7195;
    -0.4301   -0.0198    0.9026;
    -0.1489   -0.9706    0.1892;
    0.0879    0.9070   -0.4117;
    -0.7764   -0.4707   -0.4190;
    0.9850    0.1352   -0.1073;
    -0.1581   -0.3154    0.9357;
    0.8938   -0.3246    0.3096;
    0.8358   -0.4464   -0.3197;
    0.4943    0.4679    0.7327;
    -0.3095    0.9015   -0.3024;
    -0.3363   -0.8942   -0.2956;
    -0.1271   -0.9274   -0.3519;
    0.3523   -0.8717   -0.3407;
    0.7188   -0.6321    0.2895;
    -0.7447    0.0924   -0.6610;
    0.1622    0.7186    0.6762;
    -0.9406   -0.0829   -0.3293;
    -0.1229    0.9204    0.3712;
    -0.8802    0.4668    0.0856;
    -0.2062   -0.1035    0.9730;
    -0.4861   -0.7586   -0.4338;
    -0.6138    0.7851    0.0827;
    0.8476    0.0504    0.5282;
    0.3236    0.4698   -0.8213;
    -0.7053   -0.6935    0.1473;
    0.1511    0.3778    0.9135;
    0.6011    0.5847    0.5448;
    0.3610    0.3183    0.8766;
    0.9432    0.3304    0.0341;
    0.2423   -0.8079   -0.5372;
    0.4431   -0.1578    0.8825;
    0.6204    0.5320   -0.5763;
    -0.2806   -0.5376   -0.7952;
    -0.5279   -0.8071    0.2646;
    -0.4214   -0.6159    0.6656;
    0.6759   -0.5995   -0.4288;
    0.5670    0.8232   -0.0295;
    -0.0874    0.4284   -0.8994;
    0.8780   -0.0192   -0.4782;
    0.0166    0.8421    0.5391;
    -0.7741    0.2931   -0.5610;
    0.9636   -0.0579   -0.2611;
    0         0   -1.0000;
    -0.5924         0   -0.8056;
    0.7191    0.1575    0.6768;
    0.9151    0.3479   -0.2040;
    -0.5535   -0.2437   -0.7964;
    0.0844   -0.9609    0.2636;
    -0.9512    0.3015   -0.0651;
    0.4225   -0.8984   -0.1202;
    -0.5916    0.6396   -0.4909;
    -0.3172   -0.8818    0.3489;
    0.1988    0.6687   -0.7164;
    0.2735   -0.3047    0.9123;
    -0.9714    0.1171   -0.2066;
    0.5215    0.4013   -0.7530;
    0.3978    0.9131    0.0897;
    -0.2680   -0.8196   -0.5063;
    0.6824    0.6532    0.3281;
    -0.4748    0.7261    0.4973;
    -0.4504    0.4036   -0.7964;
    0.5551    0.8034    0.2153;
    -0.0455    0.2169   -0.9751;
    -0.0483   -0.5845   -0.8099;
    0.1909    0.1544    0.9694;
    -0.8383   -0.5084   -0.1969;
    0.2464   -0.1148   -0.9623;
    0.7458   -0.6318   -0.2114;
    0.0080    0.9831    0.1828;
    0.2630   -0.5386    0.8005;
    0.0507   -0.6425    0.7646;
    -0.4476    0.8877   -0.1081;
    0.5627   -0.7710   -0.2982;
    0.3790   -0.7774    0.5020;
    0.6217   -0.4586    0.6350;
    0.1506   -0.8688    0.4718;
    0.4579   -0.2131   -0.8631;
    0.8349    0.2124   -0.5077;
    -0.7682    0.1732    0.6163;
    -0.0997    0.7168    0.6901;
    -0.0386    0.2146    0.9759;
    -0.9312   -0.1655    0.3249;
    -0.9151   -0.3053   -0.2634;
    -0.8081   -0.5289    0.2593;
    0.3632    0.9225   -0.1305;
    -0.2709    0.3327    0.9033;
    0.1942    0.9790    0.0623;
    -0.6302    0.7641   -0.1377;
    0.6948    0.3137   -0.6471;
    0.6596    0.6452   -0.3854;
    0.9454   -0.2713   -0.1805;
    0.2586    0.7957   -0.5477;
    0.3576   -0.6511   -0.6695;
    0.8490    0.5275   -0.0328;
    -0.3830   -0.2499    0.8893;
    -0.8804    0.2392    0.4095;
    -0.4321    0.4475    0.7829;
    0.5821    0.1656   -0.7961;
    -0.3963   -0.6637   -0.6344;
    0.7222    0.6855    0.0929;
    -0.2130    0.9650    0.1527;
    -0.4737   -0.7367    0.4825;
    0.9956   -0.0891   -0.0278;
    0.5178   -0.7899    0.3287;
    0.8906   -0.1431    0.4317;
    -0.2431    0.9670   -0.0764;
    0.6812    0.3807    0.6254;
    0.1091    0.5141   -0.8507;
    0.2206   -0.7274    0.6498;
    -0.8359   -0.2674   -0.4794;
    -0.9873   -0.1103   -0.1147;
    -0.7471   -0.0659    0.6615;
    -0.6119    0.2508   -0.7502;
    0.6191   -0.0776   -0.7815;
    -0.7663    0.4739   -0.4339;
    0.5699   -0.5369   -0.6220;
    -0.0232    0.9989   -0.0401;
    -0.0671    0.4207    0.9047;
    0.2145   -0.5538   -0.8045;
    -0.8554    0.4894   -0.1698;
    0.7912    0.4194   -0.4450;
    0.2341   -0.0754    0.9693;
    0.7725   -0.6346    0.0216;
    -0.0228   -0.7946    0.6067;
    -0.7461    0.3966    0.5348;
    0.4045    0.0837    0.9107;
    0.4364   -0.6084    0.6629;
    -0.6177    0.3175    0.7195;
    0.4301    0.0198   -0.9026;
    0.1489    0.9706   -0.1892;
    -0.0879   -0.9070    0.4117;
    0.7764    0.4707    0.4190;
    -0.9850   -0.1352    0.1073;
    0.1581    0.3154   -0.9357;
    -0.8938    0.3246   -0.3096;
    -0.8358    0.4464    0.3197;
    -0.4943   -0.4679   -0.7327;
    0.3095   -0.9015    0.3024;
    0.3363    0.8942    0.2956;
    0.1271    0.9274    0.3519;
    -0.3523    0.8717    0.3407;
    -0.7188    0.6321   -0.2895;
    0.7447   -0.0924    0.6610;
    -0.1622   -0.7186   -0.6762;
    0.9406    0.0829    0.3293;
    0.1229   -0.9204   -0.3712;
    0.8802   -0.4668   -0.0856;
    0.2062    0.1035   -0.9730;
    0.4861    0.7586    0.4338;
    0.6138   -0.7851   -0.0827;
    -0.8476   -0.0504   -0.5282;
    -0.3236   -0.4698    0.8213;
    0.7053    0.6935   -0.1473;
    -0.1511   -0.3778   -0.9135;
    -0.6011   -0.5847   -0.5448;
    -0.3610   -0.3183   -0.8766;
    -0.9432   -0.3304   -0.0341;
    -0.2423    0.8079    0.5372;
    -0.4431    0.1578   -0.8825;
    -0.6204   -0.5320    0.5763;
    0.2806    0.5376    0.7952;
    0.5279    0.8071   -0.2646;
    0.4214    0.6159   -0.6656;
    -0.6759    0.5995    0.4288;
    -0.5670   -0.8232    0.0295;
    0.0874   -0.4284    0.8994;
    -0.8780    0.0192    0.4782;
    -0.0166   -0.8421   -0.5391;
    0.7741   -0.2931    0.5610;
    -0.9636    0.0579    0.2611];
end

function dirs = radialsampling(dir, n)

% compute Equator Points

dt = 2*pi/n;
theta = 0:dt:(2*pi-dt);

dirs = [cos(theta)', sin(theta)', 0*theta']';

v = [-dir(2), dir(1), 0];
s = sqrt(sum(v.^2));
c = dir(3);
V = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
R = eye(3) + V + V*V * (1-c)/s^2;

dirs = R*dirs;

end

