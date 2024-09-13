clear;
clc;

%   完全匹配层的厚度(网格数)
PML_l = 10;
%   指定的震源的中心频率(Hz)  【地震勘探】
f0 = 25;

dx = 10;    dz = 10;
dt = 5e-4;

t0 = 3/(2*f0);

nx_main_region = 256 - 2*PML_l;     %   X方向总点数
nz_main_region = 256 - 2*PML_l;     %   Z方向总点数

t_wavefield = 0.325;

mode = "wavefield snapshot";
if mode == "wavefield snapshot"
    % t_total = 1.0;
    t_total = t_wavefield;
    nt = floor(t_wavefield/dt);
elseif mode == "synthetic seismogram"
    t_total = 4.0;
    nt = floor(t_total/dt);
end

PML_x = PML_l*dx;   PML_z = PML_l*dz;

x = (-PML_l)*dx:dx:(nx_main_region + PML_l - 1)*dx;    % x方向上的离散化，在两端预留PML层所需的宽度，以便于后续程序实现
nx = length(x);     %   x方向的网格数量
x_main_region = 0:dx:(nx_main_region - 1)*dx;    %   计算域在x方向上的长度（忽略了PML区域，保留波场模拟和合成地震记录的有效区域）

z = (-PML_l)*dz:dz:(nz_main_region + PML_l - 1)*dz;    % z方向上的离散化，在两端预留PML层所需的宽度，以便于后续程序实现
nz = length(z);     %   z方向的网格数量
z_main_region = 0:dz:(nz_main_region - 1)*dz;    %   计算域在x方向上的长度

x0 = ((nx_main_region - 1)/2)*dx;    z0 = (((nz_main_region - 1)/2))*dz;
nx0 = floor(x0/dx);     nz0 = floor(z0/dz);

%   中心频率对应的角频率(rad/s)
omega0 = 2*pi*f0;

%   自然状态下的密度(kg m/s)
Rho_N = zeros(nz_main_region,nx_main_region);           %   2.5e3;
%   等温纵波速度（自然状态）(m/s)
Vp_iso_N = zeros(nz_main_region,nx_main_region);        %   3e3;
%   横波速度（自然状态）(m/s)
Vs_N = zeros(nz_main_region,nx_main_region);            %   1.8e3;
%   线热膨胀系数（自然状态）(K^(-1))
Alpha_N = zeros(nz_main_region,nx_main_region);         %   3e-6;
%   热导率（自然状态）(m kg s^(-3) K^(-1))
K_N = zeros(nz_main_region,nx_main_region);             %   5;
%   参考绝对温度（自然状态）（摄氏度）
T_N = zeros(nz_main_region,nx_main_region);             %   300;
%   单位体积热容（自然状态）(m^2 s^(-2) K^(-1))
cV_N = zeros(nz_main_region,nx_main_region);             %   110;
%   拉梅常数（自然状态）
Mu_N = Rho_N.*Vs_N.^2;          
Lambda_iso_N = Rho_N.*Vp_iso_N.^2 - 2*Mu_N; 
%   三阶弹性常数（自然状态）
A = zeros(nz_main_region,nx_main_region);        
B = zeros(nz_main_region,nx_main_region);         
C = zeros(nz_main_region,nx_main_region);
%   围岩压力
P = zeros(nz_main_region,nx_main_region);
%   体热膨胀系数（自然状态）
Beta_N = (3*Lambda_iso_N + 2*Mu_N).*Alpha_N;
%   热弛豫时间（自然状态）
tau_N = 0;

[Vp_iso_N,Vs_N,Rho_N,K_N,Alpha_N,T_N,cV_N,Beta_N,Lambda_iso_N,Mu_N] = Add_PML_layers_CTE(nz_main_region,nx_main_region,PML_l,Vp_iso_N,Vs_N,Rho_N,K_N,Alpha_N,T_N,cV_N,Beta_N,Lambda_iso_N,Mu_N);

[nz,nx] = size(Vp_iso_N);
nz_PML = nz;    nx_PML = nx;

x0 = x0 + PML_l*dx;     z0 = z0 + PML_l*dz;
nx0 = floor(x0/dx) + 1;     nz0 = floor(z0/dz) + 1;

P = zeros(nz_PML,nx_PML);

Rho_N(1:nz_PML,1:nx_PML) = 2.5e3;
Vp_iso_N(1:nz_PML,1:nx_PML) = 3e3;
Vs_N(1:nz_PML,1:nx_PML) = 1.8e3;

Alpha_N(1:nz_PML,1:nx_PML) = 0;

K_N(1:nz_PML,1:nx_PML) = 5;
T_N(1:nz_PML,1:nx_PML) = 300;
cV_N(1:nz_PML,1:nx_PML) = 110;
A(1:nz_PML,1:nx_PML) = -1200e9;        
B(1:nz_PML,1:nx_PML) = -450e9;         
C(1:nz_PML,1:nx_PML) = -370e9;

P(1:nz_PML,1:nx_PML) = 0;

Mu_N = Rho_N.*Vs_N.^2; 
Lambda_iso_N = Rho_N.*Vp_iso_N.^2 - 2*Mu_N;  
Beta_N = (3*Lambda_iso_N + 2*Mu_N).*Alpha_N;

K_bulk_iso_N = Lambda_iso_N + (2/3).*Mu_N;
Young_E_N = Mu_N.*(3*Lambda_iso_N + 2*Mu_N)./(Lambda_iso_N + Mu_N);
Poisson_ratio_N = Lambda_iso_N./(2.*(Lambda_iso_N + Mu_N));

b_N = Beta_N.*sqrt(T_N./(Rho_N.*cV_N));
a_N = sqrt(K_N./cV_N);
Vp_adi_N = sqrt(Vp_iso_N.^2 + b_N.^2);

%   计算指定频率下的纵波速度
M0_N = (1i.*omega0.*a_N.^2)./(1 + 1i.*omega0.*tau_N);
Vc_E = sqrt((Vp_adi_N.^2 + M0_N + sqrt((Vp_adi_N.^2 + M0_N).^2 - 4.*M0_N.*Vp_iso_N.^2))/2);
Vc_T = sqrt((Vp_adi_N.^2 + M0_N - sqrt((Vp_adi_N.^2 + M0_N).^2 - 4.*M0_N.*Vp_iso_N.^2))/2);
Vp_E = 1./real(1./Vc_E);
Vp_T = 1./real(1./Vc_T);

initial_stress = 0;
thermoelasticity = "CTE";
media = "normal";

[C11_iso_N,C12_iso_N,C13_iso_N,C14_iso_N,C15_iso_N,C16_iso_N,...
 C21_iso_N,C22_iso_N,C23_iso_N,C24_iso_N,C25_iso_N,C26_iso_N,...
 C31_iso_N,C32_iso_N,C33_iso_N,C34_iso_N,C35_iso_N,C36_iso_N,...
 C41_iso_N,C42_iso_N,C43_iso_N,C44_iso_N,C45_iso_N,C46_iso_N,...
 C51_iso_N,C52_iso_N,C53_iso_N,C54_iso_N,C55_iso_N,C56_iso_N,...
 C61_iso_N,C62_iso_N,C63_iso_N,C64_iso_N,C65_iso_N,C66_iso_N] = TEI_isothermal_isotropic_elastic_coefficients(nz_PML,nx_PML,PML_l,...
                                                                                                              Vp_iso_N,Vs_N,Rho_N,...
                                                                                                              K_N,Alpha_N,T_N,cV_N,Beta_N,...
                                                                                                              Lambda_iso_N,Mu_N,A,B,C,...
                                                                                                              thermoelasticity,initial_stress,P);


C123_iso_N = 2*C;       C112_iso_N = 2*B + 2*C;   C111_iso_N = 2*C + 6*B + 2*A;
C144_iso_N = C112_iso_N - C123_iso_N;     C456_iso_N = (C111_iso_N - C123_iso_N - 6*C144_iso_N)/8;    
C155_iso_N = C144_iso_N + 2*C456_iso_N;
C255_iso_N = C144_iso_N;    C366_iso_N = C144_iso_N;    C223_iso_N = C112_iso_N;    
C133_iso_N = C112_iso_N;    C232_iso_N = C112_iso_N;    C113_iso_N = C112_iso_N;    
C122_iso_N = C112_iso_N;    C233_iso_N = C112_iso_N;    %   2*C + 2*B
C244_iso_N= C155_iso_N;    C344_iso_N = C155_iso_N;    C166_iso_N = C155_iso_N;    
C266_iso_N = C155_iso_N;    C355_iso_N = C155_iso_N;    C222_iso_N = C111_iso_N;    
C333_iso_N = C111_iso_N;  

K11_N = K_N;      K22_N = K_N;      K33_N = K_N;
K11_I = K11_N;      K22_I = K22_N;      K33_I = K33_N;    

Alpha11_N = Alpha_N;                  Alpha12_N = 0;                          Alpha13_N = 0; 
Alpha21_N = Alpha12_N;                Alpha22_N = Alpha_N;                    Alpha23_N = 0; 
Alpha31_N = Alpha13_N;                Alpha32_N = Alpha23_N;                  Alpha33_N = Alpha_N;

Alpha11_I = Alpha11_N;                Alpha12_I = Alpha12_N;                  Alpha13_I = Alpha13_N; 
Alpha21_I = Alpha12_I;                Alpha22_I = Alpha22_N;                  Alpha23_I = Alpha23_N; 
Alpha31_I = Alpha13_I;                Alpha32_I = Alpha23_I;                  Alpha33_I = Alpha33_N; 

%%%%%%%  热应力张量  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M11_N = -(C11_iso_N.*Alpha11_N + C12_iso_N.*Alpha22_N + C13_iso_N.*Alpha33_N);
M22_N = -(C12_iso_N.*Alpha11_N + C22_iso_N.*Alpha22_N + C23_iso_N.*Alpha33_N);
M33_N = -(C13_iso_N.*Alpha11_N + C23_iso_N.*Alpha22_N + C33_iso_N.*Alpha33_N);
%%%%%%% 体积热容 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_E_N = cV_N./Rho_N;
%%%%%%%  热扩散率 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a11_N = sqrt(K11_N./cV_N);
a22_N = sqrt(K22_N./cV_N);
a33_N = sqrt(K33_N./cV_N);

T_I = zeros(nz_PML,nx_PML);
T_I(1:nz_PML,1:nx_PML) = T_N;
%   初始温度增量    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dtheta_I = T_I - T_N;

exx_initial = zeros(nz_PML,nx_PML);
eyy_initial = zeros(nz_PML,nx_PML);
ezz_initial = zeros(nz_PML,nx_PML);
exy_initial = zeros(nz_PML,nx_PML);
eyz_initial = zeros(nz_PML,nx_PML);
exz_initial = zeros(nz_PML,nx_PML);

sigma_xx_initial = zeros(nz_PML,nx_PML);
sigma_yy_initial = zeros(nz_PML,nx_PML);
sigma_zz_initial = zeros(nz_PML,nx_PML);
sigma_xy_initial = zeros(nz_PML,nx_PML);
sigma_yz_initial = zeros(nz_PML,nx_PML);
sigma_xz_initial = zeros(nz_PML,nx_PML);

enn_initial = exx_initial + eyy_initial + ezz_initial;

%%%%    密度    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rho_I = Rho_N.*(1 - enn_initial);
cV_I = Rho_I.*c_E_N;
J = (1 + exx_initial + ezz_initial);

a11_I = sqrt(K11_I./cV_I);
a22_I = sqrt(K22_I./cV_I);
a33_I = sqrt(K33_I./cV_I);

%%%%%% 等效声弹性刚度系数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C11_iso_eqv_I = C11_iso_N.*(1 + 3*exx_initial - eyy_initial - ezz_initial) + C111_iso_N.*exx_initial + C112_iso_N.*eyy_initial + C113_iso_N.*ezz_initial;
C22_iso_eqv_I = C22_iso_N.*(1 + 3*eyy_initial - exx_initial - ezz_initial) + C122_iso_N.*exx_initial + C222_iso_N.*eyy_initial + C223_iso_N.*ezz_initial;
C33_iso_eqv_I = C33_iso_N.*(1 + 3*ezz_initial - exx_initial - eyy_initial) + C133_iso_N.*exx_initial + C233_iso_N.*eyy_initial + C333_iso_N.*ezz_initial;

C13_iso_eqv_I = C13_iso_N.*(1 + exx_initial - eyy_initial + ezz_initial) + C113_iso_N.*exx_initial + C123_iso_N.*eyy_initial + C133_iso_N.*ezz_initial;
C23_iso_eqv_I = C23_iso_N.*(1 - exx_initial + eyy_initial + ezz_initial) + C123_iso_N.*exx_initial + C223_iso_N.*eyy_initial + C233_iso_N.*ezz_initial;
C12_iso_eqv_I = C12_iso_N.*(1 + exx_initial + eyy_initial - ezz_initial) + C112_iso_N.*exx_initial + C122_iso_N.*eyy_initial + C123_iso_N.*ezz_initial;

C44_iso_eqv_I = C44_iso_N.*(1 - exx_initial + eyy_initial + ezz_initial) + C144_iso_N.*exx_initial + C244_iso_N.*eyy_initial + C344_iso_N.*ezz_initial;
C55_iso_eqv_I = C55_iso_N.*(1 + exx_initial - eyy_initial + ezz_initial) + C155_iso_N.*exx_initial + C255_iso_N.*eyy_initial + C355_iso_N.*ezz_initial;
C66_iso_eqv_I = C66_iso_N.*(1 + exx_initial + eyy_initial - ezz_initial) + C166_iso_N.*exx_initial + C266_iso_N.*eyy_initial + C366_iso_N.*ezz_initial;

C14_iso_eqv_I = (C13_iso_N + C12_iso_N + 2*C144_iso_N).*eyz_initial;
C15_iso_eqv_I = (2*C55_iso_N + C11_iso_N + C13_iso_N + 2*C155_iso_N).*exz_initial;
C16_iso_eqv_I = (2*C66_iso_N + C11_iso_N + C12_iso_N + 2*C166_iso_N).*exy_initial;

C34_iso_eqv_I = (2*C44_iso_N + C33_iso_N + C23_iso_N + 2*C344_iso_N).*eyz_initial;
C35_iso_eqv_I = (2*C55_iso_N + C33_iso_N + C13_iso_N + 2*C355_iso_N).*exz_initial;
C36_iso_eqv_I = (C23_iso_N + C13_iso_N + 2*C366_iso_N).*exy_initial;

C24_iso_eqv_I = (2*C44_iso_N + C22_iso_N + C23_iso_N + 2*C244_iso_N).*eyz_initial; 
C26_iso_eqv_I = (2*C66_iso_N + C22_iso_N + C12_iso_N + 2*C266_iso_N).*exy_initial; 
C25_iso_eqv_I = (C23_iso_N + C12_iso_N + 2*C255_iso_N).*exz_initial;   

C45_iso_eqv_I = (C55_iso_N + C44_iso_N + 2*C456_iso_N).*exy_initial;  
C46_iso_eqv_I = (C66_iso_N + C44_iso_N + 2*C456_iso_N).*exz_initial; 
C56_iso_eqv_I = (C66_iso_N + C55_iso_N + 2*C456_iso_N).*eyz_initial; 

M11_eqv_I = -(C11_iso_eqv_I.*Alpha11_I + C12_iso_eqv_I.*Alpha22_I + C13_iso_eqv_I.*Alpha33_I);
M22_eqv_I = -(C12_iso_eqv_I.*Alpha11_I + C22_iso_eqv_I.*Alpha22_I + C23_iso_eqv_I.*Alpha33_I);
M33_eqv_I = -(C13_iso_eqv_I.*Alpha11_I + C23_iso_eqv_I.*Alpha22_I + C33_iso_eqv_I.*Alpha33_I);

%%%%%% 声弹性波传播系数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A11_iso_I = sigma_xx_initial + C11_iso_eqv_I;
A22_iso_I = sigma_yy_initial + C22_iso_eqv_I;
A33_iso_I = sigma_zz_initial + C33_iso_eqv_I;
A12_iso_I = C12_iso_eqv_I;
A13_iso_I = C13_iso_eqv_I;
A23_iso_I = C23_iso_eqv_I;
A44_iso_I = sigma_zz_initial + C44_iso_eqv_I;
A55_iso_I = sigma_xx_initial + C55_iso_eqv_I;
A66_iso_I = sigma_yy_initial + C66_iso_eqv_I;
A77_iso_I = sigma_yy_initial + C44_iso_eqv_I;
A88_iso_I = sigma_zz_initial + C55_iso_eqv_I;
A99_iso_I = sigma_xx_initial + C66_iso_eqv_I;
A14_iso_I = C14_iso_eqv_I;
A15_iso_I = C15_iso_eqv_I;
A16_iso_I = sigma_xy_initial + C16_iso_eqv_I;
A17_iso_I = C14_iso_eqv_I;
A18_iso_I = sigma_xz_initial + C15_iso_eqv_I;
A19_iso_I = C16_iso_eqv_I;
A24_iso_I = sigma_yz_initial + C24_iso_eqv_I;
A25_iso_I = C25_iso_eqv_I;
A26_iso_I = C26_iso_eqv_I;
A27_iso_I = C24_iso_eqv_I;
A28_iso_I = C25_iso_eqv_I;
A29_iso_I = sigma_xy_initial + C26_iso_eqv_I;
A34_iso_I = C34_iso_eqv_I;
A35_iso_I = sigma_xz_initial + C35_iso_eqv_I;
A36_iso_I = C36_iso_eqv_I;
A37_iso_I = sigma_yz_initial + C34_iso_eqv_I;
A38_iso_I = C35_iso_eqv_I;
A39_iso_I = C36_iso_eqv_I;
A45_iso_I = C45_iso_eqv_I;
A46_iso_I = C46_iso_eqv_I;
A47_iso_I = C44_iso_eqv_I;
A48_iso_I = C45_iso_eqv_I;
A49_iso_I = sigma_xz_initial + C46_iso_eqv_I;
A56_iso_I = C56_iso_eqv_I;
A57_iso_I = sigma_xy_initial + C45_iso_eqv_I;
A58_iso_I = C55_iso_eqv_I;
A59_iso_I = C56_iso_eqv_I;
A67_iso_I = C46_iso_eqv_I;
A68_iso_I = sigma_yz_initial + C56_iso_eqv_I;
A69_iso_I = C66_iso_eqv_I;
A78_iso_I = C45_iso_eqv_I;
A79_iso_I = C46_iso_eqv_I;
A89_iso_I = C56_iso_eqv_I;

V11_iso_I = sqrt(A11_iso_I./Rho_I);   
V12_iso_I = sqrt(A12_iso_I./Rho_I);   
V13_iso_I = sqrt(A13_iso_I./Rho_I);   
V14_iso_I = sqrt(A14_iso_I./Rho_I);   
V15_iso_I = sqrt(A15_iso_I./Rho_I);   
V16_iso_I = sqrt(A16_iso_I./Rho_I);   
V17_iso_I = sqrt(A17_iso_I./Rho_I);   
V18_iso_I = sqrt(A18_iso_I./Rho_I);   
V19_iso_I = sqrt(A11_iso_I./Rho_I);
                                                          
V22_iso_I = sqrt(A22_iso_I./Rho_I);   
V23_iso_I = sqrt(A23_iso_I./Rho_I);   
V24_iso_I = sqrt(A24_iso_I./Rho_I);   
V25_iso_I = sqrt(A25_iso_I./Rho_I);   
V26_iso_I = sqrt(A26_iso_I./Rho_I);   
V27_iso_I = sqrt(A27_iso_I./Rho_I);   
V28_iso_I = sqrt(A28_iso_I./Rho_I);   
V29_iso_I = sqrt(A29_iso_I./Rho_I);
                                                                                                                    
V33_iso_I = sqrt(A33_iso_I./Rho_I);   
V34_iso_I = sqrt(A34_iso_I./Rho_I);   
V35_iso_I = sqrt(A35_iso_I./Rho_I);   
V36_iso_I = sqrt(A36_iso_I./Rho_I);   
V37_iso_I = sqrt(A37_iso_I./Rho_I);   
V38_iso_I = sqrt(A38_iso_I./Rho_I);   
V39_iso_I = sqrt(A39_iso_I./Rho_I);
                                                                                                                                                                              
V44_iso_I = sqrt(A44_iso_I./Rho_I);   
V45_iso_I = sqrt(A45_iso_I./Rho_I);   
V46_iso_I = sqrt(A46_iso_I./Rho_I);   
V47_iso_I = sqrt(A47_iso_I./Rho_I);   
V48_iso_I = sqrt(A48_iso_I./Rho_I);   
V49_iso_I = sqrt(A49_iso_I./Rho_I);

V55_iso_I = sqrt(A55_iso_I./Rho_I);
V56_iso_I = sqrt(A56_iso_I./Rho_I);   
V57_iso_I = sqrt(A57_iso_I./Rho_I);   
V58_iso_I = sqrt(A58_iso_I./Rho_I);   
V59_iso_I = sqrt(A59_iso_I./Rho_I);

V66_iso_I = sqrt(A66_iso_I./Rho_I);   
V67_iso_I = sqrt(A67_iso_I./Rho_I);   
V68_iso_I = sqrt(A68_iso_I./Rho_I);   
V69_iso_I = sqrt(A69_iso_I./Rho_I);
V77_iso_I = sqrt(A77_iso_I./Rho_I);   
V78_iso_I = sqrt(A78_iso_I./Rho_I);   
V79_iso_I = sqrt(A79_iso_I./Rho_I);
V88_iso_I = sqrt(A88_iso_I./Rho_I);   
V89_iso_I = sqrt(A89_iso_I./Rho_I);
V99_iso_I = sqrt(A99_iso_I./Rho_I);

C11_adi_N = C11_iso_N + (M11_N.^2.*T_N./cV_N);
C12_adi_N = C12_iso_N + (M11_N.*M22_N.*T_N./cV_N);
C13_adi_N = C13_iso_N + (M11_N.*M33_N.*T_N./cV_N);
C22_adi_N = C22_iso_N + (M22_N.^2.*T_N./cV_N);
C23_adi_N = C23_iso_N + (M22_N.*M33_N.*T_N./cV_N);
C33_adi_N = C33_iso_N + (M33_N.^2.*T_N./cV_N);
    
C11_adi_eqv_I = C11_iso_eqv_I + (M11_eqv_I.^2.*T_I./cV_I);
C12_adi_eqv_I = C12_iso_eqv_I + (M11_eqv_I.*M22_eqv_I.*T_I./cV_I);
C13_adi_eqv_I = C13_iso_eqv_I + (M11_eqv_I.*M33_eqv_I.*T_I./cV_I);
C22_adi_eqv_I = C22_iso_eqv_I + (M22_eqv_I.^2.*T_I./cV_I);
C23_adi_eqv_I = C23_iso_eqv_I + (M22_eqv_I.*M33_eqv_I.*T_I./cV_I);
C33_adi_eqv_I = C33_iso_eqv_I + (M33_eqv_I.^2.*T_I./cV_I);

%%  热声弹性理论中的绝热系数  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A11_adi_I = sigma_xx_initial + C11_adi_eqv_I;
A22_adi_I = sigma_yy_initial + C22_adi_eqv_I;
A33_adi_I = sigma_zz_initial + C33_adi_eqv_I;
A12_adi_I = C12_adi_eqv_I;
A13_adi_I = C13_adi_eqv_I;
A23_adi_I = C23_adi_eqv_I;
A44_iso_I = sigma_zz_initial + C44_iso_eqv_I;
A55_iso_I = sigma_xx_initial + C55_iso_eqv_I;
A66_iso_I = sigma_yy_initial + C66_iso_eqv_I;
A77_iso_I = sigma_yy_initial + C44_iso_eqv_I;
A88_iso_I = sigma_zz_initial + C55_iso_eqv_I;
A99_iso_I = sigma_xx_initial + C66_iso_eqv_I;
A14_iso_I = C14_iso_eqv_I;
A15_iso_I = C15_iso_eqv_I;
A16_iso_I = sigma_xy_initial + C16_iso_eqv_I;
A17_iso_I = C14_iso_eqv_I;
A18_iso_I = sigma_xz_initial + C15_iso_eqv_I;
A19_iso_I = C16_iso_eqv_I;
A24_iso_I = sigma_yz_initial + C24_iso_eqv_I;
A25_iso_I = C25_iso_eqv_I;
A26_iso_I = C26_iso_eqv_I;
A27_iso_I = C24_iso_eqv_I;
A28_iso_I = C25_iso_eqv_I;
A29_iso_I = sigma_xy_initial + C26_iso_eqv_I;
A34_iso_I = C34_iso_eqv_I;
A35_iso_I = sigma_xz_initial + C35_iso_eqv_I;
A36_iso_I = C36_iso_eqv_I;
A37_iso_I = sigma_yz_initial + C34_iso_eqv_I;
A38_iso_I = C35_iso_eqv_I;
A39_iso_I = C36_iso_eqv_I;
A45_iso_I = C45_iso_eqv_I;
A46_iso_I = C46_iso_eqv_I;
A47_iso_I = C44_iso_eqv_I;
A48_iso_I = C45_iso_eqv_I;
A49_iso_I = sigma_xz_initial + C46_iso_eqv_I;
A56_iso_I = C56_iso_eqv_I;
A57_iso_I = sigma_xy_initial + C45_iso_eqv_I;
A58_iso_I = C55_iso_eqv_I;
A59_iso_I = C56_iso_eqv_I;
A67_iso_I = C46_iso_eqv_I;
A68_iso_I = sigma_yz_initial + C56_iso_eqv_I;
A69_iso_I = C66_iso_eqv_I;
A78_iso_I = C45_iso_eqv_I;
A79_iso_I = C46_iso_eqv_I;
A89_iso_I = C56_iso_eqv_I;

%%%%%%%% 热弛豫时间 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau11_N = K11_N ./(cV_N.*Vp_iso_N.^2);          %   0.0051
tau22_N = K22_N ./(cV_N.*Vp_iso_N.^2);
tau33_N = K33_N ./(cV_N.*Vp_iso_N.^2);

tau_ave_N = (tau11_N + tau22_N + tau33_N)./3;

tau11_I = K11_I ./(cV_I.*V11_iso_I.^2);          %   0.0051
tau22_I = K22_I ./(cV_I.*V22_iso_I.^2);
tau33_I = K33_I ./(cV_I.*V33_iso_I.^2);

tau_ave_I = (tau11_I + tau22_I + tau33_I)./3;

Factor_TEI_stability = (pi^2*K_N)./((dx)^2.*tau_ave_N.*c_E_N);

%   傅里叶伪谱法因子
dfx = 1/(nx*dx);        dfz = 1/(nz*dz);

%%%%%%%%%%   场量的预定义数组  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   应力分量
sigma_xx_last = zeros(nz_PML,nx_PML);      sigma_xx_present = zeros(nz_PML,nx_PML);       sigma_xx_next = zeros(nz_PML,nx_PML);
sigma_zz_last = zeros(nz_PML,nx_PML);      sigma_zz_present = zeros(nz_PML,nx_PML);       sigma_zz_next = zeros(nz_PML,nx_PML);
sigma_xz_last = zeros(nz_PML,nx_PML);      sigma_xz_present = zeros(nz_PML,nx_PML);       sigma_xz_next = zeros(nz_PML,nx_PML);

T_PK1_xx_last = zeros(nz_PML,nx_PML);       T_PK1_xx_present = zeros(nz_PML,nx_PML);        T_PK1_xx_next = zeros(nz_PML,nx_PML);
T_PK1_zz_last = zeros(nz_PML,nx_PML);       T_PK1_zz_present = zeros(nz_PML,nx_PML);        T_PK1_zz_next = zeros(nz_PML,nx_PML);
T_PK1_xz_last = zeros(nz_PML,nx_PML);       T_PK1_xz_present = zeros(nz_PML,nx_PML);        T_PK1_xz_next = zeros(nz_PML,nx_PML);
T_PK1_zx_last = zeros(nz_PML,nx_PML);       T_PK1_zx_present = zeros(nz_PML,nx_PML);        T_PK1_zx_next = zeros(nz_PML,nx_PML);

T_PK1_xx_stage1_last = zeros(nz_PML,nx_PML);       T_PK1_xx_stage1_present = zeros(nz_PML,nx_PML);        T_PK1_xx_stage1_next = zeros(nz_PML,nx_PML);
T_PK1_zz_stage1_last = zeros(nz_PML,nx_PML);       T_PK1_zz_stage1_present = zeros(nz_PML,nx_PML);        T_PK1_zz_stage1_next = zeros(nz_PML,nx_PML);
T_PK1_xz_stage1_last = zeros(nz_PML,nx_PML);       T_PK1_xz_stage1_present = zeros(nz_PML,nx_PML);        T_PK1_xz_stage1_next = zeros(nz_PML,nx_PML);
T_PK1_zx_stage1_last = zeros(nz_PML,nx_PML);       T_PK1_zx_stage1_present = zeros(nz_PML,nx_PML);        T_PK1_zx_stage1_next = zeros(nz_PML,nx_PML);

T_PK1_xx_stage2_last = zeros(nz_PML,nx_PML);       T_PK1_xx_stage2_present = zeros(nz_PML,nx_PML);        T_PK1_xx_stage2_next = zeros(nz_PML,nx_PML);
T_PK1_zz_stage2_last = zeros(nz_PML,nx_PML);       T_PK1_zz_stage2_present = zeros(nz_PML,nx_PML);        T_PK1_zz_stage2_next = zeros(nz_PML,nx_PML);
T_PK1_xz_stage2_last = zeros(nz_PML,nx_PML);       T_PK1_xz_stage2_present = zeros(nz_PML,nx_PML);        T_PK1_xz_stage2_next = zeros(nz_PML,nx_PML);
T_PK1_zx_stage2_last = zeros(nz_PML,nx_PML);       T_PK1_zx_stage2_present = zeros(nz_PML,nx_PML);        T_PK1_zx_stage2_next = zeros(nz_PML,nx_PML);

T_PK1_xx_stage3_last = zeros(nz_PML,nx_PML);       T_PK1_xx_stage3_present = zeros(nz_PML,nx_PML);        T_PK1_xx_stage3_next = zeros(nz_PML,nx_PML);
T_PK1_zz_stage3_last = zeros(nz_PML,nx_PML);       T_PK1_zz_stage3_present = zeros(nz_PML,nx_PML);        T_PK1_zz_stage3_next = zeros(nz_PML,nx_PML);
T_PK1_xz_stage3_last = zeros(nz_PML,nx_PML);       T_PK1_xz_stage3_present = zeros(nz_PML,nx_PML);        T_PK1_xz_stage3_next = zeros(nz_PML,nx_PML);
T_PK1_zx_stage3_last = zeros(nz_PML,nx_PML);       T_PK1_zx_stage3_present = zeros(nz_PML,nx_PML);        T_PK1_zx_stage3_next = zeros(nz_PML,nx_PML);

T_PK1_xx_stage4_last = zeros(nz_PML,nx_PML);       T_PK1_xx_stage4_present = zeros(nz_PML,nx_PML);        T_PK1_xx_stage4_next = zeros(nz_PML,nx_PML);
T_PK1_zz_stage4_last = zeros(nz_PML,nx_PML);       T_PK1_zz_stage4_present = zeros(nz_PML,nx_PML);        T_PK1_zz_stage4_next = zeros(nz_PML,nx_PML);
T_PK1_xz_stage4_last = zeros(nz_PML,nx_PML);       T_PK1_xz_stage4_present = zeros(nz_PML,nx_PML);        T_PK1_xz_stage4_next = zeros(nz_PML,nx_PML);
T_PK1_zx_stage4_last = zeros(nz_PML,nx_PML);       T_PK1_zx_stage4_present = zeros(nz_PML,nx_PML);        T_PK1_zx_stage4_next = zeros(nz_PML,nx_PML);

theta_stage1_present = zeros(nz_PML,nx_PML);        theta_stage1_next = zeros(nz_PML,nx_PML);
theta_stage2_present = zeros(nz_PML,nx_PML);        theta_stage2_next = zeros(nz_PML,nx_PML);
theta_stage3_present = zeros(nz_PML,nx_PML);        theta_stage3_next = zeros(nz_PML,nx_PML);
theta_stage4_present = zeros(nz_PML,nx_PML);        theta_stage4_next = zeros(nz_PML,nx_PML);


%   速度分量
v_x_last_half = zeros(nz_PML,nx_PML);         v_x_next_half = zeros(nz_PML,nx_PML);   
v_z_last_half = zeros(nz_PML,nx_PML);         v_z_next_half = zeros(nz_PML,nx_PML);   

v_x_stage1_last_half = zeros(nz_PML,nx_PML);        v_x_stage1_next_half = zeros(nz_PML,nx_PML);
v_x_stage2_last_half = zeros(nz_PML,nx_PML);        v_x_stage2_next_half = zeros(nz_PML,nx_PML);
v_x_stage3_last_half = zeros(nz_PML,nx_PML);        v_x_stage3_next_half = zeros(nz_PML,nx_PML);
v_x_stage4_last_half = zeros(nz_PML,nx_PML);        v_x_stage4_next_half = zeros(nz_PML,nx_PML);

v_z_stage1_last_half = zeros(nz_PML,nx_PML);        v_z_stage1_next_half = zeros(nz_PML,nx_PML);
v_z_stage2_last_half = zeros(nz_PML,nx_PML);        v_z_stage2_next_half = zeros(nz_PML,nx_PML);
v_z_stage3_last_half = zeros(nz_PML,nx_PML);        v_z_stage3_next_half = zeros(nz_PML,nx_PML);
v_z_stage4_last_half = zeros(nz_PML,nx_PML);        v_z_stage4_next_half = zeros(nz_PML,nx_PML);



%   位移分量
u_x_last = zeros(nz_PML,nx_PML);    u_x_present = zeros(nz_PML,nx_PML);     u_x_next = zeros(nz_PML,nx_PML);
u_z_last = zeros(nz_PML,nx_PML);    u_z_present = zeros(nz_PML,nx_PML);     u_z_next = zeros(nz_PML,nx_PML);

%   温度增量
theta_last = zeros(nz_PML,nx_PML);     theta_present = zeros(nz_PML,nx_PML);      theta_next = zeros(nz_PML,nx_PML);
phi_last = zeros(nz_PML,nx_PML);         phi_present = zeros(nz_PML,nx_PML);        phi_next = zeros(nz_PML,nx_PML);

%   热流分量的空间导数
q_x_last_half = zeros(nz_PML,nx_PML);       q_x_next_half = zeros(nz_PML,nx_PML);
q_z_last_half = zeros(nz_PML,nx_PML);       q_z_next_half = zeros(nz_PML,nx_PML);

q_x_stage1_last_half = zeros(nz_PML,nx_PML);       q_x_stage1_next_half = zeros(nz_PML,nx_PML);
q_z_stage1_last_half = zeros(nz_PML,nx_PML);       q_z_stage1_next_half = zeros(nz_PML,nx_PML);
q_x_stage2_last_half = zeros(nz_PML,nx_PML);       q_x_stage2_next_half = zeros(nz_PML,nx_PML);
q_z_stage2_last_half = zeros(nz_PML,nx_PML);       q_z_stage2_next_half = zeros(nz_PML,nx_PML);
q_x_stage3_last_half = zeros(nz_PML,nx_PML);       q_x_stage3_next_half = zeros(nz_PML,nx_PML);
q_z_stage3_last_half = zeros(nz_PML,nx_PML);       q_z_stage3_next_half = zeros(nz_PML,nx_PML);
q_x_stage4_last_half = zeros(nz_PML,nx_PML);       q_x_stage4_next_half = zeros(nz_PML,nx_PML);
q_z_stage4_last_half = zeros(nz_PML,nx_PML);       q_z_stage4_next_half = zeros(nz_PML,nx_PML);

phi_last_half = zeros(nz_PML,nx_PML);         phi_next_half = zeros(nz_PML,nx_PML);
phi_ave = zeros(nz_PML,nx_PML);
%   空间导数  
dv_x_dx = zeros(nz_PML,nx_PML);     dv_x_dz = zeros(nz_PML,nx_PML);
dv_z_dx = zeros(nz_PML,nx_PML);     dv_z_dz = zeros(nz_PML,nx_PML);

dsigma_xx_dx = zeros(nz_PML,nx_PML);        dsigma_xz_dx = zeros(nz_PML,nx_PML);
dsigma_xz_dz = zeros(nz_PML,nx_PML);        dsigma_zz_dz = zeros(nz_PML,nx_PML);

dtheta_dx = zeros(nz_PML,nx_PML);       d2theta_dx2 = zeros(nz_PML,nx_PML);
dtheta_dz = zeros(nz_PML,nx_PML);       d2theta_dz2 = zeros(nz_PML,nx_PML);

%   中间场量
sigma_xx_mid = zeros(nz_PML,nx_PML);        sigma_zz_mid = zeros(nz_PML,nx_PML); 
phi_mid = zeros(nz_PML,nx_PML); 

sigma_xx_mid_Cerjan = zeros(nz_PML,nx_PML);      sigma_zz_mid_Cerjan = zeros(nz_PML,nx_PML);
phi_mid_Cerjan = zeros(nz_PML,nx_PML); 

Ricker_wavelet = zeros(nz_PML,nx_PML);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   合成地震记录与质点振动的时间记录
seismogram_20_v_x = zeros(floor(t_total/dt),nx);      seismogram_20_v_z = seismogram_20_v_x;
time_record_v_z = zeros(1,floor(t_total/dt));
time_record_theta = zeros(1,floor(t_total/dt));

%   所需的时间离散化   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt = dt:dt:t_total;
nt_total = length(tt);

%   显示在地震记录上的时间步长固定为1ms，以与地震勘探标准相匹配
dt_label = 1e-3;
%   实际时间步长dt与显示时间步长之比
dt_ratio = dt/dt_label;

%   显示时间步数
nt_label = t_total/dt_label;
tt_label = dt_label:dt_label:t_total;

%   边界条件类型
absorbing_boundary_condition = "MCFS-NPML";

%%%%%%%%%% MCFS-NPML 参数预定义 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   衰减函数类型
attenuation_type = "exp";
%   指数
n_power = 4;
%   理论反射系数的计算(Zhang and Shen,2010)
reflection_coefficient = 10.^(-((log(PML_l) - 1)./log(2)) - 3);
%   稳定性因子       %   对于高频情况要选取小一些！
P_stable = 0.02;
%   衰减函数的因子
gamma_damping = 6.5e-2;      delta_damping = 9.4e-2;
%   尺度因子
beta_x = ones(nz_PML,nx_PML);   beta_z = ones(nz_PML,nx_PML);
%   尺度因子的极值
beta_0_x = 5;       beta_0_z = 5;

P_beta = n_power;
%   频移因子  % 注意！该因子的指数要适当地选取！这和主频密切相关
eta_z = zeros(nz_PML,nx_PML);   eta_x = zeros(nz_PML,nx_PML);
P_eta = n_power;

%   频移因子的极值
eta_0 = 3;

%   衰减函数
alpha_x_PML_basic = zeros(nz_PML,nx_PML);   alpha_z_PML_basic = zeros(nz_PML,nx_PML);
alpha_x_PML_exp = zeros(nz_PML,nx_PML);     alpha_z_PML_exp = zeros(nz_PML,nx_PML);
K_main_damping_x = zeros(nz_PML,nx_PML);    K_main_damping_z = zeros(nz_PML,nx_PML);
alpha_x_PML = zeros(nz_PML,nx_PML);     alpha_z_PML = zeros(nz_PML,nx_PML);
alpha_x_PML_optimized = zeros(nz_PML,nx_PML);     alpha_z_PML_optimized = zeros(nz_PML,nx_PML);

%   CPML
b_x_CPML = zeros(nz_PML,nx_PML);        b_z_CPML = zeros(nz_PML,nx_PML);
a_x_CPML = zeros(nz_PML,nx_PML);        a_z_CPML = zeros(nz_PML,nx_PML);

%   能量
Energy_All = zeros(nz_PML,nx_PML);

% %   预定义最大速度
Vmax = 0;

% %   预定义最小速度
Vmin = V55_iso_I;

for i = 1:nz_PML
    for j = 1:nx_PML
        if Vp_E(i,j) > Vmax
            Vmax = Vp_E(i,j);
        end
        if Vp_E(i,j) < Vmin
            Vmin = Vp_E(i,j);
        end
    end
end

%   Nyquist频率
k_x_nyquist = pi/(dx);       k_z_nyquist = pi/(dz);

%   求传播矩阵特征值
Eigenvalue_max = 0;

ONES = ones(nz_PML,nx_PML);

%%%%%% 预定义MCFS-NPML场量数组 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_xx_x_MCFS_NPML_last = zeros(nz_PML,nx_PML);     sigma_xx_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);      sigma_xx_x_MCFS_NPML_next = zeros(nz_PML,nx_PML);
sigma_zz_z_MCFS_NPML_last = zeros(nz_PML,nx_PML);     sigma_zz_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);      sigma_zz_z_MCFS_NPML_next = zeros(nz_PML,nx_PML);
sigma_xz_x_MCFS_NPML_last = zeros(nz_PML,nx_PML);     sigma_xz_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);      sigma_xz_x_MCFS_NPML_next = zeros(nz_PML,nx_PML);
sigma_xz_z_MCFS_NPML_last = zeros(nz_PML,nx_PML);     sigma_xz_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);      sigma_xz_z_MCFS_NPML_next = zeros(nz_PML,nx_PML);

sigma_xx_xx_MCFS_NPML_last = zeros(nz_PML,nx_PML);      sigma_xx_xx_MCFS_NPML_present = zeros(nz_PML,nx_PML);       sigma_xx_xx_MCFS_NPML_next = zeros(nz_PML,nx_PML);
sigma_zz_zz_MCFS_NPML_last = zeros(nz_PML,nx_PML);      sigma_zz_zz_MCFS_NPML_present = zeros(nz_PML,nx_PML);       sigma_zz_zz_MCFS_NPML_next = zeros(nz_PML,nx_PML);
sigma_xz_xz_MCFS_NPML_last = zeros(nz_PML,nx_PML);      sigma_xz_xz_MCFS_NPML_present = zeros(nz_PML,nx_PML);       sigma_xz_xz_MCFS_NPML_next = zeros(nz_PML,nx_PML);


T_PK1_xx_x_MCFS_NPML_last = zeros(nz_PML,nx_PML);       T_PK1_xx_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);        T_PK1_xx_x_MCFS_NPML_next = zeros(nz_PML,nx_PML);
T_PK1_zz_z_MCFS_NPML_last = zeros(nz_PML,nx_PML);       T_PK1_zz_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);        T_PK1_zz_z_MCFS_NPML_next = zeros(nz_PML,nx_PML);
T_PK1_xz_z_MCFS_NPML_last = zeros(nz_PML,nx_PML);       T_PK1_xz_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);        T_PK1_xz_z_MCFS_NPML_next = zeros(nz_PML,nx_PML);
T_PK1_zx_x_MCFS_NPML_last = zeros(nz_PML,nx_PML);       T_PK1_zx_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);        T_PK1_zx_x_MCFS_NPML_next = zeros(nz_PML,nx_PML);

T_PK1_xx_stage1_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);     T_PK1_xx_stage1_x_MCFS_NPML_next = zeros(nz_PML,nx_PML);
T_PK1_xx_stage2_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);     T_PK1_xx_stage2_x_MCFS_NPML_next = zeros(nz_PML,nx_PML);
T_PK1_xx_stage3_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);     T_PK1_xx_stage3_x_MCFS_NPML_next = zeros(nz_PML,nx_PML);
T_PK1_xx_stage4_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);     T_PK1_xx_stage4_x_MCFS_NPML_next = zeros(nz_PML,nx_PML);

T_PK1_zz_stage1_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);     T_PK1_zz_stage1_z_MCFS_NPML_next = zeros(nz_PML,nx_PML);
T_PK1_zz_stage2_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);     T_PK1_zz_stage2_z_MCFS_NPML_next = zeros(nz_PML,nx_PML);
T_PK1_zz_stage3_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);     T_PK1_zz_stage3_z_MCFS_NPML_next = zeros(nz_PML,nx_PML);
T_PK1_zz_stage4_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);     T_PK1_zz_stage4_z_MCFS_NPML_next = zeros(nz_PML,nx_PML);

T_PK1_xz_stage1_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);     T_PK1_xz_stage1_z_MCFS_NPML_next = zeros(nz_PML,nx_PML);
T_PK1_xz_stage2_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);     T_PK1_xz_stage2_z_MCFS_NPML_next = zeros(nz_PML,nx_PML);
T_PK1_xz_stage3_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);     T_PK1_xz_stage3_z_MCFS_NPML_next = zeros(nz_PML,nx_PML);
T_PK1_xz_stage4_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);     T_PK1_xz_stage4_z_MCFS_NPML_next = zeros(nz_PML,nx_PML);

T_PK1_zx_stage1_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);     T_PK1_zx_stage1_x_MCFS_NPML_next = zeros(nz_PML,nx_PML);
T_PK1_zx_stage2_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);     T_PK1_zx_stage2_x_MCFS_NPML_next = zeros(nz_PML,nx_PML);
T_PK1_zx_stage3_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);     T_PK1_zx_stage3_x_MCFS_NPML_next = zeros(nz_PML,nx_PML);
T_PK1_zx_stage4_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);     T_PK1_zx_stage4_x_MCFS_NPML_next = zeros(nz_PML,nx_PML);

dq_x_x_MCFS_NPML_dx = zeros(nz_PML,nx_PML);                     dq_z_z_MCFS_NPML_dz = zeros(nz_PML,nx_PML);
dq_x_stage1_x_MCFS_NPML_dx = zeros(nz_PML,nx_PML);         dq_z_stage1_z_MCFS_NPML_dz = zeros(nz_PML,nx_PML);
dq_x_stage2_x_MCFS_NPML_dx = zeros(nz_PML,nx_PML);         dq_z_stage2_z_MCFS_NPML_dz = zeros(nz_PML,nx_PML);
dq_x_stage3_x_MCFS_NPML_dx = zeros(nz_PML,nx_PML);         dq_z_stage3_z_MCFS_NPML_dz = zeros(nz_PML,nx_PML);
dq_x_stage4_x_MCFS_NPML_dx = zeros(nz_PML,nx_PML);         dq_z_stage4_z_MCFS_NPML_dz = zeros(nz_PML,nx_PML);

theta_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);            theta_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);  
theta_stage1_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);            theta_stage1_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);
theta_stage2_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);            theta_stage2_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);
theta_stage3_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);            theta_stage3_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);
theta_stage4_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);            theta_stage4_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);

T_PK1_xx_xx_MCFS_NPML_last = zeros(nz_PML,nx_PML);       T_PK1_xx_xx_MCFS_NPML_present = zeros(nz_PML,nx_PML);        T_PK1_xx_xx_MCFS_NPML_next = zeros(nz_PML,nx_PML);
T_PK1_zz_zz_MCFS_NPML_last = zeros(nz_PML,nx_PML);       T_PK1_zz_zz_MCFS_NPML_present = zeros(nz_PML,nx_PML);        T_PK1_zz_zz_MCFS_NPML_next = zeros(nz_PML,nx_PML);

v_x_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       v_x_x_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
v_z_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       v_z_z_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
v_x_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       v_x_z_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
v_z_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       v_z_x_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);

v_x_stage1_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       v_z_stage1_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);
v_x_stage2_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       v_z_stage2_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);
v_x_stage3_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       v_z_stage3_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);
v_x_stage4_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       v_z_stage4_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);

v_x_stage1_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       v_z_stage1_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);
v_x_stage2_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       v_z_stage2_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);
v_x_stage3_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       v_z_stage3_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);
v_x_stage4_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       v_z_stage4_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);

q_x_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);         q_x_x_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
q_z_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);         q_z_z_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
q_x_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);         q_x_z_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
q_z_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);         q_z_x_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);

q_x_x_MCFS_NPML_ave = zeros(nz_PML,nx_PML);
q_z_z_MCFS_NPML_ave = zeros(nz_PML,nx_PML);

    
q_x_stage1_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       q_x_stage1_x_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
q_x_stage2_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       q_x_stage2_x_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
q_x_stage3_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       q_x_stage3_x_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
q_x_stage4_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       q_x_stage4_x_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);

q_z_stage1_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       q_z_stage1_z_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
q_z_stage2_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       q_z_stage2_z_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
q_z_stage3_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       q_z_stage3_z_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
q_z_stage4_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       q_z_stage4_z_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);

q_x_stage1_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       q_x_stage1_z_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
q_x_stage2_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       q_x_stage2_z_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
q_x_stage3_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       q_x_stage3_z_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
q_x_stage4_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       q_x_stage4_z_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);

q_z_stage1_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       q_z_stage1_x_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
q_z_stage2_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       q_z_stage2_x_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
q_z_stage3_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       q_z_stage3_x_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);
q_z_stage4_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);       q_z_stage4_x_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);

q_x_mid = zeros(nz_PML,nx_PML);  
q_z_mid = zeros(nz_PML,nx_PML);  

dv_x_x_MCFS_NPML_dx = zeros(nz_PML,nx_PML);     dv_z_z_MCFS_NPML_dz = zeros(nz_PML,nx_PML);
dv_x_z_MCFS_NPML_dz = zeros(nz_PML,nx_PML);     dv_z_x_MCFS_NPML_dx = zeros(nz_PML,nx_PML);

theta_xx_MCFS_NPML_last = zeros(nz_PML,nx_PML);     theta_xx_MCFS_NPML_present = zeros(nz_PML,nx_PML);      theta_xx_MCFS_NPML_next = zeros(nz_PML,nx_PML);
theta_zz_MCFS_NPML_last = zeros(nz_PML,nx_PML);     theta_zz_MCFS_NPML_present = zeros(nz_PML,nx_PML);      theta_zz_MCFS_NPML_next = zeros(nz_PML,nx_PML);

phi_x_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);     phi_x_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);   
phi_z_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);     phi_z_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);  

value_dv_x_x_MCFS_NPML_dx = zeros(nz_PML,nx_PML);       value_dv_z_z_MCFS_NPML_dz = zeros(nz_PML,nx_PML);
value_dv_x_z_MCFS_NPML_dz = zeros(nz_PML,nx_PML);       value_dv_z_x_MCFS_NPML_dx = zeros(nz_PML,nx_PML);

dsigma_xx_x_MCFS_NPML_dx_next_half = zeros(nz_PML,nx_PML);        dsigma_zz_z_MCFS_NPML_dz_next_half = zeros(nz_PML,nx_PML);
dsigma_xz_x_MCFS_NPML_dx_next_half = zeros(nz_PML,nx_PML);        dsigma_xz_z_MCFS_NPML_dz_next_half = zeros(nz_PML,nx_PML);

dsigma_xx_x_MCFS_NPML_dx_last_half = zeros(nz_PML,nx_PML);        dsigma_zz_z_MCFS_NPML_dz_last_half = zeros(nz_PML,nx_PML);
dsigma_xz_x_MCFS_NPML_dx_last_half = zeros(nz_PML,nx_PML);        dsigma_xz_z_MCFS_NPML_dz_last_half = zeros(nz_PML,nx_PML);

dsigma_xx_xx_MCFS_NPML_dx = zeros(nz_PML,nx_PML);        dsigma_zz_zz_MCFS_NPML_dz = zeros(nz_PML,nx_PML);
dsigma_xz_xz_MCFS_NPML_dx = zeros(nz_PML,nx_PML);        dsigma_xz_xz_MCFS_NPML_dz = zeros(nz_PML,nx_PML);

dT_PK1_xx_x_MCFS_NPML_dx = zeros(nz_PML,nx_PML);       dT_PK1_zz_z_MCFS_NPML_dz = zeros(nz_PML,nx_PML);
dT_PK1_zx_x_MCFS_NPML_dx = zeros(nz_PML,nx_PML);       dT_PK1_xz_z_MCFS_NPML_dz = zeros(nz_PML,nx_PML);
d2T_PK1_xx_xx_MCFS_NPML_dx2 = zeros(nz_PML,nx_PML);                 d2T_PK1_zz_zz_MCFS_NPML_dz2 = zeros(nz_PML,nx_PML);
d2T_PK1_zx_xz_MCFS_NPML_dx_dz = zeros(nz_PML,nx_PML);              d2T_PK1_xz_zx_MCFS_NPML_dz_dx = zeros(nz_PML,nx_PML);

dPsi_x_x_MCFS_NPML_dx = zeros(nz_PML,nx_PML);       dPsi_z_z_MCFS_NPML_dz = zeros(nz_PML,nx_PML);

dtheta_x_MCFS_NPML_dx = zeros(nz_PML,nx_PML);       dtheta_z_MCFS_NPML_dz = zeros(nz_PML,nx_PML);

d2theta_xx_MCFS_NPML_dx2 = zeros(nz_PML,nx_PML);     d2theta_zz_MCFS_NPML_dz2 = zeros(nz_PML,nx_PML);

u_x_x_MCFS_NPML_last = zeros(nz_PML,nx_PML);        u_x_x_MCFS_NPML_present = zeros(nz_PML,nx_PML);         u_x_x_MCFS_NPML_next = zeros(nz_PML,nx_PML);
u_z_z_MCFS_NPML_last = zeros(nz_PML,nx_PML);        u_z_z_MCFS_NPML_present = zeros(nz_PML,nx_PML);         u_z_z_MCFS_NPML_next = zeros(nz_PML,nx_PML);

dtheta_x_dx_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);      dtheta_x_dx_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);        
dtheta_z_dz_MCFS_NPML_last_half = zeros(nz_PML,nx_PML);      dtheta_z_dz_MCFS_NPML_next_half = zeros(nz_PML,nx_PML);     

dtheta_x_dx_last_half = zeros(nz_PML,nx_PML);        dtheta_z_dz_last_half = zeros(nz_PML,nx_PML);
dtheta_x_dx_next_half = zeros(nz_PML,nx_PML);        dtheta_z_dz_next_half = zeros(nz_PML,nx_PML);

%   C-N格式
v_z_z_MCFS_NPML_ave = zeros(nz_PML,nx_PML);         v_z_x_MCFS_NPML_ave = zeros(nz_PML,nx_PML);
v_x_x_MCFS_NPML_ave = zeros(nz_PML,nx_PML);         v_x_z_MCFS_NPML_ave = zeros(nz_PML,nx_PML);

%   能量
phi_Helmholtz_present = zeros(nz_PML,nx_PML);
phi_Helmholtz_next = zeros(nz_PML,nx_PML);

%%%%%%% 预定义时间记录 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_record_v_z_fz = zeros(500,3);
time_record_v_x_fz = zeros(500,3);
time_record_q_z_fz = zeros(500,3);
time_record_q_x_fz = zeros(500,3);
time_record_theta_fz = zeros(500,3);

time_record_v_z_R = zeros(nt,3);
% time_record_theta_fz = zeros(nt,3);
time_record_theta_R = zeros(nt,3);

Phi_MID_SOURCE = zeros(1,nt);
Phi_MID_BOUNDARY = zeros(1,nt);
dT_PRESENT_SOURCE = zeros(1,nt);
dT_PRESENT_BOUNDARY = zeros(1,nt);
V_Z_LAST_HALF_SOURCE = zeros(1,nt);
V_Z_LAST_HALF_BOUNDARY = zeros(1,nt);
V_Z_LAST_HALF_Z_MCFS_NPML = zeros(1,nt);

%%%%%%%%%%%% 计算波速 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega = 2*pi*f0;

M11_Maxwell_I = (1i.*omega.*a11_I.^2)./(-1+1i.*omega*(tau11_I));
M22_Maxwell_I = (1i.*omega.*a22_I.^2)./(-1+1i.*omega*(tau22_I));
M33_Maxwell_I = (1i.*omega.*a33_I.^2)./(-1+1i.*omega*(tau33_I));

b11_eqv_I = sqrt(M11_eqv_I.^2.*T_I./(Rho_I.*cV_I));
b22_eqv_I = sqrt(M22_eqv_I.^2.*T_I./(Rho_I.*cV_I));
b33_eqv_I = sqrt(M33_eqv_I.^2.*T_I./(Rho_I.*cV_I));

V11_qPE_adi_I = sqrt(V11_iso_I.^2 + b11_eqv_I.^2);
V22_qPE_adi_I = sqrt(V22_iso_I.^2 + b22_eqv_I.^2);
V33_qPE_adi_I = sqrt(V33_iso_I.^2 + b33_eqv_I.^2);

V11_complex_elastic_I = sqrt((V11_qPE_adi_I.^2 + M11_Maxwell_I + sqrt((V11_qPE_adi_I.^2 + M11_Maxwell_I).^2 - 4.*M11_Maxwell_I.*V11_iso_I.^2))/2);
V22_complex_elastic_I = sqrt((V22_qPE_adi_I.^2 + M22_Maxwell_I + sqrt((V22_qPE_adi_I.^2 + M22_Maxwell_I).^2 - 4.*M22_Maxwell_I.*V22_iso_I.^2))/2);
V33_complex_elastic_I = sqrt((V33_qPE_adi_I.^2 + M33_Maxwell_I + sqrt((V33_qPE_adi_I.^2 + M33_Maxwell_I).^2 - 4.*M33_Maxwell_I.*V33_iso_I.^2))/2);

V11_complex_thermal_I = sqrt((V11_qPE_adi_I.^2 + M11_Maxwell_I - sqrt((V11_qPE_adi_I.^2 + M11_Maxwell_I).^2 - 4.*M11_Maxwell_I.*V11_iso_I.^2))/2);
V22_complex_thermal_I = sqrt((V22_qPE_adi_I.^2 + M22_Maxwell_I - sqrt((V22_qPE_adi_I.^2 + M22_Maxwell_I).^2 - 4.*M22_Maxwell_I.*V22_iso_I.^2))/2);
V33_complex_thermal_I = sqrt((V33_qPE_adi_I.^2 + M33_Maxwell_I - sqrt((V33_qPE_adi_I.^2 + M33_Maxwell_I).^2 - 4.*M33_Maxwell_I.*V33_iso_I.^2))/2);

V11_phase_elastic_I = 1./real(1./V11_complex_elastic_I);
V22_phase_elastic_I = 1./real(1./V22_complex_elastic_I);
V33_phase_elastic_I = 1./real(1./V33_complex_elastic_I);

V11_phase_thermal_I = 1./real(1./V11_complex_thermal_I);
V22_phase_thermal_I = 1./real(1./V22_complex_thermal_I);
V33_phase_thermal_I = 1./real(1./V33_complex_thermal_I);

trace = 50;         %   50个地震道
d_trace = 10000/trace;          %   trace_length

nt_ricker = 512;

t_wavelet = (-round(nt_ricker/2) + 1)*dt:dt:round(nt_ricker/2)*dt;

R_wavelet = Ricker_Wavelet(dt,f0,t0,nt_ricker);
R_wavelet = R_wavelet';

size_ricker_wavelet = length(R_wavelet);

previous_two_way_travelling_time = zeros(nz_PML - 1,trace);

Reflection = zeros(nz_PML - 2*PML_l - 1,trace);

two_way_travelling_time = zeros(nz_PML - 2*PML_l,trace);
one_way_travelling_time = zeros(nz_PML - 2*PML_l,trace);

total_travelling_time = 0;

delta_t = dz./V33_phase_elastic_I;

sum_delta_t = zeros(nz_PML,nx_PML);
V_qPE_ave = zeros(nz_PML,nx_PML);

sum_dz = zeros(nz_PML,nx_PML);

for i = PML_l + 1:nz_PML - PML_l - 1
    for j = PML_l + 1:nx_PML - PML_l - 1
        sum_dz(i,j) = (i - PML_l).*dz;
        sum_delta_t(i,j) = sum_delta_t(i - 1,j) + dz./V33_phase_elastic_I(i,j);
        V_qPE_ave(i,j) = sum_dz(i,j)./sum_delta_t(i,j);
    end
end


% convolution_model_synthetic_seismogram_flag = "true";
convolution_model_synthetic_seismogram_flag = "false";

if convolution_model_synthetic_seismogram_flag == "true"

    %   地表阵列
    for k = 1:trace
        frac1 = k/(trace);
        progressbar(frac1);   
        fprintf('k = %d\n',k);
        ix_trace = ((k - 1).*d_trace)./dx + PML_l + 1;
    
        for i = PML_l + 1:nz_PML - PML_l
            if i > PML_l && i < nz_PML - PML_l 
                two_way_travelling_time(i - PML_l,k) = sqrt((abs(k - 25 - 1)*d_trace).^2 + 4*((i - PML_l)*dz)^2)./(V_qPE_ave(i,ix_trace).*dt);    %     ms
                one_way_travelling_time(i - PML_l,k) = two_way_travelling_time(i - PML_l,k)/2;
    
                Reflection(i - PML_l,k) = (Rho_I(i + 1,ix_trace)*V33_phase_elastic_I(i + 1,ix_trace) - Rho_I(i,ix_trace)*V33_phase_elastic_I(i,ix_trace))./ ...
                                          (Rho_I(i + 1,ix_trace)*V33_phase_elastic_I(i + 1,ix_trace) + Rho_I(i,ix_trace)*V33_phase_elastic_I(i,ix_trace));
    
                t_Reflection(round(two_way_travelling_time(i - PML_l,k)),k) = Reflection(i - PML_l,k); 
    
            end
        end
        synthetic(:,k) = conv(t_Reflection(:,k),R_wavelet);            %   褶积          %   256391
        LENGTH_SYN = size_ricker_wavelet/2:1:length(synthetic(:,k));
        LENGTH_SYN_L = LENGTH_SYN - size_ricker_wavelet/2 + 1;
        syn_l = i + 1 - nt_ricker/2;
        synthetic_l(LENGTH_SYN_L,k) = synthetic(LENGTH_SYN_L,k);      %   振幅最大值对应反射系数值
        for i = nt_ricker/2:1:(length(LENGTH_SYN) - nt_ricker/2)
            synthetic_n(i + 1 - nt_ricker/2,k) = synthetic_l(i,k);
        end
    end
    
    figure(2)
    wigb(synthetic_n,1);
    set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    xlabel('trace');
    ylabel('t(ms)');
    title('synthetic','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    xlim([0 51]);
    
    %   DAS VSP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   水平偏移距
    x_migration_offset = 7000;    %  震源位置为(0m,0m)，这样可以尽可能多地接收向上反射波                   
    nx_migration_offset = x_migration_offset/dx + 1;
    %   DAS空间采样间隔可以定为20m
    dz_trace_DAS_VSP = 10*dz; 
    
    nz_DAS_VSP = 7000/dz_trace_DAS_VSP + 1;         %   50个地震道
    nx_DAS_VSP = 10000/dx + 1;
    
    nt_ricker_DAS_VSP = 512;
    
    R_wavelet_DAS_VSP = Ricker_Wavelet(dt,f0,t0,nt_ricker_DAS_VSP);
    R_wavelet_DAS_VSP = R_wavelet_DAS_VSP';
    
    size_ricker_wavelet_DAS_VSP = length(R_wavelet_DAS_VSP);
    RR(1,1,:) = R_wavelet_DAS_VSP(:);
    
    total_travelling_time_DAS_VSP = 0;
    %   各层旅行时
    sum_dz_DAS_VSP      = zeros(nz_PML - 2*PML_l - 1,nx_PML - 2*PML_l);
    sum_delta_t_DAS_VSP = zeros(nz_PML - 2*PML_l - 1,nx_PML - 2*PML_l);
    delta_t_DAS_VSP = zeros(nz_PML - 2*PML_l - 1,nx_PML - 2*PML_l);
    V_qPE_ave_DAS_VSP = zeros(nz_PML - 2*PML_l,nx_PML - 2*PML_l);
    iz_trace_DAS_VSP = zeros(nz_PML - 2*PML_l,nx_PML - 2*PML_l);
    
    for i = PML_l + 1:nz_PML - PML_l
        for j = PML_l + 1:nx_PML - PML_l    %       横向不用刻意调整
            kk = i - PML_l;      ll = j - PML_l;    
            sum_dz_DAS_VSP(kk,ll) = (i - PML_l).*dz;
            delta_t_DAS_VSP(kk,ll) = dz./V33_phase_elastic_I(i,j);
            sum_delta_t_DAS_VSP(kk,ll) = sum_delta_t_DAS_VSP(kk,ll) + delta_t_DAS_VSP(kk,ll);
            if kk <= nz_PML - 2*PML_l - 1
                sum_delta_t_DAS_VSP(kk + 1,ll) = sum_delta_t_DAS_VSP(kk,ll);
                iz_trace_DAS_VSP(kk + 1,ll) = (kk + 1)*dz./dz_trace_DAS_VSP;
            end
            iz_trace_DAS_VSP(kk,ll) = (kk - 1)*dz./dz_trace_DAS_VSP + 1;
            if kk <= nz_PML - 2*PML_l - 1
                iz_trace_DAS_VSP(kk + 1,ll) = (kk + 1)*dz./dz_trace_DAS_VSP;
            end
            V_qPE_ave_DAS_VSP(kk,ll) = sum_dz_DAS_VSP(kk,ll)./sum_delta_t_DAS_VSP(kk,ll);
        end 
    end
    
    iz_PML = (1 + PML_l:nz_PML - PML_l);
    ix_PML = (1 + PML_l:nx_PML - PML_l);
    iz = (1:nz_PML - 2*PML_l);
    ix = (1:nx_PML - 2*PML_l);
    
    V_qPE_ave_DAS_VSP_temp = zeros(nz_PML - 2*PML_l,nx_PML - 2*PML_l);
    V_qPE_ave_DAS_VSP_temp(2:nz_PML - 2*PML_l,ix) = V_qPE_ave_DAS_VSP(1:nz_PML - 2*PML_l - 1,ix);
    V_qPE_ave_DAS_VSP = V_qPE_ave_DAS_VSP_temp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   DAS-VSP阵列
    nz_trace_DAS_VSP = 1:(nz_DAS_VSP - 1);    
    Z_CRP = zeros(nz_DAS_VSP,nx_PML - 2*PML_l);
    DTT_DAS_VSP = zeros(nz_DAS_VSP - 1,nx_PML - 2*PML_l);
    Reflection_DAS_VSP = zeros(nz_DAS_VSP - 1,nx_PML - 2*PML_l);
    H_interface_DAS_VSP = zeros(nz_DAS_VSP - 1,nx_PML - 2*PML_l);
    upward_reflection_travelling_time_DAS_VSP_CRP = zeros(nz_DAS_VSP - 1,nx_PML - 2*PML_l);
    downward_travelling_time_DAS_VSP_CRP = zeros(nz_DAS_VSP - 1,nx_PML - 2*PML_l);
    
    %   索引的设置需要注意！不要以DAS-VSP的道作为索引
    RATIO_X_H = zeros(nz_DAS_VSP,nx_PML - 2*PML_l);
    
    DTT_DAS_VSP = zeros(nz_DAS_VSP,nx_PML - 2*PML_l);
    UP_RTT_DAS_VSP = zeros(nz_DAS_VSP,nx_PML - 2*PML_l);
    DOWN_RTT_DAS_VSP = zeros(nz_DAS_VSP,nx_PML - 2*PML_l);
    D_TT_DAS_VSP = zeros(nz_DAS_VSP,nx_PML - 2*PML_l);
    X_CRP = zeros(nz_DAS_VSP,nx_PML - 2*PML_l);
    Z_CRP = zeros(nz_DAS_VSP,nx_PML - 2*PML_l);
    
    for i = 1:nz_DAS_VSP     %  层深度   %    打破思维定势！多次反射的话，总深度要更大才对！
        frac1 = i/(nz_DAS_VSP);
        progressbar(frac1);   
        fprintf('kz = %d\n',i);
        %   直达波
        Z_CRP(i,ix) = (i - 1)*dz_trace_DAS_VSP;        %   DAS视深度
        X_CRP(i,ix) = (ix - 1)*dx;
        if i == 1
            V_qPE_ave_DAS_VSP(i,ix) = V33_phase_elastic_I(iz_PML(i),ix_PML);
            %   在顶部界面传输的直达波
            DTT_DAS_VSP(i,ix) = sqrt(Z_CRP(i,ix).^2 + X_CRP(i,ix).^2)./(V_qPE_ave_DAS_VSP(i,ix).*dt);
        elseif i > 1
            DTT_DAS_VSP(i,ix) = sqrt(Z_CRP(i,ix).^2 + X_CRP(i,ix).^2)./(V_qPE_ave_DAS_VSP(i,ix).*dt);
        end
        TWTT_DAS_VSP(i,ix) = sqrt(X_CRP(i,ix).^2 + 4*Z_CRP(i,ix).^2)./(V_qPE_ave_DAS_VSP(i,ix).*dt);    %     ms
        %   水平X和深度Z都是可变的，但是震源和井的位置是固定的
        Reflection_DAS_VSP(i,ix) = (Rho_I(iz_PML(10*(i - 1) + 1) + 1,ix_PML).*V33_phase_elastic_I(iz_PML(10*(i - 1) + 1) + 1,ix_PML) - ...
                                        Rho_I(iz_PML(10*(i - 1) + 1),ix_PML).*V33_phase_elastic_I(iz_PML(10*(i - 1) + 1),ix_PML))./ ...
                                       (Rho_I(iz_PML(10*(i - 1) + 1) + 1,ix_PML).*V33_phase_elastic_I(iz_PML(10*(i - 1) + 1) + 1,ix_PML) + ...
                                        Rho_I(iz_PML(10*(i - 1) + 1),ix_PML).*V33_phase_elastic_I(iz_PML(10*(i - 1) + 1),ix_PML));
        for j = 1:nx_DAS_VSP - 1
            if V33_phase_elastic_I(iz_PML,ix_PML(j)) > V33_phase_elastic_I(iz_PML,ix_PML(j) + 1)
                Reflection_DAS_VSP(i,j + 1) = (Rho_I(iz_PML,ix_PML(j) + 1).*V33_phase_elastic_I(iz_PML,ix_PML(j) + 1) - ...
                                               Rho_I(iz_PML,ix_PML(j)).*V33_phase_elastic_I(iz_PML,ix_PML(j)))./ ...
                                              (Rho_I(iz_PML,ix_PML(j) + 1).*V33_phase_elastic_I(iz_PML,ix_PML(j) + 1) + ...
                                               Rho_I(iz_PML,ix_PML(j)).*V33_phase_elastic_I(iz_PML,ix_PML(j)));
    
            elseif V33_phase_elastic_I(iz_PML,ix_PML(j)) < V33_phase_elastic_I(iz_PML,ix_PML(j) + 1)
                Reflection_DAS_VSP(i,j) = (Rho_I(iz_PML,ix_PML(j) + 1).*V33_phase_elastic_I(iz_PML,ix_PML(j) + 1) - ...
                                           Rho_I(iz_PML,ix_PML(j)).*V33_phase_elastic_I(iz_PML,ix_PML(j)))./ ...
                                          (Rho_I(iz_PML,ix_PML(j) + 1).*V33_phase_elastic_I(iz_PML,ix_PML(j) + 1) + ...
                                           Rho_I(iz_PML,ix_PML(j)).*V33_phase_elastic_I(iz_PML,ix_PML(j)));
            end
        end
    end
    
    for i = 1:nz_DAS_VSP
        for j = 1:nx_DAS_VSP
            H_DTT_DAS_VSP(i,j) = Z_CRP(i,j);
            if j == nx_migration_offset
                D_TT_DAS_VSP(i,j) = sqrt(x_migration_offset.^2 + H_DTT_DAS_VSP(i,j).^2)./(V_qPE_ave_DAS_VSP(i,j).*dt);
                t_DTT_DAS_VSP(i,round(D_TT_DAS_VSP(i,j))) = Reflection_DAS_VSP(i,j);
            end
        end
        CONV_DTT_DAS_VSP = conv(t_DTT_DAS_VSP(i,:),R_wavelet_DAS_VSP);       % 1835 + 512
        SYN_DTT_DAS_VSP(i,1:length(CONV_DTT_DAS_VSP)) = CONV_DTT_DAS_VSP;
        LENGTH_SYN_DTT_DAS_VSP = size_ricker_wavelet_DAS_VSP/2:1:length(SYN_DTT_DAS_VSP(i,:)) + size_ricker_wavelet_DAS_VSP/2 - 1;
        for j = nt_ricker_DAS_VSP/2:1:(length(LENGTH_SYN_DTT_DAS_VSP))
            SYN_DTT_DAS_VSP_n(i,j + 1 - nt_ricker_DAS_VSP/2) = SYN_DTT_DAS_VSP(i,j);
        end
    end
    
    figure(111)
    wigb_VSP(SYN_DTT_DAS_VSP_n,1);
    set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    xlabel('t(ms)');
    ylabel('trace');
    title('synthetic direct travelling DAS VSP','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   考验几何能力——几何地震学
    X_UP_surface = zeros(nz_DAS_VSP,nx_DAS_VSP);    %   地表反射接收点水平距离
    DX_RHS_surface = zeros(nz_DAS_VSP,nx_DAS_VSP);
    X_CRP = zeros(nz_DAS_VSP,nx_DAS_VSP);
    DX_LHS_surface = zeros(nz_DAS_VSP,nx_DAS_VSP);
    DX_LHS_CRP = zeros(nz_DAS_VSP,nx_DAS_VSP);
    Z_UP_RTT_DAS_VSP = zeros(nz_DAS_VSP,nx_DAS_VSP);
    % RATIO_UP_Z_X = zeros(nz_DAS_VSP,nx_DAS_VSP);
    DZ_add = zeros(nz_DAS_VSP,nx_DAS_VSP);
    H_UP = zeros(nz_DAS_VSP,nx_DAS_VSP);
    H_UP_REAL = zeros(nz_DAS_VSP,nx_DAS_VSP);
    
    kx1 = zeros(nz_DAS_VSP,nx_DAS_VSP);
    kz1 = zeros(nz_DAS_VSP,nx_DAS_VSP);
    
    kx2 = zeros(nz_DAS_VSP,nx_DAS_VSP);
    kz2 = zeros(nz_DAS_VSP,nx_DAS_VSP);
    
    X_DIRECT_PART = zeros(nz_DAS_VSP,nx_DAS_VSP);
    Z_DIRECT_PART = zeros(nz_DAS_VSP,nx_DAS_VSP);
    X_REFLECTION_PART = zeros(nz_DAS_VSP,nx_DAS_VSP);
    Z_REFLECTION_PART = zeros(nz_DAS_VSP,nx_DAS_VSP);
    X_DIRECT_PART_NEW = zeros(nz_DAS_VSP,nx_DAS_VSP);
    Z_DIRECT_PART_NEW = zeros(nz_DAS_VSP,nx_DAS_VSP);
    X_REFLECTION_PART_NEW = zeros(nz_DAS_VSP,nx_DAS_VSP);
    Z_REFLECTION_PART_NEW = zeros(nz_DAS_VSP,nx_DAS_VSP);
    Z_REFLECTION_PART_REAL = zeros(nz_DAS_VSP,nx_DAS_VSP);
    
    nx_reflection_point_lim = floor(nx_migration_offset/2) + 1;
    RATIO_Z_X = zeros(nz_DAS_VSP,nx_DAS_VSP);
    DZ_H_UP_RTT_DAS_VSP = zeros(nz_DAS_VSP,nx_DAS_VSP);
    DZ_H_UP_RTT_DAS_VSP_REAL = zeros(nz_DAS_VSP,nx_DAS_VSP);
    
    %%%%%%%%%%%%% 以通常思维进行编程 上行波%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:nz_DAS_VSP
        for j = 1:nx_DAS_VSP
            %   i,j为反射点横纵坐标
            X_CRP(i,j) = (j - 1)*dx; 
            Z_CRP(i,j) = (i - 1)*dz_trace_DAS_VSP;
            X_SOURCE = X_CRP(1,1);
            Z_SOURCE = Z_CRP(1,1);
            %   不要混淆位置和位移！
            if i <= nz_DAS_VSP && j >= nx_reflection_point_lim && j < nx_migration_offset
                DX_DIRECT_PART(i,j) = X_CRP(i,j) - X_SOURCE;
                DZ_DIRECT_PART(i,j) = Z_CRP(i,j) - Z_SOURCE;
                RATIO_Z_X(i,j) = DZ_DIRECT_PART(i,j)./DX_DIRECT_PART(i,j);
                UP_DIRECT_PART(i,j) = sqrt(DX_DIRECT_PART(i,j).^2 + DZ_DIRECT_PART(i,j).^2)./(V_qPE_ave_DAS_VSP(i,j).*dt);
                DX_REFLECTION_PART(i,j) = x_migration_offset - X_CRP(i,j);
                DZ_REFLECTION_PART(i,j) = DX_REFLECTION_PART(i,j).*RATIO_Z_X(i,j);
                X_DAS_VSP(i,j) = X_CRP(i,j) + DX_REFLECTION_PART(i,j);
                Z_DAS_VSP(i,j) = Z_CRP(i,j) - DZ_REFLECTION_PART(i,j);
                kz_DAS_VSP(i,j) = round(Z_DAS_VSP(i,j)/dz_trace_DAS_VSP) + 1;
    
                DX1(i,j) = DX_DIRECT_PART(i,j) - DX_REFLECTION_PART(i,j);
                DX2(i,j) = 2*DX_REFLECTION_PART(i,j);
                DZ1(i,j) = Z_DAS_VSP(i,j);
                DZ2(i,j) = DZ_REFLECTION_PART(i,j);
    
                DTT_DAS_VSP(i,j) = sqrt(DZ_DIRECT_PART(i,j).^2 + DX_DIRECT_PART(i,j).^2)./(V_qPE_ave_DAS_VSP(i,j).*dt);
                UP_RTT_DAS_VSP(i,j) = sqrt((2*DZ_DIRECT_PART(i,j) - Z_DAS_VSP(i,j)).^2 + X_DAS_VSP(i,j).^2)./(V_qPE_ave_DAS_VSP(i,j).*dt);
                UPWARD_RTT_DAS_VSP(i,j) = UP_RTT_DAS_VSP(i,j) - DTT_DAS_VSP(i,j);
    
                t_UP_RTT_DAS_VSP(kz_DAS_VSP(i,j),round(UP_RTT_DAS_VSP(i,j))) = Reflection_DAS_VSP(i,j);
                t_UPWARD_RTT_DAS_VSP(kz_DAS_VSP(i,j),round(DTT_DAS_VSP(i,j) + UPWARD_RTT_DAS_VSP(i,j))) = Reflection_DAS_VSP(i,j);
            end    
        end
    end
    max_kz_DAS_VSP = max(max(kz_DAS_VSP));
    for iz_DAS_VSP = 1:max_kz_DAS_VSP
        CONV_UPWARD_RTT_DAS_VSP = conv(t_UPWARD_RTT_DAS_VSP(iz_DAS_VSP,:),R_wavelet_DAS_VSP);
        SYN_UPWARD_RTT_DAS_VSP(iz_DAS_VSP,1:length(CONV_UPWARD_RTT_DAS_VSP)) = CONV_UPWARD_RTT_DAS_VSP;
        LENGTH_SYN_UPWARD_RTT_DAS_VSP = size_ricker_wavelet_DAS_VSP/2:1:length(SYN_UPWARD_RTT_DAS_VSP(iz_DAS_VSP,:)) + size_ricker_wavelet_DAS_VSP/2 - 1;
        for m = nt_ricker_DAS_VSP/2:1:(length(LENGTH_SYN_UPWARD_RTT_DAS_VSP))
            SYN_UPWARD_RTT_DAS_VSP_n(iz_DAS_VSP,m + 1 - nt_ricker_DAS_VSP/2) = SYN_UPWARD_RTT_DAS_VSP(iz_DAS_VSP,m);
        end
    end
    
    figure(122)
    wigb_VSP(SYN_UPWARD_RTT_DAS_VSP_n,1);
    set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    xlabel('t(ms)');
    ylabel('trace');
    title('synthetic upward travelling DAS VSP','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    
    
    %%%%%%%%%%%%% 以通常思维进行编程 一次下行波%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nx_reflection_point_min = floor(nx_migration_offset/3) + 1;
    
    for i = 1:nz_DAS_VSP
        for j = 1:nx_DAS_VSP
            %   i,j为反射点横纵坐标
            X_CRP(i,j) = (j - 1)*dx; 
            Z_CRP(i,j) = (i - 1)*dz_trace_DAS_VSP;
    
            X_CRP_SURFACE(i,j) = 2*X_CRP(i,j);
    
            X_SOURCE = X_CRP(1,1);
            Z_SOURCE = Z_CRP(1,1);        
    
            if i <= nz_DAS_VSP && j >= nx_reflection_point_min && j < nx_reflection_point_lim
                DX_DIRECT_PART(i,j) = X_CRP(i,j) - X_SOURCE;
                DZ_DIRECT_PART(i,j) = Z_CRP(i,j) - Z_SOURCE;
                RATIO_Z_X(i,j) = DZ_DIRECT_PART(i,j)./DX_DIRECT_PART(i,j);
                DOWN_DIRECT_PART(i,j) = sqrt(DX_DIRECT_PART(i,j).^2 + DZ_DIRECT_PART(i,j).^2)./(V_qPE_ave_DAS_VSP(i,j).*dt);
    
                DX_UPWARD_PART(i,j) = DX_DIRECT_PART(i,j);
                DZ_UPWARD_PART(i,j) = DZ_DIRECT_PART(i,j);
    
                DX_DOWNWARD_PART(i,j) = x_migration_offset - (DX_DIRECT_PART(i,j) + DX_UPWARD_PART(i,j));
                DZ_DOWNWARD_PART(i,j) = DX_DOWNWARD_PART(i,j).*RATIO_Z_X(i,j);
    
                X_DAS_VSP(i,j) = x_migration_offset;
                Z_DAS_VSP(i,j) = DZ_DOWNWARD_PART(i,j);
                kz_DAS_VSP(i,j) = round(Z_DAS_VSP(i,j)/dz_trace_DAS_VSP) + 1;
    
                DOWN_RTT_DAS_VSP(i,j) = sqrt((2*DZ_DIRECT_PART(i,j) + DZ_DOWNWARD_PART(i,j)).^2 + X_DAS_VSP(i,j).^2)./(V_qPE_ave_DAS_VSP(i,j).*dt);
                t_DOWN_RTT_DAS_VSP(kz_DAS_VSP(i,j),round(DOWN_RTT_DAS_VSP(i,j))) = Reflection_DAS_VSP(i,j);
    
            end
        end
    end
    
    max_kz_DAS_VSP = max(max(kz_DAS_VSP));
    for iz_DAS_VSP = 1:max_kz_DAS_VSP
        CONV_DOWN_RTT_DAS_VSP = conv(t_DOWN_RTT_DAS_VSP(iz_DAS_VSP,:),R_wavelet_DAS_VSP);
        SYN_DOWN_RTT_DAS_VSP(iz_DAS_VSP,1:length(CONV_DOWN_RTT_DAS_VSP)) = CONV_DOWN_RTT_DAS_VSP;
        LENGTH_SYN_DOWN_RTT_DAS_VSP = size_ricker_wavelet_DAS_VSP/2:1:length(SYN_DOWN_RTT_DAS_VSP(iz_DAS_VSP,:)) + size_ricker_wavelet_DAS_VSP/2 - 1;
        for m = nt_ricker_DAS_VSP/2:1:(length(LENGTH_SYN_DOWN_RTT_DAS_VSP))
            SYN_DOWN_RTT_DAS_VSP_n(iz_DAS_VSP,m + 1 - nt_ricker_DAS_VSP/2) = SYN_DOWN_RTT_DAS_VSP(iz_DAS_VSP,m);
        end
    end
    
    figure(133)
    wigb_VSP(SYN_DOWN_RTT_DAS_VSP_n,1);
    set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    xlabel('t(ms)');
    ylabel('trace');
    title('synthetic downward travelling DAS VSP','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    
    %%%%%%  整体VSP-DAS:直达 + 上行 + 一次下行  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:nz_DAS_VSP
        for j = 1:nx_DAS_VSP
            %   i,j为反射点横纵坐标
            X_CRP(i,j) = (j - 1)*dx; 
            Z_CRP(i,j) = (i - 1)*dz_trace_DAS_VSP;
            X_CRP_SURFACE(i,j) = 2*X_CRP(i,j);
    
            X_SOURCE = X_CRP(1,1);
            Z_SOURCE = Z_CRP(1,1);        
    
            if i <= nz_DAS_VSP && j >= nx_reflection_point_min && j < nx_reflection_point_lim
                DX_DIRECT_PART(i,j) = X_CRP(i,j) - X_SOURCE;
                DZ_DIRECT_PART(i,j) = Z_CRP(i,j) - Z_SOURCE;
                RATIO_Z_X(i,j) = DZ_DIRECT_PART(i,j)./DX_DIRECT_PART(i,j);
                DOWN_DIRECT_PART(i,j) = sqrt(DX_DIRECT_PART(i,j).^2 + DZ_DIRECT_PART(i,j).^2)./(V_qPE_ave_DAS_VSP(i,j).*dt);
    
                DX_UPWARD_PART(i,j) = DX_DIRECT_PART(i,j);
                DZ_UPWARD_PART(i,j) = DZ_DIRECT_PART(i,j);
    
                DX_DOWNWARD_PART(i,j) = x_migration_offset - (DX_DIRECT_PART(i,j) + DX_UPWARD_PART(i,j));
                DZ_DOWNWARD_PART(i,j) = DX_DOWNWARD_PART(i,j).*RATIO_Z_X(i,j);
    
                X_DAS_VSP(i,j) = x_migration_offset;
                Z_DAS_VSP(i,j) = DZ_DOWNWARD_PART(i,j);
                kz_DAS_VSP(i,j) = round(Z_DAS_VSP(i,j)/dz_trace_DAS_VSP) + 1;
    
                RTT_DAS_VSP(i,j) = sqrt((2*DZ_DIRECT_PART(i,j) + DZ_DOWNWARD_PART(i,j)).^2 + X_DAS_VSP(i,j).^2)./(V_qPE_ave_DAS_VSP(i,j).*dt);
                t_RTT_DAS_VSP(kz_DAS_VSP(i,j),round(RTT_DAS_VSP(i,j))) = Reflection_DAS_VSP(i,j);
    
                % t_DOWN_DAS_VSP_UP_RTT(kz_DAS_VSP(i,j),round(DOWN_RTT_DAS_VSP(i,j))) = t_DOWN_RTT_DAS_VSP(i,round(DOWN_RTT_DAS_VSP(i,j)));
            elseif i <= nz_DAS_VSP && j >= nx_reflection_point_lim && j < nx_migration_offset
                DX_DIRECT_PART(i,j) = X_CRP(i,j) - X_SOURCE;
                DZ_DIRECT_PART(i,j) = Z_CRP(i,j) - Z_SOURCE;
    
                RATIO_Z_X(i,j) = DZ_DIRECT_PART(i,j)./DX_DIRECT_PART(i,j);
    
                UP_DIRECT_PART(i,j) = sqrt(DX_DIRECT_PART(i,j).^2 + DZ_DIRECT_PART(i,j).^2)./(V_qPE_ave_DAS_VSP(i,j).*dt);
    
                DX_REFLECTION_PART(i,j) = x_migration_offset - X_CRP(i,j);
                DZ_REFLECTION_PART(i,j) = DX_REFLECTION_PART(i,j).*RATIO_Z_X(i,j);
    
                X_DAS_VSP(i,j) = X_CRP(i,j) + DX_REFLECTION_PART(i,j);
                Z_DAS_VSP(i,j) = Z_CRP(i,j) - DZ_REFLECTION_PART(i,j);
                kz_DAS_VSP(i,j) = round(Z_DAS_VSP(i,j)/dz_trace_DAS_VSP) + 1;
    
                DTT_DAS_VSP(i,j) = sqrt(DZ_DIRECT_PART(i,j).^2 + DX_DIRECT_PART(i,j).^2)./(V_qPE_ave_DAS_VSP(i,j).*dt);
                RTT_DAS_VSP(i,j) = sqrt((2*DZ_DIRECT_PART(i,j) - Z_DAS_VSP(i,j)).^2 + X_DAS_VSP(i,j).^2)./(V_qPE_ave_DAS_VSP(i,j).*dt);
    
                t_RTT_DAS_VSP(kz_DAS_VSP(i,j),round(RTT_DAS_VSP(i,j))) = Reflection_DAS_VSP(i,j);
            end 
        end
    end
    for iz_DAS_VSP = 1:max_kz_DAS_VSP
        CONV_RTT_DAS_VSP = conv(t_RTT_DAS_VSP(iz_DAS_VSP,:),R_wavelet_DAS_VSP);
        SYN_RTT_DAS_VSP(iz_DAS_VSP,1:length(CONV_RTT_DAS_VSP)) = CONV_RTT_DAS_VSP;
        LENGTH_SYN_RTT_DAS_VSP = size_ricker_wavelet_DAS_VSP/2:1:length(SYN_RTT_DAS_VSP(iz_DAS_VSP,:)) + size_ricker_wavelet_DAS_VSP/2 - 1;
        for m = nt_ricker_DAS_VSP/2:1:(length(LENGTH_SYN_RTT_DAS_VSP))
            SYN_RTT_DAS_VSP_n(iz_DAS_VSP,m + 1 - nt_ricker_DAS_VSP/2) = SYN_RTT_DAS_VSP(iz_DAS_VSP,m);
        end
    end
    
    figure(144)
    wigb_VSP(SYN_RTT_DAS_VSP_n,1);
    set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    xlabel('t(ms)');
    ylabel('trace');
    title('synthetic DAS VSP','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    
    figure(555)
    imagesc(x_main_region,z_main_region,V33_qPE_adi_I);
    xlabel('nx');
    ylabel('nz');
end

%%%%%%%%%%% 残差完全匹配层预定义数组 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    %   速度变量更新
    v_x_last_half = zeros(nz_PML,nx_PML);
    v_x_next_half = zeros(nz_PML,nx_PML);      
    v_z_last_half = zeros(nz_PML,nx_PML);
    v_z_next_half = zeros(nz_PML,nx_PML);

    v_x_stage1 = zeros(nz_PML,nx_PML);
    v_x_stage2 = zeros(nz_PML,nx_PML);
    v_x_stage3 = zeros(nz_PML,nx_PML);
    v_x_stage4 = zeros(nz_PML,nx_PML);

    v_z_stage1 = zeros(nz_PML,nx_PML);
    v_z_stage2 = zeros(nz_PML,nx_PML);
    v_z_stage3 = zeros(nz_PML,nx_PML);
    v_z_stage4 = zeros(nz_PML,nx_PML);

    v_z_z_last_half_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    v_z_z_next_half_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    v_z_x_last_half_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    v_z_x_next_half_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    v_x_z_last_half_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    v_x_z_next_half_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    v_x_x_last_half_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    v_x_x_next_half_MCFS_RPML_II = zeros(nz_PML,nx_PML);


    %   速度残差变量更新
    epsilon_v_x_x_1_last_half = zeros(nz_PML,nx_PML);
    epsilon_v_x_x_1_next_half = zeros(nz_PML,nx_PML);
    epsilon_v_x_x_2_last_half = zeros(nz_PML,nx_PML);
    epsilon_v_x_x_2_next_half = zeros(nz_PML,nx_PML);
    epsilon_v_x_z_1_last_half = zeros(nz_PML,nx_PML);
    epsilon_v_x_z_1_next_half = zeros(nz_PML,nx_PML);
    epsilon_v_x_z_2_last_half = zeros(nz_PML,nx_PML);
    epsilon_v_x_z_2_next_half = zeros(nz_PML,nx_PML);
    epsilon_v_z_x_1_last_half = zeros(nz_PML,nx_PML);
    epsilon_v_z_x_1_next_half = zeros(nz_PML,nx_PML);
    epsilon_v_z_x_2_last_half = zeros(nz_PML,nx_PML);
    epsilon_v_z_x_2_next_half = zeros(nz_PML,nx_PML);
    epsilon_v_z_z_1_last_half = zeros(nz_PML,nx_PML);
    epsilon_v_z_z_1_next_half = zeros(nz_PML,nx_PML);
    epsilon_v_z_z_2_last_half = zeros(nz_PML,nx_PML);
    epsilon_v_z_z_2_next_half = zeros(nz_PML,nx_PML);

    %   应力变量更新
    T_PK1_xx_present = zeros(nz_PML,nx_PML);
    T_PK1_xx_next = zeros(nz_PML,nx_PML);
    T_PK1_zz_present = zeros(nz_PML,nx_PML);
    T_PK1_zz_next = zeros(nz_PML,nx_PML);
    T_PK1_xz_present = zeros(nz_PML,nx_PML);
    T_PK1_xz_next = zeros(nz_PML,nx_PML);
    T_PK1_zx_present = zeros(nz_PML,nx_PML);
    T_PK1_zx_next = zeros(nz_PML,nx_PML);

    T_PK1_xx_stage1 = zeros(nz_PML,nx_PML);
    T_PK1_xx_stage2 = zeros(nz_PML,nx_PML); 
    T_PK1_xx_stage3 = zeros(nz_PML,nx_PML); 
    T_PK1_xx_stage4 = zeros(nz_PML,nx_PML);    
    T_PK1_zz_stage1 = zeros(nz_PML,nx_PML);
    T_PK1_zz_stage2 = zeros(nz_PML,nx_PML); 
    T_PK1_zz_stage3 = zeros(nz_PML,nx_PML); 
    T_PK1_zz_stage4 = zeros(nz_PML,nx_PML);    
    T_PK1_xz_stage1 = zeros(nz_PML,nx_PML);
    T_PK1_xz_stage2 = zeros(nz_PML,nx_PML); 
    T_PK1_xz_stage3 = zeros(nz_PML,nx_PML); 
    T_PK1_xz_stage4 = zeros(nz_PML,nx_PML);    
    T_PK1_zx_stage1 = zeros(nz_PML,nx_PML);
    T_PK1_zx_stage2 = zeros(nz_PML,nx_PML); 
    T_PK1_zx_stage3 = zeros(nz_PML,nx_PML); 
    T_PK1_zx_stage4 = zeros(nz_PML,nx_PML);    

    T_PK1_xx_x_present_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    T_PK1_xx_x_next_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    T_PK1_xx_z_present_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    T_PK1_xx_z_next_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    T_PK1_zz_x_present_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    T_PK1_zz_x_next_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    T_PK1_zz_z_present_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    T_PK1_zz_z_next_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    T_PK1_xz_x_present_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    T_PK1_xz_x_next_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    T_PK1_xz_z_present_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    T_PK1_xz_z_next_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    T_PK1_zx_x_present_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    T_PK1_zx_x_next_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    T_PK1_zx_z_present_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    T_PK1_zx_z_next_MCFS_RPML_II = zeros(nz_PML,nx_PML);


    %   应力残差变量更新
    epsilon_T_PK1_xx_x_1_present = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_xx_x_1_next = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_xx_x_2_present = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_xx_x_2_next = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_xx_z_1_present = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_xx_z_1_next = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_xx_z_2_present = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_xx_z_2_next = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_zz_x_1_present = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_zz_x_1_next = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_zz_x_2_present = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_zz_x_2_next = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_zz_z_1_present = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_zz_z_1_next = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_zz_z_2_present = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_zz_z_2_next = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_zx_x_1_present = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_zx_x_1_next = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_zx_x_2_present = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_zx_x_2_next = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_zx_z_1_present = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_zx_z_1_next = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_zx_z_2_present = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_zx_z_2_next = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_xz_x_1_present = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_xz_x_1_next = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_xz_x_2_present = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_xz_x_2_next = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_xz_z_1_present = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_xz_z_1_next = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_xz_z_2_present = zeros(nz_PML,nx_PML);
    epsilon_T_PK1_xz_z_2_next = zeros(nz_PML,nx_PML);

    %   热流变量更新
    q_x_last_half = zeros(nz_PML,nx_PML);
    q_x_next_half = zeros(nz_PML,nx_PML);
    q_z_last_half = zeros(nz_PML,nx_PML);
    q_z_next_half = zeros(nz_PML,nx_PML);

    q_x_mid_stage1 = zeros(nz_PML,nx_PML);
    q_x_mid_stage2 = zeros(nz_PML,nx_PML);
    q_x_mid_stage3 = zeros(nz_PML,nx_PML);
    q_x_mid_stage4 = zeros(nz_PML,nx_PML);

    q_z_mid_stage1 = zeros(nz_PML,nx_PML);
    q_z_mid_stage2 = zeros(nz_PML,nx_PML);
    q_z_mid_stage3 = zeros(nz_PML,nx_PML);
    q_z_mid_stage4 = zeros(nz_PML,nx_PML);

    q_z_z_last_half_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    q_z_z_next_half_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    q_z_x_last_half_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    q_z_x_next_half_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    q_x_z_last_half_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    q_x_z_next_half_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    q_x_x_last_half_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    q_x_x_next_half_MCFS_RPML_II = zeros(nz_PML,nx_PML);

    %   热流残差变量更新
    epsilon_q_x_x_1_last_half = zeros(nz_PML,nx_PML);
    epsilon_q_x_x_1_mid = zeros(nz_PML,nx_PML);
    epsilon_q_x_x_1_next_half = zeros(nz_PML,nx_PML);
    epsilon_q_x_x_2_last_half = zeros(nz_PML,nx_PML);
    epsilon_q_x_x_2_mid = zeros(nz_PML,nx_PML);
    epsilon_q_x_x_2_next_half = zeros(nz_PML,nx_PML);
    epsilon_q_x_z_1_last_half = zeros(nz_PML,nx_PML);
    epsilon_q_x_z_1_mid = zeros(nz_PML,nx_PML);
    epsilon_q_x_z_1_next_half = zeros(nz_PML,nx_PML);
    epsilon_q_x_z_2_last_half = zeros(nz_PML,nx_PML);
    epsilon_q_x_z_2_mid = zeros(nz_PML,nx_PML);
    epsilon_q_x_z_2_next_half = zeros(nz_PML,nx_PML);
    epsilon_q_z_x_1_last_half = zeros(nz_PML,nx_PML);
    epsilon_q_z_x_1_mid = zeros(nz_PML,nx_PML);
    epsilon_q_z_x_1_next_half = zeros(nz_PML,nx_PML);
    epsilon_q_z_x_2_last_half = zeros(nz_PML,nx_PML);
    epsilon_q_z_x_2_mid = zeros(nz_PML,nx_PML);
    epsilon_q_z_x_2_next_half = zeros(nz_PML,nx_PML);
    epsilon_q_z_z_1_last_half = zeros(nz_PML,nx_PML);
    epsilon_q_z_z_1_mid = zeros(nz_PML,nx_PML);
    epsilon_q_z_z_1_next_half = zeros(nz_PML,nx_PML);
    epsilon_q_z_z_2_last_half = zeros(nz_PML,nx_PML);
    epsilon_q_z_z_2_mid = zeros(nz_PML,nx_PML);
    epsilon_q_z_z_2_next_half = zeros(nz_PML,nx_PML);

    %   温度变量更新
    theta_present = zeros(nz_PML,nx_PML);
    theta_next = zeros(nz_PML,nx_PML);

    theta_stage1 = zeros(nz_PML,nx_PML);
    theta_stage2 = zeros(nz_PML,nx_PML);
    theta_stage3 = zeros(nz_PML,nx_PML);
    theta_stage4 = zeros(nz_PML,nx_PML);

    theta_x_present_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    theta_z_present_MCFS_RPML_II = zeros(nz_PML,nx_PML);

    theta_x_next_MCFS_RPML_II = zeros(nz_PML,nx_PML);
    theta_z_next_MCFS_RPML_II = zeros(nz_PML,nx_PML);

    %   温度残差变量更新
    epsilon_theta_x_1_present = zeros(nz_PML,nx_PML);
    epsilon_theta_x_1_next = zeros(nz_PML,nx_PML);
    epsilon_theta_x_2_present = zeros(nz_PML,nx_PML);
    epsilon_theta_x_2_next = zeros(nz_PML,nx_PML);
    epsilon_theta_z_1_present = zeros(nz_PML,nx_PML);
    epsilon_theta_z_1_next = zeros(nz_PML,nx_PML);
    epsilon_theta_z_2_present = zeros(nz_PML,nx_PML);
    epsilon_theta_z_2_next = zeros(nz_PML,nx_PML);

    %   自由能变量更新
    phi_total = zeros(nt_total,1);
    phi_interior = zeros(nt_total,1);
    phi_Helmholtz_present = zeros(nz_PML,nx_PML);
    phi_Helmholtz_next = zeros(nz_PML,nx_PML);

    %   空间导数
    dq_x_x_MCFS_RPML_II_dx = zeros(nz_PML,nx_PML);
    dq_x_z_MCFS_RPML_II_dz = zeros(nz_PML,nx_PML);
    dq_z_z_MCFS_RPML_II_dz = zeros(nz_PML,nx_PML);
    dq_z_x_MCFS_RPML_II_dx = zeros(nz_PML,nx_PML);