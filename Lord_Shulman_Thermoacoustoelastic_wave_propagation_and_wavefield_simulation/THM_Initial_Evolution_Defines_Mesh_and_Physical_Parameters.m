%   初始温度演化程序，无需PML边界条件

clear;
clc;

%   网格空间步长
dx = 20;        dz = 20;

dt = 5e-4;              %   dt一定要小于等于1e-4

%   非局部匹配层【重点】 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NLML_l = 20;

nx_main_region = (10000/dx) + 1;     %   X方向总点数
nz_main_region = (7000/dx) + 1;     %   Z方向总点数     

x = -NLML_l*dx:dx:(nx_main_region + NLML_l - 1)*dx;    % x方向上的离散化
nx_NLML = length(x);     %   x方向的网格数量
x_main_region = 0:dx:(nx_main_region - 1)*dx;

z = -NLML_l*dz:dz:(nz_main_region + NLML_l - 1)*dz;
nz_NLML = length(z);     %   z方向的网格数量
z_main_region = 0:dz:(nz_main_region - 1)*dz;

dfx = 1/(nx_NLML*dx);        dfz = 1/(nz_NLML*dz);
%%%%%%%% 数组预定义 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vp = zeros(nz_NLML,nx_NLML);
Vs = zeros(nz_NLML,nx_NLML);
rho = zeros(nz_NLML,nx_NLML);

lambda = zeros(nz_NLML,nx_NLML);
mu_0 = zeros(nz_NLML,nx_NLML);
Q_0 = zeros(nz_NLML,nx_NLML);

eta_viscoelastic = zeros(nz_NLML,nx_NLML);

%%%% 物性参数分布 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parisio et al.(2019); Farina et al.(2022)

%   断裂岩石设定为典型的结晶岩
K_host_dry_bulk = zeros(nz_NLML,nx_NLML);         %   实际上就是K
K_host_dry_bulk(1:nz_NLML,1:nx_NLML) = 40e9;

save K_host_dry_bulk_pre_invasion.mat K_host_dry_bulk

mu_host_dry = zeros(nz_NLML,nx_NLML);
mu_host_dry(1:nz_NLML,1:nx_NLML) = 24e9;
save mu_host_dry_pre_invasion.mat mu_host_dry

alpha_host_Biot = zeros(nz_NLML,nx_NLML);
alpha_host_Biot(1:nz_NLML,1:nx_NLML) = 0.5;
save alpha_host_Biot_pre_invasion.mat alpha_host_Biot

K_host_s_bulk = zeros(nz_NLML,nx_NLML);
K_host_s_bulk(1:nz_NLML,1:nx_NLML) = 80e9;

save K_host_s_bulk_pre_invasion.mat K_host_s_bulk

rho_host_s = zeros(nz_NLML,nx_NLML);
rho_host_s(1:nz_NLML,1:nx_NLML) = 2.7e3;

save rho_host_s_pre_invasion.mat rho_host_s

phi_host = zeros(nz_NLML,nx_NLML);
phi_host(1:nz_NLML,1:nx_NLML) = 0.01;

save phi_host_pre_invasion.mat phi_host

%%%%% 杨氏模量和体积模量 %%%%%%%%%%%%%%%%%%%%
E_youngs_host = zeros(nz_NLML,nx_NLML);
E_youngs_host(1:nz_NLML,1:nx_NLML) = 60e9;

save E_youngs_host_pre_invasion.mat E_youngs_host

nu_poisson_host = zeros(nz_NLML,nx_NLML);
nu_poisson_host(1:nz_NLML,1:nx_NLML) = 0.25;

save nu_poisson_host_pre_invasion.mat nu_poisson_host

K_host_dry_bulk = E_youngs_host./(3.*(1 - 2.*nu_poisson_host));
K_host_s_bulk = K_host_dry_bulk./(1 - alpha_host_Biot);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   拉梅常数
lambda_host = K_host_dry_bulk - (2/3).*mu_host_dry;

save lambda_host_pre_invasion.mat lambda_host
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rho_host_s_frame = (1 - phi_host).*rho_host_s;

%   断层岩石物理

%   储层热物性
cs_s = zeros(nz_NLML,nx_NLML);
cs_s(1:nz_NLML,1:nx_NLML) = 950;

alpha_s = zeros(nz_NLML,nx_NLML);
alpha_s(1:nz_NLML,1:nx_NLML) = 1e-5;

save alpha_s_pre_invasion.mat alpha_s


K_T_s = zeros(nz_NLML,nx_NLML);
K_T_s(1:nz_NLML,1:nx_NLML) = 3;

% %   储层Arrhenius参数
% n_stress = 3.7;
% E_Arrhenius = 59e3;
% A_Arrhenius = 1.3e-9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   初始化渗透率

K_p_intact = (4.979e-11).*phi_host.^3.11;         %   完好无损岩石渗透率    3.0001e-17，很小
K_p_fracture= (1.143e-11).*phi_host.^0.64;          %   完全破裂岩石渗透率

omega_fracture = 0;                                         %   破裂因子

%   在完好岩石渗透率和完全破裂岩石渗透率之间插值
K_p = exp((1 - omega_fracture).*log(K_p_intact) + omega_fracture.*log(K_p_fracture));       %   实际渗透率(Lamur et al.,2017)

%   在超临界地热储层勘探过程中，我们令omega_fracture = 0，因为不考虑注水开发。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   地温，假设地面温度为300K
%   地温梯度

%%%% 初始迭代，5次就差不多了 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nk = 5;

P_hydrostatic_plot = zeros(1,nk);
rho_w_plot = zeros(1,nk);


T_supercritical_bottom = 550;       %   超临界地热储层底边界，摄氏度
temperature_gradient = 30;         %   

ZEROS = zeros(nz_NLML,nx_NLML);

T = ZEROS;
p = ZEROS;

% rho_w = zeros(nz,nx);
rho_w_vapor = ZEROS;
rho_w_liquid = ZEROS;
h_w_vapor = ZEROS;
h_w_liquid = ZEROS;
s_w_vapor = ZEROS;
s_w_liquid = ZEROS;
vapour_volume_fraction = ZEROS;
vapour_fraction = ZEROS;
rho_w = ZEROS;
rho_w_p_variation = ZEROS;
rho_w_T_variation = ZEROS;
specific_volume_w_vapor = ZEROS;
specific_volume_w_liquid = ZEROS;

eta_w = ZEROS;      cp_w = ZEROS;       cv_w = ZEROS;       K_T_w = ZEROS;

beta_w = ZEROS;
alpha_w = ZEROS;

%   垂向压力

sigma_top_atmosphere = 1.01325e5;       %   标准大气压
sigma_v = zeros(nz_NLML,nx_NLML);     sigma_H = zeros(nz_NLML,nx_NLML);     sigma_h = zeros(nz_NLML,nx_NLML);
sigma_v_s = zeros(nz_NLML,nx_NLML);     sigma_H_s = zeros(nz_NLML,nx_NLML);     sigma_h_s = zeros(nz_NLML,nx_NLML);   

bar_sigma_v_s = ZEROS;
P_hydrostatic = ZEROS;
bar_P_hydrostatic = ZEROS;

rho_w(1:nz_NLML,1:nx_NLML) = 1e3;         %   通常情况下的水密度



% eta_w(1:nz,1:nx,1) = 8.9e-4;      %   水的黏度
% rho = rho_w.*phi_host + (1 - phi_host).*rho_host_s;

g_z = 9.8;     %  定义成向量分量！！！      
g_x = 0;
%%%%%%% 所以引出一个问题：地球上的笛卡尔坐标系（二维情况）怎么选取？%%%%%%%

%   当考虑到地球尺度，存在垂向应力和水平向应力，同时重力加速度向量的时候，最好将笛卡尔坐标系z轴取为引力加速度G（指向地心）方向！
%   Gz分量，即g，为重力加速度
%   Gx分量，为切向（地表）加速度

h_waitbar = waitbar(0,'please wait ... ');

% Mode = "elastic";

Mode = "thermoacoustoelastic";


for k = 1:nk

    str = ['正在初始化温度、压力和超临界水性质......',num2str((k/nk).*100),'%'];
    waitbar(k/nk,h_waitbar,str)    
    for i = 1:nz_NLML
        for j = 1:nx_NLML
    
            if i >= 1 && i <= NLML_l        %   吸收边界位置
                T(i,j) = 300;
                if Mode == "elastic" || Mode == "thermoelastic"
                    P_hydrostatic(i,j) = sigma_top_atmosphere;
                    bar_P_hydrostatic(i,j) = (1e-5).*P_hydrostatic(i,j);
                elseif Mode == "acoustoelastic" || Mode == "thermoacoustoelastic"
                    P_hydrostatic(i,j) = sigma_top_atmosphere;
                    bar_P_hydrostatic(i,j) = (1e-5).*P_hydrostatic(i,j);
                end
            end

            %   主计算域
            if Mode == "elastic" || Mode == "thermoelastic"
                if i >= NLML_l + 1 && i <= nz_NLML - NLML_l
                    if Mode == "elastic"
                        T(i,j) = 300;
                    elseif Mode == "thermoelastic"
                        T(i,j) = 300 + (((i - NLML_l - 1).*dz)./1000).*temperature_gradient;
                    end

                    P_hydrostatic(i,j) = sigma_top_atmosphere;
                    bar_P_hydrostatic(i,j) = (1e-5).*P_hydrostatic(i,j);
    
                    if j >= (4000/dx) + NLML_l + 1 && j <= (6000/dx) + NLML_l
                        T(nz_NLML - NLML_l - (2000/dz):nz_NLML - NLML_l - (1000/dz),j) = 550 + 273.15;
                    end
                end         
            elseif Mode == "acoustoelastic" || Mode == "thermoacoustoelastic"
                if i >= NLML_l + 1 && i <= nz_NLML - NLML_l
                    if Mode == "acoustoelastic"
                        T(i,j) = 300;
                    elseif Mode == "thermoacoustoelastic"
                        T(i,j) = 300 + (((i - NLML_l - 1).*dz)./1000).*temperature_gradient;
                    end

                    P_hydrostatic(i,j) = sigma_top_atmosphere + rho_w(i,j).*g_z.*z(i);
                    bar_P_hydrostatic(i,j) = (1e-5).*P_hydrostatic(i,j);
    
                    if j >= (4000/dx) + NLML_l + 1 && j <= (6000/dx) + NLML_l
                        T(nz_NLML - NLML_l - (2000/dz):nz_NLML - NLML_l - (1000/dz),j) = 550 + 273.15;
                    end
                end
            end
            if Mode == "elastic" || Mode == "thermoelastic"
                if i >= nz_NLML - NLML_l + 1 && i <= nz_NLML
                    % T(nz_NLML - NLML_l + 1:nz_NLML,j) = 300 + (((nz_NLML - 2*NLML_l - 1).*dz)./1000).*temperature_gradient;
                    T(i,j) = 300 + (((i - NLML_l - 1).*dz)./1000).*temperature_gradient;
                    P_hydrostatic(i,j) = sigma_top_atmosphere;
                    bar_P_hydrostatic(i,j) = (1e-5).*P_hydrostatic(i,j);
                end
            elseif Mode == "acoustoelastic" || Mode == "thermoacoustoelastic"
                if i >= nz_NLML - NLML_l + 1 && i <= nz_NLML
                    % T(nz_NLML - NLML_l + 1:nz_NLML,j) = 300 + (((nz_NLML - 2*NLML_l - 1).*dz)./1000).*temperature_gradient;
                    T(i,j) = 300 + (((i - NLML_l - 1).*dz)./1000).*temperature_gradient;
                    P_hydrostatic(i,j) = sigma_top_atmosphere + rho_w(i,j).*g_z.*z(i);
                    bar_P_hydrostatic(i,j) = (1e-5).*P_hydrostatic(i,j);
                end
            end

%             p(nz,j) = XSteam('psat_T',T(nz,j) - 273.15);    
            rho_w(i,j)                    = XSteam('rho_pT',bar_P_hydrostatic(i,j),T(i,j) - 273.15);      %       压力(bar)转换到MPa，温度转换为摄氏度！
            rho_w_p_variation(i,j)        = XSteam('rho_pT',bar_P_hydrostatic(i,j) + (1e-5).*bar_P_hydrostatic(i,j),T(i,j) - 273.15);
            rho_w_T_variation(i,j)        = XSteam('rho_pT',bar_P_hydrostatic(i,j),T(i,j) - 273.15 + (1e-5).*(T(i,j) - 273.15));

            eta_w(i,j) = XSteam('my_pT',bar_P_hydrostatic(i,j),T(i,j) - 273.15);        %       计算水的黏度
            cp_w(i,j) = XSteam('Cp_pT',bar_P_hydrostatic(i,j),T(i,j) - 273.15);        %       计算水的等压比热容 注意，单位为kJ
            cv_w(i,j) = XSteam('Cv_pT',bar_P_hydrostatic(i,j),T(i,j) - 273.15);        %       计算水的等容比热容
            K_T_w(i,j) = XSteam('tc_pT',bar_P_hydrostatic(i,j),T(i,j) - 273.15);        %       计算水的热导率

            %%%%% 有限差分法计算可压缩性和热导率 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            beta_w(i,j)   =  (1./rho_w(i,j)).*((rho_w_p_variation(i,j) - rho_w(i,j))./((1e-5).*bar_P_hydrostatic(i,j)));
            alpha_w(i,j) = -(1./rho_w(i,j)).*((rho_w_T_variation(i,j) - rho_w(i,j))./((1e-5).*(T(i,j) - 273.15)));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %   流体体积模量
            K_w(i,j) = 1./beta_w(i,j);

            % if Mode == "elastic"
            %     rho_w(i,j) = 1e3;
            % end

        end
    end
%     P_hydrostatic_plot(k) = P_hydrostatic(701,1).*(1e-6);
%     rho_w_plot(k) = rho_w(701,1);
end

%%%%%% Biot模量 （理想孔隙介质假设）%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = K_w.*K_host_s_bulk.^2./(K_w.*(K_host_s_bulk - K_host_dry_bulk) + phi_host.*K_host_s_bulk.*(K_host_s_bulk - K_w));

%%%%%% 不排水体积模量（考虑流体影响的体积模量）%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_host_undrained = K_host_dry_bulk + K_w.*(K_host_s_bulk - K_host_dry_bulk).^2./(K_w.*(K_host_s_bulk - K_host_dry_bulk) + phi_host.*K_host_s_bulk.*(K_host_s_bulk - K_w));





delete(h_waitbar);

save T_pre_invasion.mat T 
save P_pre_invasion.mat P_hydrostatic
save rho_w_pre_invasion.mat rho_w
save eta_w_pre_invasion.mat eta_w
save cv_w_pre_invasion.mat cv_w
save cp_w_pre_invasion.mat cp_w
save K_T_w_pre_invasion.mat K_T_w
save beta_w_pre_invasion.mat beta_w
save alpha_w_pre_invasion.mat alpha_w

save K_w_pre_invasion.mat K_w

save K_host_undrained.mat K_host_undrained




figure(1)
imagesc(x_main_region,z_main_region,T(NLML_l + 1:nz_NLML - NLML_l,NLML_l + 1:nx_NLML - NLML_l));
axis equal;
% contourf(x,z,T,[300 823.15],'ShowText','on');
set(gca,'YDir','reverse');
set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
title('Temperature(T)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
colormap(jet);
shading interp;
xlabel('x(m)');
ylabel('z(m)');
xlim([0 10000]);
ylim([0 7000]);

annotation('textbox',[0.15,0.4,0.25,0.5],'LineStyle','none','LineWidth',2,'String',strings,'Fontname','Times New Roman','fontweight','bold');

figure(2)
imagesc(x_main_region,z_main_region,P_hydrostatic(NLML_l + 1:nz_NLML - NLML_l,NLML_l + 1:nx_NLML - NLML_l));
axis equal;
% contourf(x,z,P_hydrostatic,[101325 5.9394e7],'ShowText','on');
set(gca,'YDir','reverse');
set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
title('Pressure(MPa)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
colormap(cool);
shading interp;
xlabel('x(m)');
ylabel('z(m)');
xlim([0 10000]);
ylim([0 7000]);

annotation('textbox',[0.15,0.4,0.25,0.5],'LineStyle','none','LineWidth',2,'String',strings,'Fontname','Times New Roman','fontweight','bold');




%%%%%%%%%%%%%%%%% 检查迭代次数效果 %%%%%%%%%%%%%%%%%%%%%%%
% ik = 1:nk;

% figure(1)
% plot(ik,P_hydrostatic_plot,'b');
% set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
% title('P','Fontname','Times New Roman','fontweight','bold','FontSize',10);
% xlabel('iteration');
% ylabel('P(MPa)');
% 
% figure(2)
% plot(ik,rho_w_plot,'b');
% set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
% title('\rho_w','Fontname','Times New Roman','fontweight','bold','FontSize',10);
% xlabel('iteration');
% ylabel('\rho_w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   设置模型之时，除了人为设置的超临界侵入体超过了超临界温度374摄氏度(647.15K)之外，其余地方都没有超过超临界温度
  
%%%%% 计算复合性质 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho = (1 - phi_host).*rho_host_s + phi_host.*rho_w;
c_s = ((1 - phi_host).*rho_host_s.*cs_s  + phi_host.*rho_w.*((1e3).*cv_w))./rho;
K_T = (1 - phi_host).*K_T_s  + phi_host.*K_T_w;

save rho_pre_invasion.mat rho
save c_s_pre_invasion.mat c_s
save K_T_pre_invasion.mat K_T

for i = 1:nz_NLML
    for j = 1:nx_NLML
        sigma_v(i,j) = sigma_top_atmosphere + rho(i,j).*g_z.*z(i);   %   岩石力学——垂向应力
        sigma_H(i,j) = 0.65.*sigma_v(i,j);    %   岩石力学——水平应力
    end
end

save sigma_v_pre_invasion.mat sigma_v
save sigma_H_pre_invasion.mat sigma_H

figure(3)
imagesc(x_main_region,z_main_region,sigma_v(NLML_l + 1:nz_NLML - NLML_l,NLML_l + 1:nx_NLML - NLML_l));
axis equal;
% contourf(x,z,P_hydrostatic,[101325 5.9394e7],'ShowText','on');
set(gca,'YDir','reverse');
set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
title('\sigma_v(MPa)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
colormap(cool);
shading interp;
xlabel('x(m)');
ylabel('z(m)');
xlim([0 10000]);
ylim([0 7000]);

annotation('textbox',[0.15,0.4,0.25,0.5],'LineStyle','none','LineWidth',2,'String',strings,'Fontname','Times New Roman','fontweight','bold');


v_permeability_x_next_half = ZEROS;
v_permeability_z_next_half = ZEROS;



%%%%%%%% 预定义场变量 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   岩石力学

sigma_zz_present = sigma_v;
sigma_yy_present = sigma_H;
sigma_xx_present = sigma_H;
sigma_zx_present = ZEROS;
sigma_xz_present = ZEROS;





%%%%%%%  岩石力学的垂向应力和水平应力，需要转换为弹性力学/地震波动力学中的正应力和剪应力 %%%
%   注意，重力(gravity)加速度向量是"竖直向下的"，而非"指向地心的"！
%   只有在极点和赤道位置上，引力(gravitation)和重力方向才是重合的！
%   Parisio(2019)的文章中定义的是重力加速度，但这还是有问题的。
%   在本文的情境中可以将二者视为等价。

%   尽管有横向应力，但是横向应力并不随空间而变化（至少在侵入前），近似为均匀分布。
%   重力加速度向量指向z方向。




