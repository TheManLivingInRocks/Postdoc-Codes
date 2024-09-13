clear;
clc;
warning off;
tic;

slCharacterEncoding = 'UTF-8';
%%%%%%% 网格参数和物性参数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Defines_Mesh_and_Physical_Parameters_Elastic();

pole = 1;
%%%%%%%% MCFS-NPML参数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[eta_x,eta_z,beta_x,beta_z,...
 alpha_x_PML,alpha_z_PML,...
 alpha_x_PML_optimized,alpha_z_PML_optimized] = MCFS_NPML_Boundary(dx,dz,dt,nz_PML,nx_PML,PML_l,PML_x,PML_z,...
                                                                            f0,pole,eta_0,beta_0_x,beta_0_z,...
                                                                            reflection_coefficient,n_power,V55_iso_I, ...
                                                                            P_eta,P_beta,gamma_damping,delta_damping,P_stable,attenuation_type,media);

%%%%%%%% 合成地震记录数组 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seismogram_1_v_x = zeros(nt,nx_PML);
seismogram_1_v_z = zeros(nt,nx_PML);

%%%%%%%% 合成温度记录 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thermogram_1_theta = zeros(nt,nx_PML);

seismogram_1_v_x_DAS_VSP = zeros(nt,nz_PML);
seismogram_1_v_z_DAS_VSP = zeros(nt,nz_PML);
thermogram_1_theta_DTS_VSP = zeros(nt,nz_PML);

%%%%%%%% 时间记录 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_record_v_z_fz_1500ms = zeros(nt,3);
time_record_v_x_fz_1500ms = zeros(nt,3);
time_record_q_z_fz_1500ms = zeros(nt,3);
time_record_q_x_fz_1500ms = zeros(nt,3);
time_record_theta_fz_1500ms = zeros(nt,3);

%   能量数组
phi_total = zeros(nt,1);
phi_interior = zeros(nt,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   网格类型  %%%%%%%%%%%%%%%%%%%%%%%%%
grid_type = "RSG";

%%%%%%%% 震源类型 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source_type = "explosive";     %  爆炸震源
a_cerjan = 0.04;
Cerjan_Absorbing_Boundary_Condition = Cerjan_Absorbing_Boundary(a_cerjan,nx_PML,nz_PML,PML_l);


synthetic_seismogram_flag = "false";    %   合成地震记录
time_record_flag = "false";             %   时间记录

%%%%%%%% 时间步进流程 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:nt
    t    = dt*n;
    fprintf('time = %.10f\n',t);

    frac1 = n/(nt);
    progressbar(frac1); 
    if t >= t0
        if source_type == "vertical"
            for i = 1:nz_PML
                for j = 1:nx_PML
                    Ricker_wavelet(i,j) = (1 - 2*(pi*f0*(t - t0)).^2).*exp(-(pi*f0*(t - t0)).^2)*exp(-0.25*((j - nx0).^2 + (i - nz0).^2));
                    v_z_last_half(i,j) = v_z_next_half(i,j) - Ricker_wavelet(i,j)./Rho_I(i,j);        % - 2*t0
                    v_z_z_MCFS_NPML_last_half(i,j) = (beta_z(i,j)/dt - (eta_z(i,j)*beta_z(i,j) + alpha_z_PML_optimized(i,j))/2).^(-1)*...
                                                     ((beta_z(i,j)/dt + (eta_z(i,j)*beta_z(i,j) + alpha_z_PML_optimized(i,j))/2)*v_z_z_MCFS_NPML_next_half(i,j) - ...
                                                      (eta_z(i,j)/2 + 1/dt)*v_z_next_half(i,j) - (eta_z(i,j)/2 - 1/dt)*v_z_last_half(i,j));
                    v_z_x_MCFS_NPML_last_half(i,j) = (beta_x(i,j)/dt - (eta_x(i,j)*beta_x(i,j) + alpha_x_PML_optimized(i,j))/2)^(-1)*...
                                                     ((beta_x(i,j)/dt + (eta_x(i,j)*beta_x(i,j) + alpha_x_PML_optimized(i,j))/2)*v_z_x_MCFS_NPML_next_half(i,j) - ...
                                                      (eta_x(i,j)/2 + 1/dt)*v_z_next_half(i,j) - (eta_x(i,j)/2 - 1/dt)*v_z_last_half(i,j));      
                end
            end        
        elseif source_type == "explosive"
            for i = 1:nz_PML
                for j = 1:nx_PML
                    Ricker_wavelet(i,j) = (1 - 2*(pi*f0*(t - t0)).^2).*exp(-(pi*f0*(t - t0)).^2)*exp(-0.25*((j - nx0).^2 + (i - nz0).^2));
                    T_PK1_xx_present(i,j) = T_PK1_xx_next(i,j) + Ricker_wavelet(i,j);
                    T_PK1_zz_present(i,j) = T_PK1_zz_next(i,j) + Ricker_wavelet(i,j);
                    T_PK1_zz_z_MCFS_NPML_present(i,j) = (beta_z(i,j)/dt - (eta_z(i,j)*beta_z(i,j) + alpha_z_PML_optimized(i,j))/2).^(-1)*...
                                                                          ((beta_z(i,j)/dt + (eta_z(i,j)*beta_z(i,j) + alpha_z_PML_optimized(i,j))/2)*T_PK1_zz_stage1_z_MCFS_NPML_next(i,j) - ...
                                                                           (eta_z(i,j)/2 + 1/dt)*T_PK1_zz_stage1_next(i,j) - (eta_z(i,j)/2 - 1/dt)*T_PK1_zz_present(i,j));
                    T_PK1_xx_x_MCFS_NPML_present(i,j) = (beta_x(i,j)/dt - (eta_x(i,j)*beta_x(i,j) + alpha_x_PML_optimized(i,j))/2).^(-1)*...
                                                                          ((beta_x(i,j)/dt + (eta_x(i,j)*beta_x(i,j) + alpha_x_PML_optimized(i,j))/2)*T_PK1_xx_stage1_x_MCFS_NPML_next(i,j) - ...
                                                                           (eta_x(i,j)/2 + 1/dt)*T_PK1_xx_stage1_next(i,j) - (eta_x(i,j)/2 - 1/dt)*T_PK1_xx_present(i,j)); 
                end
            end
        end
    end

    [dT_PK1_xx_x_MCFS_NPML_dx,dT_PK1_zz_z_MCFS_NPML_dz,...
     dT_PK1_xz_z_MCFS_NPML_dz,dT_PK1_zx_x_MCFS_NPML_dx] = Spatial_Derivative_RSG_PSO("second-order tensor",T_PK1_xx_x_MCFS_NPML_present,T_PK1_zz_z_MCFS_NPML_present,...
                                                                                        T_PK1_xz_z_MCFS_NPML_present,T_PK1_zx_x_MCFS_NPML_present,dx,dz,dfx,dfz);

    v_x_stage1_next_half = (v_x_last_half + 0.391752226571890*dt.*(1./Rho_I).*(dT_PK1_xx_x_MCFS_NPML_dx + dT_PK1_xz_z_MCFS_NPML_dz)).*Cerjan_Absorbing_Boundary_Condition;
    v_z_stage1_next_half = (v_z_last_half + 0.391752226571890*dt.*(1./Rho_I).*(dT_PK1_zx_x_MCFS_NPML_dx + dT_PK1_zz_z_MCFS_NPML_dz)).*Cerjan_Absorbing_Boundary_Condition;

        [v_x_stage1_x_MCFS_NPML_next_half,v_x_stage1_z_MCFS_NPML_next_half,...
         v_z_stage1_x_MCFS_NPML_next_half,v_z_stage1_z_MCFS_NPML_next_half] = MCFS_NPML_Iteration("vector",v_x_stage1_last_half, v_x_stage1_next_half,...
                                                                                                           v_x_stage1_x_MCFS_NPML_last_half, v_x_stage1_z_MCFS_NPML_last_half,...
                                                                                                           v_z_stage1_last_half, v_z_stage1_next_half,...
                                                                                                           v_z_stage1_x_MCFS_NPML_last_half,v_z_stage1_z_MCFS_NPML_last_half,...
                                                                                                           0,0,0,0,0,0,...
                                                                                                           dt,beta_x,beta_z,eta_x,eta_z,alpha_x_PML_optimized,alpha_z_PML_optimized);

        [dv_x_stage1_x_MCFS_NPML_dx,dv_z_stage1_z_MCFS_NPML_dz,...
            dv_x_stage1_z_MCFS_NPML_dz,dv_z_stage1_x_MCFS_NPML_dx] = Spatial_Derivative_RSG_PSO("vector",v_x_stage1_x_MCFS_NPML_next_half,v_z_stage1_z_MCFS_NPML_next_half,...
                                                                                                         v_x_stage1_z_MCFS_NPML_next_half,v_z_stage1_x_MCFS_NPML_next_half,dx,dz,dfx,dfz);

        theta_stage1_next = (theta_present + 0.391752226571890.*dt.*(- (1./cV_I).*(dq_x_x_MCFS_NPML_dx + dq_z_z_MCFS_NPML_dz) ...
                                                                     + (T_I./cV_I).*(M11_eqv_I.*dv_x_stage1_x_MCFS_NPML_dx + M33_eqv_I.*dv_z_stage1_z_MCFS_NPML_dz))).*Cerjan_Absorbing_Boundary_Condition;

        [theta_stage1_x_MCFS_NPML_next,theta_stage1_z_MCFS_NPML_next,OUTPUT3,OUTPUT4] = MCFS_NPML_Iteration("scalar",theta_stage1_present, theta_stage1_next,...
                                                                                                                     theta_stage1_x_MCFS_NPML_present, theta_stage1_z_MCFS_NPML_present,...
                                                                                                                     0,0,0,0,0,0,0,0,0,0,...
                                                                                                                     dt,beta_x,beta_z,eta_x,eta_z,alpha_x_PML_optimized,alpha_z_PML_optimized);

        [dtheta_stage1_x_MCFS_NPML_dx,dtheta_stage1_z_MCFS_NPML_dz,...
         dtheta_stage1_x_MCFS_NPML_dz,dtheta_stage1_z_MCFS_NPML_dx] = Spatial_Derivative_RSG_PSO("scalar",theta_stage1_x_MCFS_NPML_next,theta_stage1_z_MCFS_NPML_next,0,0,dx,dz,dfx,dfz);

        %%%%%  刚性热流 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        q_x_stage1_next_half = (q_x_mid + 0.391752226571890*dt.*(-1./tau11_I).*(K11_I.*dtheta_stage1_x_MCFS_NPML_dx)).*Cerjan_Absorbing_Boundary_Condition;
        q_z_stage1_next_half = (q_z_mid + 0.391752226571890*dt.*(-1./tau33_I).*(K33_I.*dtheta_stage1_z_MCFS_NPML_dz)).*Cerjan_Absorbing_Boundary_Condition;

        [q_x_stage1_x_MCFS_NPML_next_half,q_x_stage1_z_MCFS_NPML_next_half,...
         q_z_stage1_x_MCFS_NPML_next_half,q_z_stage1_z_MCFS_NPML_next_half] = MCFS_NPML_Iteration("vector",q_x_stage1_last_half, q_x_stage1_next_half,...
                                                                                                  q_x_stage1_x_MCFS_NPML_last_half, q_x_stage1_z_MCFS_NPML_last_half,...
                                                                                                  q_z_stage1_last_half, q_z_stage1_next_half,...
                                                                                                  q_z_stage1_x_MCFS_NPML_last_half, q_z_stage1_z_MCFS_NPML_last_half,...
                                                                                                  0,0,0,0,0,0,...
                                                                                                  dt,beta_x,beta_z,eta_x,eta_z,alpha_x_PML_optimized,alpha_z_PML_optimized);

        [dq_x_stage1_x_MCFS_NPML_dx,dq_z_stage1_z_MCFS_NPML_dz,...
         dq_x_stage1_z_MCFS_NPML_dz,dq_z_stage1_x_MCFS_NPML_dx] = Spatial_Derivative_RSG_PSO("vector",q_x_stage1_x_MCFS_NPML_next_half,q_z_stage1_z_MCFS_NPML_next_half,...
                                                                                                      q_x_stage1_z_MCFS_NPML_next_half,q_z_stage1_x_MCFS_NPML_next_half,dx,dz,dfx,dfz);

        T_PK1_xx_stage1_next = (T_PK1_xx_present + 0.391752226571890*(dt).*(A11_adi_I.*dv_x_stage1_x_MCFS_NPML_dx + A13_adi_I.*dv_z_stage1_z_MCFS_NPML_dz + ...
                                                                            A18_iso_I.*dv_x_stage1_z_MCFS_NPML_dz + A15_iso_I.*dv_z_stage1_x_MCFS_NPML_dx + ...
                                                                           (-M11_eqv_I./(Rho_I.*cV_I)).*(dq_x_stage1_x_MCFS_NPML_dx + dq_z_stage1_z_MCFS_NPML_dz))).*Cerjan_Absorbing_Boundary_Condition;           
        T_PK1_zz_stage1_next = (T_PK1_zz_present + 0.391752226571890*(dt).*(A13_adi_I.*dv_x_stage1_x_MCFS_NPML_dx + A33_adi_I.*dv_z_stage1_z_MCFS_NPML_dz + ...
                                                                            A38_iso_I.*dv_x_stage1_z_MCFS_NPML_dz + A35_iso_I.*dv_z_stage1_x_MCFS_NPML_dx + ...
                                                                           (-M33_eqv_I./(Rho_I.*cV_I)).*(dq_x_stage1_x_MCFS_NPML_dx + dq_z_stage1_z_MCFS_NPML_dz))).*Cerjan_Absorbing_Boundary_Condition;     
        T_PK1_xz_stage1_next = (T_PK1_xz_present + 0.391752226571890*(dt).*(A18_iso_I.*dv_x_stage1_x_MCFS_NPML_dx + A38_iso_I.*dv_z_stage1_z_MCFS_NPML_dz + ...
                                                                            A88_iso_I.*dv_x_stage1_z_MCFS_NPML_dz + A58_iso_I.*dv_z_stage1_x_MCFS_NPML_dx)).*Cerjan_Absorbing_Boundary_Condition;
        T_PK1_zx_stage1_next = (T_PK1_zx_present + 0.391752226571890*(dt).*(A15_iso_I.*dv_x_stage1_x_MCFS_NPML_dx + A35_iso_I.*dv_z_stage1_z_MCFS_NPML_dz + ...
                                                                            A58_iso_I.*dv_x_stage1_z_MCFS_NPML_dz + A55_iso_I.*dv_z_stage1_x_MCFS_NPML_dx)).*Cerjan_Absorbing_Boundary_Condition;

        T_PK1_xx_stage1_x_MCFS_NPML_next = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*...
                                            ((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*T_PK1_xx_stage1_x_MCFS_NPML_present + ...
                                             (eta_x./2 + 1./dt).*T_PK1_xx_stage1_next + (eta_x./2 - 1./dt).*T_PK1_xx_stage1_present); 
        T_PK1_zz_stage1_z_MCFS_NPML_next = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*...
                                            ((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*T_PK1_zz_stage1_z_MCFS_NPML_present + ...
                                             (eta_z./2 + 1./dt).*T_PK1_zz_stage1_next + (eta_z./2 - 1./dt).*T_PK1_zz_stage1_present);  
        T_PK1_xz_stage1_z_MCFS_NPML_next = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*...
                                            ((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*T_PK1_xz_stage1_z_MCFS_NPML_present + ...
                                             (eta_z./2 + 1./dt).*T_PK1_xz_stage1_next + (eta_z./2 - 1./dt).*T_PK1_xz_stage1_present);  
        T_PK1_zx_stage1_x_MCFS_NPML_next = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*...
                                            ((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*T_PK1_zx_stage1_x_MCFS_NPML_present + ...
                                             (eta_x./2 + 1./dt).*T_PK1_zx_stage1_next + (eta_x./2 - 1./dt).*T_PK1_zx_stage1_present);  
%%%%%% 阶段2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            dT_PK1_xx_stage1_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",T_PK1_xx_stage1_x_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            dT_PK1_zz_stage1_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",T_PK1_zz_stage1_z_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            dT_PK1_xz_stage1_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",T_PK1_xz_stage1_z_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            dT_PK1_zx_stage1_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",T_PK1_zx_stage1_x_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);

            v_x_stage2_next_half =  0.444370493651235.*v_x_last_half + 0.555629506348765.*v_x_stage1_next_half ...
                                                                     + 0.368410593050371*(dt).*(1./Rho_I).*(dT_PK1_xx_stage1_x_MCFS_NPML_dx + dT_PK1_xz_stage1_z_MCFS_NPML_dz);

            v_z_stage2_next_half =  0.444370493651235.*v_z_last_half + 0.555629506348765.*v_z_stage1_next_half ...
                                                                     + 0.368410593050371*(dt).*(1./Rho_I).*(dT_PK1_zx_stage1_x_MCFS_NPML_dx + dT_PK1_zz_stage1_z_MCFS_NPML_dz);
        
            v_x_stage2_next_half = v_x_stage2_next_half.*Cerjan_Absorbing_Boundary_Condition;
            v_z_stage2_next_half = v_z_stage2_next_half.*Cerjan_Absorbing_Boundary_Condition;

            v_x_stage2_x_MCFS_NPML_next_half = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*v_x_stage2_x_MCFS_NPML_last_half + ...
                                                             (eta_x./2 + 1./dt).*v_x_stage2_next_half + (eta_x./2 - 1./dt).*v_x_stage2_last_half);
            v_x_stage2_z_MCFS_NPML_next_half = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*v_x_stage2_z_MCFS_NPML_last_half + ...
                                                             (eta_z./2 + 1./dt).*v_x_stage2_next_half + (eta_z./2 - 1./dt).*v_x_stage2_last_half);    
            v_z_stage2_x_MCFS_NPML_next_half = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*v_z_stage2_x_MCFS_NPML_last_half + ...
                                                             (eta_x./2 + 1./dt).*v_z_stage2_next_half + (eta_x./2 - 1./dt).*v_z_stage2_last_half);
            v_z_stage2_z_MCFS_NPML_next_half = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*v_z_stage2_z_MCFS_NPML_last_half + ...
                                                             (eta_z./2 + 1./dt).*v_z_stage2_next_half + (eta_z./2 - 1./dt).*v_z_stage2_last_half);      

            dv_x_stage2_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",v_x_stage2_x_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            dv_z_stage2_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",v_z_stage2_z_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            dv_x_stage2_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",v_x_stage2_z_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            dv_z_stage2_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",v_z_stage2_x_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);

            theta_stage2_next = 0.444370493651235.*theta_present + 0.555629506348765.*theta_stage1_next ...
                                                                 + 0.368410593050371*dt.*(- (1./cV_I).*(dq_x_stage1_x_MCFS_NPML_dx + dq_z_stage1_z_MCFS_NPML_dz) ...
                                                                 + (T_I./cV_I).*(M11_eqv_I.*dv_x_stage2_x_MCFS_NPML_dx + M33_eqv_I.*dv_z_stage2_z_MCFS_NPML_dz));

            theta_stage2_next = theta_stage2_next.*Cerjan_Absorbing_Boundary_Condition;

            theta_stage2_x_MCFS_NPML_next = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*theta_stage2_x_MCFS_NPML_present + ...
                                                                        (eta_x./2 + 1./dt).*theta_stage2_next + (eta_x./2 - 1./dt).*theta_stage2_present);            
            theta_stage2_z_MCFS_NPML_next = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*theta_stage2_z_MCFS_NPML_present + ...
                                                                        (eta_z./2 + 1./dt).*theta_stage2_next + (eta_z./2 - 1./dt).*theta_stage2_present);    

            dtheta_stage2_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",theta_stage2_x_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            dtheta_stage2_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",theta_stage2_z_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);       

            q_x_stage2_next_half = 0.444370493651235.*q_x_mid + 0.555629506348765.*q_x_stage1_next_half ...
                                                              + 0.368410593050371*dt.*(-1./tau11_I).*(K11_I.*dtheta_stage2_x_MCFS_NPML_dx);
            q_z_stage2_next_half = 0.444370493651235.*q_z_mid + 0.555629506348765.*q_z_stage1_next_half ...
                                                              + 0.368410593050371*dt.*(-1./tau33_I).*(K33_I.*dtheta_stage2_z_MCFS_NPML_dz);

            q_x_stage2_next_half = q_x_stage2_next_half.*Cerjan_Absorbing_Boundary_Condition;
            q_z_stage2_next_half = q_z_stage2_next_half.*Cerjan_Absorbing_Boundary_Condition;

            q_x_stage2_x_MCFS_NPML_next_half = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*q_x_stage2_x_MCFS_NPML_last_half + ...
                                                             (eta_x./2 + 1./dt).*q_x_stage2_next_half + (eta_x./2 - 1./dt).*q_x_stage2_last_half);
            q_x_stage2_z_MCFS_NPML_next_half = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*q_x_stage2_z_MCFS_NPML_last_half + ...
                                                             (eta_z./2 + 1./dt).*q_x_stage2_next_half + (eta_z./2 - 1./dt).*q_x_stage2_last_half);  
            q_z_stage2_x_MCFS_NPML_next_half = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*q_z_stage2_x_MCFS_NPML_last_half + ...
                                                             (eta_x./2 + 1./dt).*q_z_stage2_next_half + (eta_x./2 - 1./dt).*q_z_stage2_last_half);
            q_z_stage2_z_MCFS_NPML_next_half = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*q_z_stage2_z_MCFS_NPML_last_half + ...
                                                             (eta_z./2 + 1./dt).*q_z_stage2_next_half + (eta_z./2 - 1./dt).*q_z_stage2_last_half);                 

            dq_x_stage2_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",q_x_stage2_x_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            dq_z_stage2_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",q_z_stage2_z_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);

            T_PK1_xx_stage2_next = 0.444370493651235.*T_PK1_xx_present + 0.555629506348765.*T_PK1_xx_stage1_next ... 
                                                                       + 0.368410593050371*(dt).*(A11_adi_I.*dv_x_stage2_x_MCFS_NPML_dx + A13_adi_I.*dv_z_stage2_z_MCFS_NPML_dz + ...
                                                                                                  A18_iso_I.*dv_x_stage2_z_MCFS_NPML_dz + A15_iso_I.*dv_z_stage2_x_MCFS_NPML_dx + ...
                                                                                                 (-M11_eqv_I./(Rho_I.*cV_I)).*(dq_x_stage2_x_MCFS_NPML_dx + dq_z_stage2_z_MCFS_NPML_dz));                    

            T_PK1_zz_stage2_next = 0.444370493651235.*T_PK1_zz_present + 0.555629506348765.*T_PK1_zz_stage1_next ... 
                                                                       + 0.368410593050371*(dt).*(A13_adi_I.*dv_x_stage2_x_MCFS_NPML_dx + A33_adi_I.*dv_z_stage2_z_MCFS_NPML_dz + ...
                                                                                                  A38_iso_I.*dv_x_stage2_z_MCFS_NPML_dz + A35_iso_I.*dv_z_stage2_x_MCFS_NPML_dx + ...
                                                                                                 (-M33_eqv_I./(Rho_I.*cV_I)).*(dq_x_stage2_x_MCFS_NPML_dx + dq_z_stage2_z_MCFS_NPML_dz));

            T_PK1_xz_stage2_next = 0.444370493651235.*T_PK1_xz_present + 0.555629506348765.*T_PK1_xz_stage1_next ... 
                                                                       + 0.368410593050371*(dt).*(A18_iso_I.*dv_x_stage2_x_MCFS_NPML_dx + A38_iso_I.*dv_z_stage2_z_MCFS_NPML_dz + ...
                                                                                                  A88_iso_I.*dv_x_stage2_z_MCFS_NPML_dz + A58_iso_I.*dv_z_stage2_x_MCFS_NPML_dx);                 
            
            T_PK1_zx_stage2_next = 0.444370493651235.*T_PK1_zx_present + 0.555629506348765.*T_PK1_zx_stage1_next ... 
                                                                       + 0.368410593050371*(dt).*(A15_iso_I.*dv_x_stage2_x_MCFS_NPML_dx + A35_iso_I.*dv_z_stage2_z_MCFS_NPML_dz + ...
                                                                                                  A58_iso_I.*dv_x_stage2_z_MCFS_NPML_dz + A55_iso_I.*dv_z_stage2_x_MCFS_NPML_dx);                
            
            T_PK1_xx_stage2_next = T_PK1_xx_stage2_next.*Cerjan_Absorbing_Boundary_Condition;
            T_PK1_zz_stage2_next = T_PK1_zz_stage2_next.*Cerjan_Absorbing_Boundary_Condition;
            T_PK1_xz_stage2_next = T_PK1_xz_stage2_next.*Cerjan_Absorbing_Boundary_Condition;
            T_PK1_zx_stage2_next = T_PK1_zx_stage2_next.*Cerjan_Absorbing_Boundary_Condition;


            T_PK1_xx_stage2_x_MCFS_NPML_next = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*...
                                            ((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*T_PK1_xx_stage2_x_MCFS_NPML_present + ...
                                             (eta_x./2 + 1./dt).*T_PK1_xx_stage2_next + (eta_x./2 - 1./dt).*T_PK1_xx_stage2_present); 
            T_PK1_zz_stage2_z_MCFS_NPML_next = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*...
                                            ((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*T_PK1_zz_stage2_z_MCFS_NPML_present + ...
                                             (eta_z./2 + 1./dt).*T_PK1_zz_stage2_next + (eta_z./2 - 1./dt).*T_PK1_zz_stage2_present);  
            T_PK1_xz_stage2_z_MCFS_NPML_next = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*...
                                            ((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*T_PK1_xz_stage2_z_MCFS_NPML_present + ...
                                             (eta_z./2 + 1./dt).*T_PK1_xz_stage2_next + (eta_z./2 - 1./dt).*T_PK1_xz_stage2_present);  
            T_PK1_zx_stage2_x_MCFS_NPML_next = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*...
                                            ((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*T_PK1_zx_stage2_x_MCFS_NPML_present + ...
                                             (eta_x./2 + 1./dt).*T_PK1_zx_stage2_next + (eta_x./2 - 1./dt).*T_PK1_zx_stage2_present);  

%%%%%%%%%% 阶段3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

            dT_PK1_xx_stage2_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",T_PK1_xx_stage2_x_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            dT_PK1_zz_stage2_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",T_PK1_zz_stage2_z_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            dT_PK1_xz_stage2_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",T_PK1_xz_stage2_z_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            dT_PK1_zx_stage2_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",T_PK1_zx_stage2_x_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);

            v_x_stage3_next_half =  0.620101851488403.*v_x_last_half + 0.379898148511597.*v_x_stage2_next_half ...
                                                                     + 0.251891774271694*(dt).*(1./Rho_I).*(dT_PK1_xx_stage2_x_MCFS_NPML_dx + dT_PK1_xz_stage2_z_MCFS_NPML_dz);

            v_z_stage3_next_half =  0.620101851488403.*v_z_last_half + 0.379898148511597.*v_z_stage2_next_half ...
                                                                     + 0.251891774271694*(dt).*(1./Rho_I).*(dT_PK1_zx_stage2_x_MCFS_NPML_dx + dT_PK1_zz_stage2_z_MCFS_NPML_dz);


            v_x_stage3_next_half = v_x_stage3_next_half.*Cerjan_Absorbing_Boundary_Condition;
            v_z_stage3_next_half = v_z_stage3_next_half.*Cerjan_Absorbing_Boundary_Condition;

            v_x_stage3_x_MCFS_NPML_next_half = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*v_x_stage3_x_MCFS_NPML_last_half + ...
                                                             (eta_x./2 + 1./dt).*v_x_stage3_next_half + (eta_x./2 - 1./dt).*v_x_stage3_last_half);
            v_x_stage3_z_MCFS_NPML_next_half = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*v_x_stage3_z_MCFS_NPML_last_half + ...
                                                             (eta_z./2 + 1./dt).*v_x_stage3_next_half + (eta_z./2 - 1./dt).*v_x_stage3_last_half);    
            v_z_stage3_x_MCFS_NPML_next_half = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*v_z_stage3_x_MCFS_NPML_last_half + ...
                                                             (eta_x./2 + 1./dt).*v_z_stage3_next_half + (eta_x./2 - 1./dt).*v_z_stage3_last_half);
            v_z_stage3_z_MCFS_NPML_next_half = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*v_z_stage3_z_MCFS_NPML_last_half + ...
                                                             (eta_z./2 + 1./dt).*v_z_stage3_next_half + (eta_z./2 - 1./dt).*v_z_stage3_last_half);      

            dv_x_stage3_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",v_x_stage3_x_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            dv_z_stage3_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",v_z_stage3_z_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            dv_x_stage3_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",v_x_stage3_z_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            dv_z_stage3_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",v_z_stage3_x_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);

            theta_stage3_next = 0.620101851488403.*theta_present + 0.379898148511597.*theta_stage2_next ...
                                                                                                        + 0.251891774271694*dt.*(- (1./cV_I).*(dq_x_stage2_x_MCFS_NPML_dx + dq_z_stage2_z_MCFS_NPML_dz) ...
                                                                                                                                 + (T_I./cV_I).*(M11_eqv_I.*dv_x_stage3_x_MCFS_NPML_dx + M33_eqv_I.*dv_z_stage3_z_MCFS_NPML_dz));

            theta_stage3_next = theta_stage3_next.*Cerjan_Absorbing_Boundary_Condition;

            theta_stage3_x_MCFS_NPML_next = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*theta_stage3_x_MCFS_NPML_present + ...
                                                                        (eta_x./2 + 1./dt).*theta_stage3_next + (eta_x./2 - 1./dt).*theta_stage3_present);            
            theta_stage3_z_MCFS_NPML_next = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*theta_stage3_z_MCFS_NPML_present + ...
                                                                        (eta_z./2 + 1./dt).*theta_stage3_next + (eta_z./2 - 1./dt).*theta_stage3_present);    

            dtheta_stage3_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",theta_stage3_x_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            dtheta_stage3_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",theta_stage3_z_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);       

            q_x_stage3_next_half = 0.620101851488403.*q_x_mid + 0.379898148511597.*q_x_stage2_next_half ...
                                                              + 0.251891774271694*dt.*(-1./tau11_I).*(K11_I.*dtheta_stage3_x_MCFS_NPML_dx);
            q_z_stage3_next_half = 0.620101851488403.*q_z_mid + 0.379898148511597.*q_z_stage2_next_half ...
                                                              + 0.251891774271694*dt.*(-1./tau33_I).*(K33_I.*dtheta_stage3_z_MCFS_NPML_dz);

            q_x_stage3_next_half = q_x_stage3_next_half.*Cerjan_Absorbing_Boundary_Condition;
            q_z_stage3_next_half = q_z_stage3_next_half.*Cerjan_Absorbing_Boundary_Condition;

            q_x_stage3_x_MCFS_NPML_next_half = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*q_x_stage3_x_MCFS_NPML_last_half + ...
                                                             (eta_x./2 + 1./dt).*q_x_stage3_next_half + (eta_x./2 - 1./dt).*q_x_stage3_last_half);
            q_x_stage3_z_MCFS_NPML_next_half = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*q_x_stage3_z_MCFS_NPML_last_half + ...
                                                             (eta_z./2 + 1./dt).*q_x_stage3_next_half + (eta_z./2 - 1./dt).*q_x_stage3_last_half);  
            q_z_stage3_x_MCFS_NPML_next_half = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*q_z_stage3_x_MCFS_NPML_last_half + ...
                                                             (eta_x./2 + 1./dt).*q_z_stage3_next_half + (eta_x./2 - 1./dt).*q_z_stage3_last_half);
            q_z_stage3_z_MCFS_NPML_next_half = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*q_z_stage3_z_MCFS_NPML_last_half + ...
                                                             (eta_z./2 + 1./dt).*q_z_stage3_next_half + (eta_z./2 - 1./dt).*q_z_stage3_last_half);                 

            dq_x_stage3_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",q_x_stage3_x_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            dq_z_stage3_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",q_z_stage3_z_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);


            T_PK1_xx_stage3_next = 0.620101851488403.*T_PK1_xx_present + 0.379898148511597.*T_PK1_xx_stage2_next ... 
                                                                       + 0.251891774271694*(dt).*(A11_adi_I.*dv_x_stage3_x_MCFS_NPML_dx + A13_adi_I.*dv_z_stage3_z_MCFS_NPML_dz + ...
                                                                                                  A18_iso_I.*dv_x_stage3_z_MCFS_NPML_dz + A15_iso_I.*dv_z_stage3_x_MCFS_NPML_dx + ...
                                                                                                  (-M11_eqv_I./(Rho_I.*cV_I)).*(dq_x_stage3_x_MCFS_NPML_dx + dq_z_stage3_z_MCFS_NPML_dz));

            T_PK1_zz_stage3_next = 0.620101851488403.*T_PK1_zz_present + 0.379898148511597.*T_PK1_zz_stage2_next ... 
                                                                       + 0.251891774271694*(dt).*(A13_iso_I.*dv_x_stage3_x_MCFS_NPML_dx + A33_iso_I.*dv_z_stage3_z_MCFS_NPML_dz + ...
                                                                                                  A38_iso_I.*dv_x_stage3_z_MCFS_NPML_dz + A35_iso_I.*dv_z_stage3_x_MCFS_NPML_dx + ...
                                                                                                 (-M33_eqv_I./(Rho_I.*cV_I)).*(dq_x_stage3_x_MCFS_NPML_dx + dq_z_stage3_z_MCFS_NPML_dz));   

            T_PK1_xz_stage3_next = 0.620101851488403.*T_PK1_xz_present + 0.379898148511597.*T_PK1_xz_stage2_next ... 
                                                                       + 0.251891774271694*(dt).*(A18_iso_I.*dv_x_stage3_x_MCFS_NPML_dx + A38_iso_I.*dv_z_stage3_z_MCFS_NPML_dz + ...
                                                                                                  A88_iso_I.*dv_x_stage3_z_MCFS_NPML_dz + A58_iso_I.*dv_z_stage3_x_MCFS_NPML_dx);                 
            
            T_PK1_zx_stage3_next = 0.620101851488403.*T_PK1_zx_present + 0.379898148511597.*T_PK1_zx_stage2_next ... 
                                                                       + 0.251891774271694*(dt).*(A15_iso_I.*dv_x_stage3_x_MCFS_NPML_dx + A35_iso_I.*dv_z_stage3_z_MCFS_NPML_dz + ...
                                                                                                  A58_iso_I.*dv_x_stage3_z_MCFS_NPML_dz + A55_iso_I.*dv_z_stage3_x_MCFS_NPML_dx);                
         
            T_PK1_xx_stage3_next = T_PK1_xx_stage3_next.*Cerjan_Absorbing_Boundary_Condition;
            T_PK1_zz_stage3_next = T_PK1_zz_stage3_next.*Cerjan_Absorbing_Boundary_Condition;
            T_PK1_xz_stage3_next = T_PK1_xz_stage3_next.*Cerjan_Absorbing_Boundary_Condition;
            T_PK1_zx_stage3_next = T_PK1_zx_stage3_next.*Cerjan_Absorbing_Boundary_Condition;

            T_PK1_xx_stage3_x_MCFS_NPML_next = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*...
                                                                                ((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*T_PK1_xx_stage3_x_MCFS_NPML_present + ...
                                                                                 (eta_x./2 + 1./dt).*T_PK1_xx_stage3_next + (eta_x./2 - 1./dt).*T_PK1_xx_stage3_present); 
            T_PK1_zz_stage3_z_MCFS_NPML_next = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*...
                                                                                ((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*T_PK1_zz_stage3_z_MCFS_NPML_present + ...
                                                                                 (eta_z./2 + 1./dt).*T_PK1_zz_stage3_next + (eta_z./2 - 1./dt).*T_PK1_zz_stage3_present);  
            T_PK1_xz_stage3_z_MCFS_NPML_next = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*...
                                                                                ((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*T_PK1_xz_stage3_z_MCFS_NPML_present + ...
                                                                                 (eta_z./2 + 1./dt).*T_PK1_xz_stage3_next + (eta_z./2 - 1./dt).*T_PK1_xz_stage3_present);  
            T_PK1_zx_stage3_x_MCFS_NPML_next = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*...
                                                                                ((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*T_PK1_zx_stage3_x_MCFS_NPML_present + ...
                                                                                 (eta_x./2 + 1./dt).*T_PK1_zx_stage3_next + (eta_x./2 - 1./dt).*T_PK1_zx_stage3_present);  


%%%%% 阶段4 %%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            dT_PK1_xx_stage3_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",T_PK1_xx_stage3_x_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            dT_PK1_zz_stage3_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",T_PK1_zz_stage3_z_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            dT_PK1_xz_stage3_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",T_PK1_xz_stage3_z_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            dT_PK1_zx_stage3_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",T_PK1_zx_stage3_x_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);

            v_x_stage4_next_half =  0.178079954393132.*v_x_last_half + 0.821920045606868.*v_x_stage3_next_half ...
                                                                                                           + 0.544974750228521*(dt).*(1./Rho_I).*(dT_PK1_xx_stage3_x_MCFS_NPML_dx + dT_PK1_xz_stage3_z_MCFS_NPML_dz);

            v_z_stage4_next_half =  0.178079954393132.*v_z_last_half + 0.821920045606868.*v_z_stage3_next_half ...
                                                                                                           + 0.544974750228521*(dt).*(1./Rho_I).*(dT_PK1_zx_stage3_x_MCFS_NPML_dx + dT_PK1_zz_stage3_z_MCFS_NPML_dz);
        
            v_x_stage4_next_half = v_x_stage4_next_half.*Cerjan_Absorbing_Boundary_Condition;
            v_z_stage4_next_half = v_z_stage4_next_half.*Cerjan_Absorbing_Boundary_Condition;

            v_x_stage4_x_MCFS_NPML_next_half = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*v_x_stage4_x_MCFS_NPML_last_half + ...
                                                             (eta_x./2 + 1./dt).*v_x_stage4_next_half + (eta_x./2 - 1./dt).*v_x_stage4_last_half);
            v_x_stage4_z_MCFS_NPML_next_half = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*v_x_stage4_z_MCFS_NPML_last_half + ...
                                                             (eta_z./2 + 1./dt).*v_x_stage4_next_half + (eta_z./2 - 1./dt).*v_x_stage4_last_half);    
            v_z_stage4_x_MCFS_NPML_next_half = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*v_z_stage4_x_MCFS_NPML_last_half + ...
                                                             (eta_x./2 + 1./dt).*v_z_stage4_next_half + (eta_x./2 - 1./dt).*v_z_stage4_last_half);
            v_z_stage4_z_MCFS_NPML_next_half = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*v_z_stage4_z_MCFS_NPML_last_half + ...
                                                             (eta_z./2 + 1./dt).*v_z_stage4_next_half + (eta_z./2 - 1./dt).*v_z_stage4_last_half);      

            dv_x_stage4_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",v_x_stage4_x_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            dv_z_stage4_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",v_z_stage4_z_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            dv_x_stage4_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",v_x_stage4_z_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            dv_z_stage4_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",v_z_stage4_x_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);

            theta_stage4_next = 0.178079954393132.*theta_present + 0.821920045606868.*theta_stage3_next ...
                                                                                                        + 0.544974750228521*dt.*(- (1./cV_I).*(dq_x_stage3_x_MCFS_NPML_dx + dq_z_stage3_z_MCFS_NPML_dz) ...
                                                                                                                                                    + (T_I./cV_I).*(M11_eqv_I.*dv_x_stage4_x_MCFS_NPML_dx + M33_eqv_I.*dv_z_stage4_z_MCFS_NPML_dz));

            theta_stage4_next = theta_stage4_next.*Cerjan_Absorbing_Boundary_Condition;

            theta_stage4_x_MCFS_NPML_next = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*theta_stage4_x_MCFS_NPML_present + ...
                                                                        (eta_x./2 + 1./dt).*theta_stage4_next + (eta_x./2 - 1./dt).*theta_stage4_present);            
            theta_stage4_z_MCFS_NPML_next = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*theta_stage4_z_MCFS_NPML_present + ...
                                                                        (eta_z./2 + 1./dt).*theta_stage4_next + (eta_z./2 - 1./dt).*theta_stage4_present);    

            dtheta_stage4_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",theta_stage4_x_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            dtheta_stage4_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",theta_stage4_z_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);       

            q_x_stage4_next_half = 0.178079954393132.*q_x_mid + 0.821920045606868.*q_x_stage3_next_half ...
                                                                                                          + 0.544974750228521*dt.*(-1./tau11_I).*(K11_I.*dtheta_stage4_x_MCFS_NPML_dx);
            q_z_stage4_next_half = 0.178079954393132.*q_z_mid + 0.821920045606868.*q_z_stage3_next_half ...
                                                                                                          + 0.544974750228521*dt.*(-1./tau33_I).*(K33_I.*dtheta_stage4_z_MCFS_NPML_dz);

            q_x_stage4_next_half = q_x_stage4_next_half.*Cerjan_Absorbing_Boundary_Condition;
            q_z_stage4_next_half = q_z_stage4_next_half.*Cerjan_Absorbing_Boundary_Condition;

            q_x_stage4_x_MCFS_NPML_next_half = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*q_x_stage4_x_MCFS_NPML_last_half + ...
                                                             (eta_x./2 + 1./dt).*q_x_stage4_next_half + (eta_x./2 - 1./dt).*q_x_stage4_last_half);
            q_x_stage4_z_MCFS_NPML_next_half = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*q_x_stage4_z_MCFS_NPML_last_half + ...
                                                             (eta_z./2 + 1./dt).*q_x_stage4_next_half + (eta_z./2 - 1./dt).*q_x_stage4_last_half);  
            q_z_stage4_x_MCFS_NPML_next_half = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*q_z_stage4_x_MCFS_NPML_last_half + ...
                                                             (eta_x./2 + 1./dt).*q_z_stage4_next_half + (eta_x./2 - 1./dt).*q_z_stage4_last_half);
            q_z_stage4_z_MCFS_NPML_next_half = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*q_z_stage4_z_MCFS_NPML_last_half + ...
                                                             (eta_z./2 + 1./dt).*q_z_stage4_next_half + (eta_z./2 - 1./dt).*q_z_stage4_last_half);                 

            dq_x_stage4_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",q_x_stage4_x_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            dq_z_stage4_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",q_z_stage4_z_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);


            T_PK1_xx_stage4_next = 0.178079954393132.*T_PK1_xx_present + 0.821920045606868.*T_PK1_xx_stage3_next ... 
                                                                                                                    + 0.544974750228521*(dt).*(A11_adi_I.*dv_x_stage4_x_MCFS_NPML_dx + A13_adi_I.*dv_z_stage4_z_MCFS_NPML_dz + ...
                                                                                                                                                                   A18_iso_I.*dv_x_stage4_z_MCFS_NPML_dz + A15_iso_I.*dv_z_stage4_x_MCFS_NPML_dx + ...
                                                                                                                                                                   (-M11_eqv_I./(Rho_I.*cV_I)).*(dq_x_stage4_x_MCFS_NPML_dx + dq_z_stage4_z_MCFS_NPML_dz));     

            T_PK1_zz_stage4_next = 0.178079954393132.*T_PK1_zz_present + 0.821920045606868.*T_PK1_zz_stage3_next ... 
                                                                                                                    + 0.544974750228521*(dt).*(A13_adi_I.*dv_x_stage4_x_MCFS_NPML_dx + A33_adi_I.*dv_z_stage4_z_MCFS_NPML_dz + ...
                                                                                                                                                                   A38_iso_I.*dv_x_stage4_z_MCFS_NPML_dz + A35_iso_I.*dv_z_stage4_x_MCFS_NPML_dx + ...
                                                                                                                                                                   (-M33_eqv_I./(Rho_I.*cV_I)).*(dq_x_stage4_x_MCFS_NPML_dx + dq_z_stage4_z_MCFS_NPML_dz));      

            T_PK1_xz_stage4_next = 0.178079954393132.*T_PK1_xz_present + 0.821920045606868.*T_PK1_xz_stage3_next ... 
                                                                                                                    + 0.544974750228521*(dt).*(A18_iso_I.*dv_x_stage4_x_MCFS_NPML_dx + A38_iso_I.*dv_z_stage4_z_MCFS_NPML_dz + ...
                                                                                                                                                                   A88_iso_I.*dv_x_stage4_z_MCFS_NPML_dz + A58_iso_I.*dv_z_stage4_x_MCFS_NPML_dx);                 
            
            T_PK1_zx_stage4_next = 0.178079954393132.*T_PK1_zx_present + 0.821920045606868.*T_PK1_zx_stage3_next ... 
                                                                                                                    + 0.544974750228521*(dt).*(A15_iso_I.*dv_x_stage4_x_MCFS_NPML_dx + A35_iso_I.*dv_z_stage4_z_MCFS_NPML_dz + ...
                                                                                                                                                                   A58_iso_I.*dv_x_stage4_z_MCFS_NPML_dz + A55_iso_I.*dv_z_stage4_x_MCFS_NPML_dx);                
            
            T_PK1_xx_stage4_next = T_PK1_xx_stage4_next.*Cerjan_Absorbing_Boundary_Condition;
            T_PK1_zz_stage4_next = T_PK1_zz_stage4_next.*Cerjan_Absorbing_Boundary_Condition;
            T_PK1_xz_stage4_next = T_PK1_xz_stage4_next.*Cerjan_Absorbing_Boundary_Condition;
            T_PK1_zx_stage4_next = T_PK1_zx_stage4_next.*Cerjan_Absorbing_Boundary_Condition;

            T_PK1_xx_stage4_x_MCFS_NPML_next = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*...
                                                                                ((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*T_PK1_xx_stage4_x_MCFS_NPML_present + ...
                                                                                 (eta_x./2 + 1./dt).*T_PK1_xx_stage4_next + (eta_x./2 - 1./dt).*T_PK1_xx_stage4_present); 
            T_PK1_zz_stage4_z_MCFS_NPML_next = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*...
                                                                                ((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*T_PK1_zz_stage4_z_MCFS_NPML_present + ...
                                                                                 (eta_z./2 + 1./dt).*T_PK1_zz_stage4_next + (eta_z./2 - 1./dt).*T_PK1_zz_stage4_present);  
            T_PK1_xz_stage4_z_MCFS_NPML_next = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*...
                                                                                ((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*T_PK1_xz_stage4_z_MCFS_NPML_present + ...
                                                                                 (eta_z./2 + 1./dt).*T_PK1_xz_stage4_next + (eta_z./2 - 1./dt).*T_PK1_xz_stage4_present);  
            T_PK1_zx_stage4_x_MCFS_NPML_next = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*...
                                                                                ((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*T_PK1_zx_stage4_x_MCFS_NPML_present + ...
                                                                                 (eta_x./2 + 1./dt).*T_PK1_zx_stage4_next + (eta_x./2 - 1./dt).*T_PK1_zx_stage4_present);  

%%%%%% 阶段5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            dT_PK1_xx_stage4_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",T_PK1_xx_stage4_x_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            dT_PK1_zz_stage4_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",T_PK1_zz_stage4_z_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            dT_PK1_xz_stage4_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",T_PK1_xz_stage4_z_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            dT_PK1_zx_stage4_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",T_PK1_zx_stage4_x_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);

            v_x_next_half =  0.517231671970585.*v_x_stage2_last_half + 0.096059710526147.*v_x_stage3_next_half ...
                                                                                                           + 0.063692468666290*(dt).*(1./Rho_I).*(dT_PK1_xx_stage3_x_MCFS_NPML_dx + dT_PK1_xz_stage3_z_MCFS_NPML_dz) + ...
                                                                                                           + 0.386708617503269.*v_x_stage4_next_half ...
                                                                                                           + 0.226007483236906*(dt).*(1./Rho_I).*(dT_PK1_xx_stage4_x_MCFS_NPML_dx + dT_PK1_xz_stage4_z_MCFS_NPML_dz); 

            v_z_next_half =  0.517231671970585.*v_z_stage2_last_half + 0.096059710526147.*v_z_stage3_next_half ...
                                                                                                           + 0.063692468666290*(dt).*(1./Rho_I).*(dT_PK1_zx_stage3_x_MCFS_NPML_dx + dT_PK1_zz_stage3_z_MCFS_NPML_dz) + ...
                                                                                                           + 0.386708617503269.*v_z_stage4_next_half ...
                                                                                                           + 0.226007483236906*(dt).*(1./Rho_I).*(dT_PK1_zx_stage4_x_MCFS_NPML_dx + dT_PK1_zz_stage4_z_MCFS_NPML_dz);                                                                                                            

            v_x_next_half = v_x_next_half.*Cerjan_Absorbing_Boundary_Condition;
            v_z_next_half = v_z_next_half.*Cerjan_Absorbing_Boundary_Condition;

            v_x_x_MCFS_NPML_next_half = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*v_x_x_MCFS_NPML_last_half + ...
                                                             (eta_x./2 + 1./dt).*v_x_next_half + (eta_x./2 - 1./dt).*v_x_last_half);
            v_x_z_MCFS_NPML_next_half = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*v_x_z_MCFS_NPML_last_half + ...
                                                             (eta_z./2 + 1./dt).*v_x_next_half + (eta_z./2 - 1./dt).*v_x_last_half);    
            v_z_x_MCFS_NPML_next_half = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*v_z_x_MCFS_NPML_last_half + ...
                                                             (eta_x./2 + 1./dt).*v_z_next_half + (eta_x./2 - 1./dt).*v_z_last_half);
            v_z_z_MCFS_NPML_next_half = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*v_z_z_MCFS_NPML_last_half + ...
                                                             (eta_z./2 + 1./dt).*v_z_next_half + (eta_z./2 - 1./dt).*v_z_last_half);      

            dv_x_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",v_x_x_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            dv_z_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",v_z_z_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            dv_x_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",v_x_z_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            dv_z_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",v_z_x_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);

            theta_next = 0.517231671970585.*theta_stage2_present + 0.096059710526147.*theta_stage3_next ...
                                                                                                        + 0.063692468666290*dt.*(- (1./cV_I).*(dq_x_stage3_x_MCFS_NPML_dx + dq_z_stage3_z_MCFS_NPML_dz) ...
                                                                                                                                                    + (T_I./cV_I).*(M11_eqv_I.*dv_x_stage4_x_MCFS_NPML_dx + M33_eqv_I.*dv_z_stage4_z_MCFS_NPML_dz)) ...
                                                                                                        + 0.386708617503269.*theta_stage4_next ...
                                                                                                        + 0.226007483236906*dt.*(- (1./cV_I).*(dq_x_stage4_x_MCFS_NPML_dx + dq_z_stage4_z_MCFS_NPML_dz) ...
                                                                                                                                                    + (T_I./cV_I).*(M11_eqv_I.*dv_x_x_MCFS_NPML_dx + M33_eqv_I.*dv_z_z_MCFS_NPML_dz));

            theta_next = theta_next.*Cerjan_Absorbing_Boundary_Condition;

            theta_x_MCFS_NPML_next = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*theta_x_MCFS_NPML_present + ...
                                                                        (eta_x./2 + 1./dt).*theta_next + (eta_x./2 - 1./dt).*theta_present);            
            theta_z_MCFS_NPML_next = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*theta_z_MCFS_NPML_present + ...
                                                                        (eta_z./2 + 1./dt).*theta_next + (eta_z./2 - 1./dt).*theta_present);    

            dtheta_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",theta_x_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            dtheta_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",theta_z_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);       

            q_x_next_half = 0.517231671970585.*q_x_stage2_last_half + 0.096059710526147.*q_x_stage3_next_half ...
                                                                                                          + 0.063692468666290*dt.*(-1./tau11_I).*(K11_I.*dtheta_stage4_x_MCFS_NPML_dx) ...
                                                                                                          + 0.386708617503269.*q_x_stage4_next_half ...
                                                                                                          + 0.226007483236906*dt.*(-1./tau11_I).*(K11_I.*dtheta_x_MCFS_NPML_dx);

            q_z_next_half = 0.517231671970585.*q_z_stage2_last_half + 0.096059710526147.*q_z_stage3_next_half ...
                                                                                                          + 0.063692468666290*dt.*(-1./tau33_I).*(K33_I.*dtheta_stage4_z_MCFS_NPML_dz) ...
                                                                                                          + 0.386708617503269.*q_z_stage4_next_half ...
                                                                                                          + 0.226007483236906*dt.*(-1./tau33_I).*(K33_I.*dtheta_z_MCFS_NPML_dz);

            q_x_next_half = q_x_next_half.*Cerjan_Absorbing_Boundary_Condition;
            q_z_next_half = q_z_next_half.*Cerjan_Absorbing_Boundary_Condition;

            q_x_x_MCFS_NPML_next_half = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*q_x_x_MCFS_NPML_last_half + ...
                                                             (eta_x./2 + 1./dt).*q_x_next_half + (eta_x./2 - 1./dt).*q_x_last_half);
            q_x_z_MCFS_NPML_next_half = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*q_x_z_MCFS_NPML_last_half + ...
                                                             (eta_z./2 + 1./dt).*q_x_next_half + (eta_z./2 - 1./dt).*q_x_last_half);  
            q_z_x_MCFS_NPML_next_half = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*q_z_x_MCFS_NPML_last_half + ...
                                                             (eta_x./2 + 1./dt).*q_z_next_half + (eta_x./2 - 1./dt).*q_z_last_half);
            q_z_z_MCFS_NPML_next_half = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*q_z_z_MCFS_NPML_last_half + ...
                                                             (eta_z./2 + 1./dt).*q_z_next_half + (eta_z./2 - 1./dt).*q_z_last_half);                 

            dq_x_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",q_x_x_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            dq_z_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",q_z_z_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);


            %   刚性方程组 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            q_x_mid = (exp(-dt./tau_ave_I) - 1).*q_x_next_half;
            q_z_mid = (exp(-dt./tau_ave_I) - 1).*q_z_next_half;

            %   全局场采用Cerjan边界层进行吸收，空间导数采用PML吸收
            q_x_mid = q_x_mid.*Cerjan_Absorbing_Boundary_Condition;
            q_z_mid = q_z_mid.*Cerjan_Absorbing_Boundary_Condition;

            T_PK1_xx_next = 0.517231671970585.*T_PK1_xx_stage2_present + 0.096059710526147.*T_PK1_xx_stage3_next ... 
                                                                       + 0.063692468666290*(dt).*(A11_adi_I.*dv_x_stage4_x_MCFS_NPML_dx + A13_adi_I.*dv_z_stage4_z_MCFS_NPML_dz + ...
                                                                                                  A18_iso_I.*dv_x_stage4_z_MCFS_NPML_dz + A15_iso_I.*dv_z_stage4_x_MCFS_NPML_dx + ...
                                                                                                 (-M11_eqv_I./(Rho_I.*cV_I)).*(dq_x_stage4_x_MCFS_NPML_dx + dq_z_stage4_z_MCFS_NPML_dz)) ...
                                                                       + 0.386708617503269.*T_PK1_xx_stage4_next ...
                                                                       + 0.226007483236906*(dt).*(A11_adi_I.*dv_x_x_MCFS_NPML_dx + A13_adi_I.*dv_z_z_MCFS_NPML_dz + ...
                                                                                                  A18_iso_I.*dv_x_z_MCFS_NPML_dz + A15_iso_I.*dv_z_x_MCFS_NPML_dx + ...
                                                                                                 (-M11_eqv_I./(Rho_I.*cV_I)).*(dq_x_x_MCFS_NPML_dx + dq_z_z_MCFS_NPML_dz));     

            T_PK1_zz_next = 0.517231671970585.*T_PK1_zz_stage2_present + 0.096059710526147.*T_PK1_zz_stage3_next ... 
                                                                                                                    + 0.063692468666290*(dt).*(A13_adi_I.*dv_x_stage4_x_MCFS_NPML_dx + A33_adi_I.*dv_z_stage4_z_MCFS_NPML_dz + ...
                                                                                                                                                                   A38_iso_I.*dv_x_stage4_z_MCFS_NPML_dz + A35_iso_I.*dv_z_stage4_x_MCFS_NPML_dx + ...
                                                                                                                                                                  (-M33_eqv_I./(Rho_I.*cV_I)).*(dq_x_stage4_x_MCFS_NPML_dx + dq_z_stage4_z_MCFS_NPML_dz)) ...
                                                                                                                    + 0.386708617503269.*T_PK1_zz_stage4_next ...
                                                                                                                    + 0.226007483236906*(dt).*(A13_adi_I.*dv_x_x_MCFS_NPML_dx + A33_adi_I.*dv_z_z_MCFS_NPML_dz + ...
                                                                                                                                                                   A38_iso_I.*dv_x_z_MCFS_NPML_dz + A35_iso_I.*dv_z_x_MCFS_NPML_dx + ...
                                                                                                                                                                  (-M33_eqv_I./(Rho_I.*cV_I)).*(dq_x_x_MCFS_NPML_dx + dq_z_z_MCFS_NPML_dz));       

            T_PK1_xz_next = 0.517231671970585.*T_PK1_xz_stage2_present + 0.096059710526147.*T_PK1_xz_stage3_next ... 
                                                                                                                    + 0.063692468666290*(dt).*(A18_iso_I.*dv_x_stage4_x_MCFS_NPML_dx + A38_iso_I.*dv_z_stage4_z_MCFS_NPML_dz + ...
                                                                                                                                                                   A88_iso_I.*dv_x_stage4_z_MCFS_NPML_dz + A58_iso_I.*dv_z_stage4_x_MCFS_NPML_dx) ...
                                                                                                                    + 0.386708617503269.*T_PK1_xz_stage4_next ...
                                                                                                                    + 0.226007483236906*(dt).*(A18_iso_I.*dv_x_x_MCFS_NPML_dx + A38_iso_I.*dv_z_z_MCFS_NPML_dz + ...
                                                                                                                                                                   A88_iso_I.*dv_x_z_MCFS_NPML_dz + A58_iso_I.*dv_z_x_MCFS_NPML_dx);      

            T_PK1_zx_next = 0.517231671970585.*T_PK1_zx_stage2_present + 0.096059710526147.*T_PK1_zx_stage3_next ... 
                                                                                                                    + 0.063692468666290*(dt).*(A15_iso_I.*dv_x_stage4_x_MCFS_NPML_dx + A35_iso_I.*dv_z_stage4_z_MCFS_NPML_dz + ...
                                                                                                                                                                   A58_iso_I.*dv_x_stage4_z_MCFS_NPML_dz + A55_iso_I.*dv_z_stage4_x_MCFS_NPML_dx) ...
                                                                                                                    + 0.386708617503269.*T_PK1_xz_stage4_next ...
                                                                                                                    + 0.226007483236906*(dt).*(A15_iso_I.*dv_x_x_MCFS_NPML_dx + A35_iso_I.*dv_z_z_MCFS_NPML_dz + ...
                                                                                                                                                                   A58_iso_I.*dv_x_z_MCFS_NPML_dz + A55_iso_I.*dv_z_x_MCFS_NPML_dx);      

            T_PK1_xx_next = T_PK1_xx_next.*Cerjan_Absorbing_Boundary_Condition;
            T_PK1_zz_next = T_PK1_zz_next.*Cerjan_Absorbing_Boundary_Condition;
            T_PK1_xz_next = T_PK1_xz_next.*Cerjan_Absorbing_Boundary_Condition;
            T_PK1_zx_next = T_PK1_zx_next.*Cerjan_Absorbing_Boundary_Condition;

            T_PK1_xx_x_MCFS_NPML_next = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*...
                                                                 ((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*T_PK1_xx_x_MCFS_NPML_present + ...
                                                                 (eta_x./2 + 1./dt).*T_PK1_xx_next + (eta_x./2 - 1./dt).*T_PK1_xx_present); 
            T_PK1_zz_z_MCFS_NPML_next = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*...
                                                                 ((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*T_PK1_zz_z_MCFS_NPML_present + ...
                                                                 (eta_z./2 + 1./dt).*T_PK1_zz_next + (eta_z./2 - 1./dt).*T_PK1_zz_present);  
            T_PK1_xz_z_MCFS_NPML_next = ((beta_z./dt + (eta_z.*beta_z + alpha_z_PML_optimized)./2)).^(-1).*...
                                                                 ((beta_z./dt - (eta_z.*beta_z + alpha_z_PML_optimized)./2).*T_PK1_xz_z_MCFS_NPML_present + ...
                                                                 (eta_z./2 + 1./dt).*T_PK1_xz_next + (eta_z./2 - 1./dt).*T_PK1_xz_present);  
            T_PK1_zx_x_MCFS_NPML_next = ((beta_x./dt + (eta_x.*beta_x + alpha_x_PML_optimized)./2)).^(-1).*...
                                                                 ((beta_x./dt - (eta_x.*beta_x + alpha_x_PML_optimized)./2).*T_PK1_zx_x_MCFS_NPML_present + ...
                                                                 (eta_x./2 + 1./dt).*T_PK1_zx_next + (eta_x./2 - 1./dt).*T_PK1_zx_present);  


            phi_Helmholtz_next = phi_Helmholtz_present + 0.5.*Rho_I.*(v_x_next_half.^2 + v_z_next_half.^2) + 0.5.*(abs(T_PK1_xx_next.*(dv_x_x_MCFS_NPML_dx)) + abs(T_PK1_zz_next.*(dv_z_z_MCFS_NPML_dz)));
        
        for i = 1:nz_PML
            for j = 1:nx_PML
                phi_total(n,1) = phi_total(n,1) + phi_Helmholtz_next(i,j);
            end
        end            
         %   计算内部区域总能量
        for i = PML_l + 1: nz_PML - PML_l
            for j = PML_l + 1: nx_PML - PML_l
                phi_interior(n,1) = phi_interior(n,1) + phi_Helmholtz_next(i,j);
            end
        end       

         v_x_last_half = v_x_next_half;      v_z_last_half = v_z_next_half;        

        v_x_stage1_last_half = v_x_stage1_next_half;    v_z_stage1_last_half = v_z_stage1_next_half;    
        v_x_stage2_last_half = v_x_stage2_next_half;    v_z_stage2_last_half = v_z_stage2_next_half;    
        v_x_stage3_last_half = v_x_stage3_next_half;    v_z_stage3_last_half = v_z_stage3_next_half;    
        v_x_stage4_last_half = v_x_stage4_next_half;    v_z_stage4_last_half = v_z_stage4_next_half;    

        T_PK1_xx_last = T_PK1_xx_present;       T_PK1_xx_present = T_PK1_xx_next;
        T_PK1_xx_stage1_last = T_PK1_xx_stage1_present;     T_PK1_xx_stage1_present = T_PK1_xx_stage1_next;
        T_PK1_xx_stage2_last = T_PK1_xx_stage2_present;     T_PK1_xx_stage2_present = T_PK1_xx_stage2_next;
        T_PK1_xx_stage3_last = T_PK1_xx_stage3_present;     T_PK1_xx_stage3_present = T_PK1_xx_stage3_next;
        T_PK1_xx_stage4_last = T_PK1_xx_stage4_present;     T_PK1_xx_stage4_present = T_PK1_xx_stage4_next;

        T_PK1_zz_last = T_PK1_zz_present;       T_PK1_zz_present = T_PK1_zz_next;
        T_PK1_zz_stage1_last = T_PK1_zz_stage1_present;     T_PK1_zz_stage1_present = T_PK1_zz_stage1_next;
        T_PK1_zz_stage2_last = T_PK1_zz_stage2_present;     T_PK1_zz_stage2_present = T_PK1_zz_stage2_next;
        T_PK1_zz_stage3_last = T_PK1_zz_stage3_present;     T_PK1_zz_stage3_present = T_PK1_zz_stage3_next;
        T_PK1_zz_stage4_last = T_PK1_zz_stage4_present;     T_PK1_zz_stage4_present = T_PK1_zz_stage4_next;

        T_PK1_xz_last = T_PK1_xz_present;       T_PK1_xz_present = T_PK1_xz_next;
        T_PK1_xz_stage1_last = T_PK1_xz_stage1_present;     T_PK1_xz_stage1_present = T_PK1_xz_stage1_next;
        T_PK1_xz_stage2_last = T_PK1_xz_stage2_present;     T_PK1_xz_stage2_present = T_PK1_xz_stage2_next;
        T_PK1_xz_stage3_last = T_PK1_xz_stage3_present;     T_PK1_xz_stage3_present = T_PK1_xz_stage3_next;
        T_PK1_xz_stage4_last = T_PK1_xz_stage4_present;     T_PK1_xz_stage4_present = T_PK1_xz_stage4_next;

        T_PK1_zx_last = T_PK1_zx_present;       T_PK1_zx_present = T_PK1_zx_next;
        T_PK1_zx_stage1_last = T_PK1_zx_stage1_present;     T_PK1_zx_stage1_present = T_PK1_zx_stage1_next;
        T_PK1_zx_stage2_last = T_PK1_zx_stage2_present;     T_PK1_zx_stage2_present = T_PK1_zx_stage2_next;
        T_PK1_zx_stage3_last = T_PK1_zx_stage3_present;     T_PK1_zx_stage3_present = T_PK1_zx_stage3_next;
        T_PK1_zx_stage4_last = T_PK1_zx_stage4_present;     T_PK1_zx_stage4_present = T_PK1_zx_stage4_next;

        q_x_last_half = q_x_next_half;
        q_z_last_half = q_z_next_half;

        q_x_stage1_last_half = q_x_stage1_next_half;
        q_z_stage1_last_half = q_z_stage1_next_half;

        q_x_stage2_last_half = q_x_stage2_next_half;
        q_z_stage2_last_half = q_z_stage2_next_half;

        q_x_stage3_last_half = q_x_stage3_next_half;
        q_z_stage3_last_half = q_z_stage3_next_half;

        q_x_stage4_last_half = q_x_stage4_next_half;
        q_z_stage4_last_half = q_z_stage4_next_half;

        theta_present = theta_next;

        theta_stage1_last = theta_stage1_present;
        theta_stage1_present = theta_stage1_next;
        theta_stage2_last = theta_stage2_present;
        theta_stage2_present = theta_stage2_next;
        theta_stage3_last = theta_stage3_present;
        theta_stage3_present = theta_stage3_next;
        theta_stage4_last = theta_stage4_present;
        theta_stage4_present = theta_stage4_next;

        v_x_x_MCFS_NPML_last_half = v_x_x_MCFS_NPML_next_half;
        v_z_z_MCFS_NPML_last_half = v_z_z_MCFS_NPML_next_half;
        v_x_z_MCFS_NPML_last_half = v_x_z_MCFS_NPML_next_half;
        v_z_x_MCFS_NPML_last_half = v_z_x_MCFS_NPML_next_half;

        v_x_stage1_x_MCFS_NPML_last_half = v_x_stage1_x_MCFS_NPML_next_half;
        v_z_stage1_z_MCFS_NPML_last_half = v_z_stage1_z_MCFS_NPML_next_half;
        v_x_stage1_z_MCFS_NPML_last_half = v_x_stage1_z_MCFS_NPML_next_half;
        v_z_stage1_x_MCFS_NPML_last_half = v_z_stage1_x_MCFS_NPML_next_half;

        v_x_stage2_x_MCFS_NPML_last_half = v_x_stage2_x_MCFS_NPML_next_half;
        v_z_stage2_z_MCFS_NPML_last_half = v_z_stage2_z_MCFS_NPML_next_half;
        v_x_stage2_z_MCFS_NPML_last_half = v_x_stage2_z_MCFS_NPML_next_half;
        v_z_stage2_x_MCFS_NPML_last_half = v_z_stage2_x_MCFS_NPML_next_half;

        v_x_stage3_x_MCFS_NPML_last_half = v_x_stage3_x_MCFS_NPML_next_half;
        v_z_stage3_z_MCFS_NPML_last_half = v_z_stage3_z_MCFS_NPML_next_half;
        v_x_stage3_z_MCFS_NPML_last_half = v_x_stage3_z_MCFS_NPML_next_half;
        v_z_stage3_x_MCFS_NPML_last_half = v_z_stage3_x_MCFS_NPML_next_half;

        v_x_stage4_x_MCFS_NPML_last_half = v_x_stage4_x_MCFS_NPML_next_half;
        v_z_stage4_z_MCFS_NPML_last_half = v_z_stage4_z_MCFS_NPML_next_half;
        v_x_stage4_z_MCFS_NPML_last_half = v_x_stage4_z_MCFS_NPML_next_half;
        v_z_stage4_x_MCFS_NPML_last_half = v_z_stage4_x_MCFS_NPML_next_half;

        T_PK1_xx_stage1_x_MCFS_NPML_last = T_PK1_xx_stage1_x_MCFS_NPML_present;
        T_PK1_xx_stage1_x_MCFS_NPML_present = T_PK1_xx_stage1_x_MCFS_NPML_next;

        T_PK1_xx_stage2_x_MCFS_NPML_last = T_PK1_xx_stage2_x_MCFS_NPML_present;
        T_PK1_xx_stage2_x_MCFS_NPML_present = T_PK1_xx_stage2_x_MCFS_NPML_next;

        T_PK1_xx_stage3_x_MCFS_NPML_last = T_PK1_xx_stage3_x_MCFS_NPML_present;
        T_PK1_xx_stage3_x_MCFS_NPML_present = T_PK1_xx_stage3_x_MCFS_NPML_next;

        T_PK1_xx_stage4_x_MCFS_NPML_last = T_PK1_xx_stage4_x_MCFS_NPML_present;
        T_PK1_xx_stage4_x_MCFS_NPML_present = T_PK1_xx_stage4_x_MCFS_NPML_next;

        T_PK1_xx_x_MCFS_NPML_last = T_PK1_xx_x_MCFS_NPML_present;
        T_PK1_xx_x_MCFS_NPML_present = T_PK1_xx_x_MCFS_NPML_next;

        T_PK1_zz_stage1_z_MCFS_NPML_last = T_PK1_zz_stage1_z_MCFS_NPML_present;
        T_PK1_zz_stage1_z_MCFS_NPML_present = T_PK1_zz_stage1_z_MCFS_NPML_next;

        T_PK1_zz_stage2_z_MCFS_NPML_last = T_PK1_zz_stage2_z_MCFS_NPML_present;
        T_PK1_zz_stage2_z_MCFS_NPML_present = T_PK1_zz_stage2_z_MCFS_NPML_next;

        T_PK1_zz_stage3_z_MCFS_NPML_last = T_PK1_zz_stage3_z_MCFS_NPML_present;
        T_PK1_zz_stage3_z_MCFS_NPML_present = T_PK1_zz_stage3_z_MCFS_NPML_next;

        T_PK1_zz_stage4_z_MCFS_NPML_last = T_PK1_zz_stage4_z_MCFS_NPML_present;
        T_PK1_zz_stage4_z_MCFS_NPML_present = T_PK1_zz_stage4_z_MCFS_NPML_next;

        T_PK1_zz_z_MCFS_NPML_last = T_PK1_zz_z_MCFS_NPML_present;
        T_PK1_zz_z_MCFS_NPML_present = T_PK1_zz_z_MCFS_NPML_next;

        T_PK1_xz_z_MCFS_NPML_last = T_PK1_xz_z_MCFS_NPML_present;
        T_PK1_xz_z_MCFS_NPML_present = T_PK1_xz_z_MCFS_NPML_next;

        T_PK1_xz_stage1_z_MCFS_NPML_last = T_PK1_xz_stage1_z_MCFS_NPML_present;
        T_PK1_xz_stage1_z_MCFS_NPML_present = T_PK1_xz_stage1_z_MCFS_NPML_next;

        T_PK1_xz_stage2_z_MCFS_NPML_last = T_PK1_xz_stage2_z_MCFS_NPML_present;
        T_PK1_xz_stage2_z_MCFS_NPML_present = T_PK1_xz_stage2_z_MCFS_NPML_next;

        T_PK1_xz_stage3_z_MCFS_NPML_last = T_PK1_xz_stage3_z_MCFS_NPML_present;
        T_PK1_xz_stage3_z_MCFS_NPML_present = T_PK1_xz_stage3_z_MCFS_NPML_next;

        T_PK1_xz_stage4_z_MCFS_NPML_last = T_PK1_xz_stage4_z_MCFS_NPML_present;
        T_PK1_xz_stage4_z_MCFS_NPML_present = T_PK1_xz_stage4_z_MCFS_NPML_next;

        T_PK1_zx_x_MCFS_NPML_last = T_PK1_zx_x_MCFS_NPML_present;
        T_PK1_zx_x_MCFS_NPML_present = T_PK1_zx_x_MCFS_NPML_next;

        T_PK1_zx_stage1_x_MCFS_NPML_last = T_PK1_zx_stage1_x_MCFS_NPML_present;
        T_PK1_zx_stage1_x_MCFS_NPML_present = T_PK1_zx_stage1_x_MCFS_NPML_next;

        T_PK1_zx_stage2_x_MCFS_NPML_last = T_PK1_zx_stage2_x_MCFS_NPML_present;
        T_PK1_zx_stage2_x_MCFS_NPML_present = T_PK1_zx_stage2_x_MCFS_NPML_next;

        T_PK1_zx_stage3_x_MCFS_NPML_last = T_PK1_zx_stage3_x_MCFS_NPML_present;
        T_PK1_zx_stage3_x_MCFS_NPML_present = T_PK1_zx_stage3_x_MCFS_NPML_next;

        T_PK1_zx_stage4_x_MCFS_NPML_last = T_PK1_zx_stage4_x_MCFS_NPML_present;
        T_PK1_zx_stage4_x_MCFS_NPML_present = T_PK1_zx_stage4_x_MCFS_NPML_next;

        q_x_x_MCFS_NPML_last_half = q_x_x_MCFS_NPML_next_half;
        q_z_z_MCFS_NPML_last_half = q_z_z_MCFS_NPML_next_half;

        q_x_stage1_x_MCFS_NPML_last_half = q_x_stage1_x_MCFS_NPML_next_half;
        q_z_stage1_z_MCFS_NPML_last_half = q_z_stage1_z_MCFS_NPML_next_half;

        q_x_stage2_x_MCFS_NPML_last_half = q_x_stage2_x_MCFS_NPML_next_half;
        q_z_stage2_z_MCFS_NPML_last_half = q_z_stage2_z_MCFS_NPML_next_half;

        q_x_stage3_x_MCFS_NPML_last_half = q_x_stage3_x_MCFS_NPML_next_half;
        q_z_stage3_z_MCFS_NPML_last_half = q_z_stage3_z_MCFS_NPML_next_half;

        q_x_stage4_x_MCFS_NPML_last_half = q_x_stage4_x_MCFS_NPML_next_half;
        q_z_stage4_z_MCFS_NPML_last_half = q_z_stage4_z_MCFS_NPML_next_half;

        theta_x_MCFS_NPML_present = theta_x_MCFS_NPML_next;
        theta_z_MCFS_NPML_present = theta_z_MCFS_NPML_next;

        theta_stage1_x_MCFS_NPML_last = theta_stage1_x_MCFS_NPML_present;
        theta_stage1_x_MCFS_NPML_present = theta_stage1_x_MCFS_NPML_next;
        theta_stage2_x_MCFS_NPML_last = theta_stage2_x_MCFS_NPML_present;
        theta_stage2_x_MCFS_NPML_present = theta_stage2_x_MCFS_NPML_next;
        theta_stage3_x_MCFS_NPML_last = theta_stage3_x_MCFS_NPML_present;
        theta_stage3_x_MCFS_NPML_present = theta_stage3_x_MCFS_NPML_next;
        theta_stage4_x_MCFS_NPML_last = theta_stage4_x_MCFS_NPML_present;
        theta_stage4_x_MCFS_NPML_present = theta_stage4_x_MCFS_NPML_next;

        theta_stage1_z_MCFS_NPML_last = theta_stage1_z_MCFS_NPML_present;
        theta_stage1_z_MCFS_NPML_present = theta_stage1_z_MCFS_NPML_next;
        theta_stage2_z_MCFS_NPML_last = theta_stage2_z_MCFS_NPML_present;
        theta_stage2_z_MCFS_NPML_present = theta_stage2_z_MCFS_NPML_next;
        theta_stage3_z_MCFS_NPML_last = theta_stage3_z_MCFS_NPML_present;
        theta_stage3_z_MCFS_NPML_present = theta_stage3_z_MCFS_NPML_next;
        theta_stage4_z_MCFS_NPML_last = theta_stage4_z_MCFS_NPML_present;
        theta_stage4_z_MCFS_NPML_present = theta_stage4_z_MCFS_NPML_next;

        phi_Helmholtz_present = phi_Helmholtz_next;
 
       %   合成地震记录 
        %   水平布置 （布置在地表）  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if synthetic_seismogram_flag == "true"
           seismogram_1_v_x(n,:) = v_x_last_half(PML_l + 1,:);
           seismogram_1_v_z(n,:) = v_z_last_half(PML_l + 1,:);
           thermogram_1_theta(n,:) = theta_present(PML_l + 1,:);
            %   垂直布置 VSP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           seismogram_1_v_x_DAS_VSP(n,:) = v_x_last_half(:,(7000/dx) + 1 + PML_l);
           seismogram_1_v_z_DAS_VSP(n,:) = v_z_last_half(:,(7000/dx) + 1 + PML_l);
           thermogram_1_theta_DTS_VSP(n,:) = theta_present(:,(7000/dx) + 1 + PML_l);
       end

       if time_record_flag == "true"
            %   时间记录
           time_record_v_z_fz_1500ms(n,1) = t;
           time_record_v_z_fz_1500ms(n,2) = v_z_last_half((3000/dz) + 1 + PML_l,(7000/dx) + 1 + PML_l);
           time_record_v_z_fz_1500ms(n,3) = 0;
           time_record_v_x_fz_1500ms(n,1) = t;
           time_record_v_x_fz_1500ms(n,2) = v_x_last_half((3000/dz) + 1 + PML_l,(7000/dx) + 1 + PML_l);
           time_record_v_x_fz_1500ms(n,3) = 0;
           time_record_q_z_fz_1500ms(n,1) = t;
           time_record_q_z_fz_1500ms(n,2) = q_z_last_half((3000/dz) + 1 + PML_l,(7000/dx) + 1 + PML_l);
           time_record_q_z_fz_1500ms(n,3) = 0;
           time_record_q_x_fz_1500ms(n,1) = t;
           time_record_q_x_fz_1500ms(n,2) = q_x_last_half((3000/dz) + 1 + PML_l,(7000/dx) + 1 + PML_l);
           time_record_q_x_fz_1500ms(n,3) = 0;
           time_record_theta_fz_1500ms(n,1) = t;
           time_record_theta_fz_1500ms(n,2) = theta_present((3000/dz) + 1 + PML_l,(7000/dx) + 1 + PML_l);
           time_record_theta_fz_1500ms(n,3) = 0; 
       end

    %   当时间进行到0.325秒时输出快照和地震记录
       if t == t_wavefield
            %   质点速度垂直分量的波场快照    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure('Name','TAE_vz_explosive','Units','centimeters','Position',[10,10,7.3,5.45]);
            surf(x_main_region,z_main_region,v_z_last_half(PML_l + 1:nz_PML - PML_l,PML_l + 1:nx_PML - PML_l));
            set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
            str_title = ['v_Z, (t=',num2str(t*1e3),'ms)'];
            title(str_title,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
            shading interp;
            xlabel('X(m)');
            ylabel('Z(m)');
            axis image;
            %             axis([0 dx*(256 - 4*PML_l - 1) 0 dz*(256 - 4*PML_l - 1)]); 
            colorbar;
            colormap(gray);
            strings = {'(a)'};
            annotation('textbox',[0.2,0.4,0.25,0.5],'LineStyle','none','LineWidth',2,'String',strings,'Fontname','Times New Roman','fontweight','bold');

%             %   质点速度水平分量的波场快照   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure('Name','TAE_wavefield_snapshot_vx_MCFS_NPML','Units','centimeters','Position',[10,10,7.3,5.45]);
            surf(x_main_region,z_main_region,v_x_last_half(PML_l + 1:nz_PML - PML_l,PML_l + 1:nx_PML - PML_l));
            set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
            str_title = ['v_X, (t=',num2str(t*1e3),'ms)'];
            title(str_title,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
            shading interp;
            xlabel('X(m)');
            ylabel('Z(m)');
            axis image;
            %             axis([0 dx*(256 - 4*PML_l - 1) 0 dz*(256 - 4*PML_l - 1)]); 
            colorbar;
            colormap(gray);            
            strings = {'(b)'};
            annotation('textbox',[0.2,0.4,0.25,0.5],'LineStyle','none','LineWidth',2,'String',strings,'Fontname','Times New Roman','fontweight','bold');

        elseif t == t_total
            if time_record_flag == "true"
                writematrix(time_record_v_z_fz_1500ms,['C:/Users/Yuanxie Li/OneDrive/Postdoc/Time_Record/TAE_vz_fz_time_record_',...
                                                                        num2str(n),'ms_f0_',num2str(f0),'Hz_',prestress_mode,'_',direction,'_P_',num2str(P(1,1)/1e6),'_MPa_alpha_',num2str(Alpha11_I(1,1))],'FileType','text');
                writematrix(time_record_v_x_fz_1500ms,['C:/Users/Yuanxie Li/OneDrive/Postdoc/Time_Record/TAE_vx_fz_time_record_',...
                                                                        num2str(n),'ms_f0_',num2str(f0),'Hz_',prestress_mode,'_',direction,'_P_',num2str(P(1,1)/1e6),'_MPa_alpha_',num2str(Alpha11_I(1,1))],'FileType','text');
                writematrix(time_record_q_z_fz_1500ms,['C:/Users/Yuanxie Li/OneDrive/Postdoc/Time_Record/TAE_qz_fz_time_record_',...
                                                                        num2str(n),'ms_f0_',num2str(f0),'Hz_',prestress_mode,'_',direction,'_P_',num2str(P(1,1)/1e6),'_MPa_alpha_',num2str(Alpha11_I(1,1))],'FileType','text');
                writematrix(time_record_q_x_fz_1500ms,['C:/Users/Yuanxie Li/OneDrive/Postdoc/Time_Record/TAE_qx_fz_time_record_',...
                                                                        num2str(n),'ms_f0_',num2str(f0),'Hz_',prestress_mode,'_',direction,'_P_',num2str(P(1,1)/1e6),'_MPa_alpha_',num2str(Alpha11_I(1,1))],'FileType','text');
                writematrix(time_record_theta_fz_1500ms,['C:/Users/Yuanxie Li/OneDrive/Postdoc/Time_Record/TAE_theta_fz_time_record_',...
                                                                        num2str(n),'ms_f0_',num2str(f0),'Hz_',prestress_mode,'_',direction,'_P_',num2str(P(1,1)/1e6),'_MPa_alpha_',num2str(Alpha11_I(1,1))],'FileType','text');
                writematrix(phi_total,['C:/Users/Yuanxie Li/OneDrive/Postdoc/Time_Record/TAE_phi_total_fz_time_record_',...
                                                                        num2str(n),'ms_f0_',num2str(f0),'Hz_',prestress_mode,'_',direction,'_P_',num2str(P(1,1)/1e6),'_MPa_alpha_',num2str(Alpha11_I(1,1))],'FileType','text');
                writematrix(phi_interior,['C:/Users/Yuanxie Li/OneDrive/Postdoc/Time_Record/TAE_phi_interior_fz_time_record_',...
                                                                    num2str(n),'ms_f0_',num2str(f0),'Hz_',prestress_mode,'_',direction,'_P_',num2str(P(1,1)/1e6),'_MPa_alpha_',num2str(Alpha11_I(1,1))],'FileType','text');
       
            end
       end

end       

toc;
feature('memstats');

if synthetic_seismogram_flag == "true"
    % figure(3)
    figure('Name','TAE_seismogram_vz_MCFS_NPML','Units','centimeters','Position',[10,10,7.3,5.45]);
    imagesc(x_main_region,tt,seismogram_1_v_z(:,PML_l + 1:nx_PML - PML_l));
    set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    str_title = ['Vz(P=',num2str(P(1,1)/1e6), 'MPa, t=',num2str(tt*1e3),'ms)'];
    title(str_title,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    shading interp;
    xlabel('Distance (m)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    ylabel('t(s)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    colorbar;
    colormap(gray);
    strings = {'(e)'};
    annotation('textbox',[0.2,0.4,0.25,0.5],'LineStyle','none','LineWidth',2,'String',strings,'Fontname','Times New Roman','fontweight','bold');
    
    figure('Name','TAE_seismogram_vx_MCFS_NPML','Units','centimeters','Position',[10,10,7.3,5.45]);
    imagesc(x_main_region,tt,seismogram_1_v_x(:,PML_l + 1:nx_PML - PML_l));
    set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    str_title = ['Vx(P=',num2str(P(1,1)/1e6), 'MPa, t=',num2str(tt*1e3),'ms)'];
    title(str_title,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    shading interp;
    xlabel('Distance (m)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    ylabel('t(s)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    colorbar;
    colormap(gray);
    strings = {'(f)'};
    annotation('textbox',[0.2,0.4,0.25,0.5],'LineStyle','none','LineWidth',2,'String',strings,'Fontname','Times New Roman','fontweight','bold');
    
    figure('Name','TAE_thermogram_theta_MCFS_NPML','Units','centimeters','Position',[10,10,7.3,5.45]);
    imagesc(x_main_region,tt,thermogram_1_theta(:,PML_l + 1:nx_PML - PML_l));
    set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    str_title = ['\theta(P=',num2str(P(1,1)/1e6), 'MPa, t=',num2str(tt*1e3),'ms)'];
    title(str_title,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    shading interp;
    xlabel('Distance (m)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    ylabel('t(s)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    colorbar;
    colormap(gray);
    strings = {'(g)'};
    annotation('textbox',[0.2,0.4,0.25,0.5],'LineStyle','none','LineWidth',2,'String',strings,'Fontname','Times New Roman','fontweight','bold');
    
    
    figure('Name','VSP_DAS_seismogram_vx_MCFS_NPML','Units','centimeters','Position',[10,10,7.3,5.45]);
    imagesc(tt,z_main_region,seismogram_1_v_x_DAS_VSP(:,PML_l + 1:nz_PML - PML_l)');
    set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    str_title = ['VSP-DAS Vx(P=',num2str(P(1,1)/1e6), 'MPa, t=',num2str(tt*1e3),'ms)'];
    title(str_title,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    shading interp;
    xlabel('t(s)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    ylabel('Depth(m)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    colorbar;
    colormap(gray);
    strings = {'(c)'};
    annotation('textbox',[0.2,0.4,0.25,0.5],'LineStyle','none','LineWidth',2,'String',strings,'Fontname','Times New Roman','fontweight','bold');
    
    figure('Name','VSP_DTS_thermogram_\theta_MCFS_NPML','Units','centimeters','Position',[10,10,7.3,5.45]);
    imagesc(tt,z_main_region,thermogram_1_theta_DTS_VSP(:,PML_l + 1:nz_PML - PML_l)');
    set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    str_title = ['VSP-DTS \theta(P=',num2str(P(1,1)/1e6), 'MPa, t=',num2str(tt*1e3),'ms)'];
    title(str_title,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    shading interp;
    xlabel('t(s)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    ylabel('Depth(m)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    colorbar;
    colormap(gray);
    strings = {'(d)'};
    annotation('textbox',[0.2,0.4,0.25,0.5],'LineStyle','none','LineWidth',2,'String',strings,'Fontname','Times New Roman','fontweight','bold');
end


if time_record_flag == "true"
    figure('Name','TAE_vz_time_1500ms_record_875e-5m_875e-5m','Units','centimeters','Position',[10,10,7.3,5.45]);
    plot(time_record_v_z_fz_1500ms(:,1),time_record_v_z_fz_1500ms(:,2),'m','LineWidth',0.5);  %   ,'MarkerSize',3
    hold off;
    set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    xlabel('t(s)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    ylabel('Amplitude','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    if media == "two layers"
        xlim([0.5 1.5]);
    else
        xlim([0 t_wavefield]);
    end
    title('vz(1710m,330m)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    strings = {'(a)'};
    annotation('textbox',[0.2,0.4,0.25,0.5],'LineStyle','none','LineWidth',2,'String',strings,'Fontname','Times New Roman','fontweight','bold');
    
    figure('Name','TAE_vx_time_1500ms_record_875e-5m_875e-5m','Units','centimeters','Position',[10,10,7.3,5.45]);
    plot(time_record_v_x_fz_1500ms(:,1),time_record_v_x_fz_1500ms(:,2),'m','LineWidth',0.5);  %   ,'MarkerSize',3
    hold off;
    set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    xlabel('t(s)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    ylabel('Amplitude','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    if media == "two layers"
        xlim([0.5 1.5]);
    else
        xlim([0 t_wavefield]);
    end
    title('vx(1710m,330m)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    strings = {'(a)'};
    annotation('textbox',[0.2,0.4,0.25,0.5],'LineStyle','none','LineWidth',2,'String',strings,'Fontname','Times New Roman','fontweight','bold');
    
    
    figure('Name','TAE_qz_time_1500ms_record_875e-5m_875e-5m','Units','centimeters','Position',[10,10,7.3,5.45]);
    plot(time_record_q_z_fz_1500ms(:,1),time_record_q_z_fz_1500ms(:,2),'m','LineWidth',0.5);  %   ,'MarkerSize',3
    hold off;
    set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    xlabel('t(s)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    ylabel('Amplitude','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    if media == "two layers"
        xlim([0.5 1.5]);
    else
        xlim([0 t_wavefield]);
    end
    title('qz(1710m,330m)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    strings = {'(d)'};
    annotation('textbox',[0.2,0.4,0.25,0.5],'LineStyle','none','LineWidth',2,'String',strings,'Fontname','Times New Roman','fontweight','bold');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure('Name','TAE_qx_time_1500ms_record_875e-5m_875e-5m','Units','centimeters','Position',[10,10,7.3,5.45]);
    plot(time_record_q_x_fz_1500ms(:,1),time_record_q_x_fz_1500ms(:,2),'m','LineWidth',0.5);  %   ,'MarkerSize',3
    hold off;
    set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    xlabel('t(s)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    ylabel('Amplitude','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    if media == "two layers"
        xlim([0.5 1.5]);
    else
        xlim([0 t_wavefield]);
    end
    title('qx(1710m,330m)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    strings = {'(d)'};
    annotation('textbox',[0.2,0.4,0.25,0.5],'LineStyle','none','LineWidth',2,'String',strings,'Fontname','Times New Roman','fontweight','bold');
    
    figure('Name','TAE_theta_time_1500ms_record_875e-5m_875e-5m','Units','centimeters','Position',[10,10,7.3,5.45]);
    plot(time_record_theta_fz_1500ms(:,1),time_record_theta_fz_1500ms(:,2),'m','LineWidth',0.5);  %   ,'MarkerSize',3
    hold off;
    set(gca,'Fontname','Times New Roman','fontweight','bold','FontSize',10);
    xlabel('t(s)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    ylabel('Amplitude','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    if media == "two layers"
        xlim([0.5 1.5]);
    else
        xlim([0 t_wavefield]);
    end
    title('\theta(1710m,330m)','Fontname','Times New Roman','fontweight','bold','FontSize',10);
    strings = {'(d)'};
    annotation('textbox',[0.2,0.4,0.25,0.5],'LineStyle','none','LineWidth',2,'String',strings,'Fontname','Times New Roman','fontweight','bold');
end









