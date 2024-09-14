function [eta_x,eta_z,beta_x,beta_z,...
          alpha_x_PML,alpha_z_PML,...
          alpha_x_PML_optimized,alpha_z_PML_optimized] = MCFS_NPML_Boundary(dx,dz,dt,nz_PML,nx_PML,PML_l,PML_x,PML_z,...
                                                                            f0,pole,eta_0,beta_0_x,beta_0_z,...
                                                                            reflection_coefficient,n_power,Vmin, ...
                                                                            P_eta,P_beta,gamma_damping,delta_damping,P_stable,attenuation_type,media)

            beta_z = ones(nz_PML,nx_PML);      beta_x = ones(nz_PML,nx_PML);
            eta_z = zeros(nz_PML,nx_PML);       eta_x = zeros(nz_PML,nx_PML);
            K_main_damping_x = zeros(nz_PML,nx_PML);      K_main_damping_z = zeros(nz_PML,nx_PML);
            alpha_x_PML_basic = zeros(nz_PML,nx_PML);     alpha_x_PML_exp = zeros(nz_PML,nx_PML);   
            alpha_z_PML_basic = zeros(nz_PML,nx_PML);     alpha_z_PML_exp = zeros(nz_PML,nx_PML);
            alpha_x_PML = zeros(nz_PML,nx_PML);       alpha_z_PML = zeros(nz_PML,nx_PML);
            alpha_x_PML_optimized = zeros(nz_PML,nx_PML);     alpha_z_PML_optimized = zeros(nz_PML,nx_PML);
            b_x_CPML = zeros(nz_PML,nx_PML);      b_z_CPML = zeros(nz_PML,nx_PML);
            a_x_CPML = zeros(nz_PML,nx_PML);      a_z_CPML = zeros(nz_PML,nx_PML);  

    for i = 1:nz_PML
        for j = 1:nx_PML
            %%%%%%% 閸栧搫鐓?1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%   闂勩倓绨＄粻�??褰囩悰鏉垮櫤閸戣姤鏆? 鏉╂ɑ鐪伴崣鏍︾啊cpml鐠侊紕鐣婚崡椋幮濇い閫涜厬閻ㄥ嫮閮撮弫鐧產x aaz bbx bbz %%%%%%%%
    
            if i >= 1 && i <= PML_l && j >= 1 && j <= PML_l             %   ���Ͻ�

                z_PML = dz*(PML_l - i);     x_PML = dx*(PML_l - j);     
                eta_x(i,j) = eta_0*pi*f0*(1 - (x_PML/PML_x).^P_eta);
                eta_z(i,j) = eta_0*pi*f0*(1 - (z_PML/PML_z).^P_eta);

                beta_x(i,j) = 1 + (beta_0_x - 1).*(x_PML/PML_x).^P_beta;    % beta_0閸欐�???1-10
                beta_z(i,j) = 1 + (beta_0_z - 1).*(z_PML/PML_z).^P_beta;



                K_main_damping_x(i,j) = log(1./reflection_coefficient).*((n_power+1).*Vmin(i,j)./(2*PML_x));
                K_main_damping_z(i,j) = log(1./reflection_coefficient).*((n_power+1).*Vmin(i,j)./(2*PML_z));
    
                if attenuation_type == "exp"
                    alpha_x_PML_basic(i,j) = (x_PML/PML_x).^n_power;
                elseif attenuation_type == "trigonometric"
                    alpha_x_PML_basic(i,j) = 1 - cos(pi*(x_PML)/(2*PML_x));    %    余弦�?
                end
    
                alpha_x_PML_exp(i,j) = gamma_damping.*exp(-(delta_damping.*PML_x./x_PML));

                if attenuation_type == "exp"                
                    alpha_z_PML_basic(i,j) = (z_PML/PML_z).^n_power;
                elseif attenuation_type == "trigonometric"                    
                    alpha_z_PML_basic(i,j) = 1 - cos(pi*(z_PML)/(2*PML_z));
                end
    
                alpha_z_PML_exp(i,j) = gamma_damping.*exp(-(delta_damping.*PML_z./z_PML));
                %   缂佸繐鍚?鐞涙澘鍣洪崙鑺ユ�? 
                alpha_x_PML(i,j) = K_main_damping_x(i,j).*alpha_x_PML_basic(i,j);
                alpha_z_PML(i,j) = K_main_damping_z(i,j).*alpha_z_PML_basic(i,j);
                %   娴兼ê瀵叉潻鍥╂畱鐞涙澘鍣洪崙鑺ユ�?
                alpha_x_PML_optimized(i,j) = K_main_damping_x(i,j).*(alpha_x_PML_basic(i,j) + alpha_x_PML_exp(i,j));
                alpha_z_PML_optimized(i,j) = K_main_damping_z(i,j).*(alpha_z_PML_basic(i,j) + alpha_z_PML_exp(i,j));
    



                % if media == "two layers"            %%%%    ���ɱ���߽�
                %     alpha_x_PML_optimized(i,j) = 0;
                %     alpha_z_PML_optimized(i,j) = 0;
                % end

                %%%%%% CPML %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                b_x_CPML(i,j) = exp(-(alpha_x_PML(i,j)./beta_x(i,j) + eta_x(i,j)).*dt);
                b_z_CPML(i,j) = exp(-(alpha_z_PML(i,j)./beta_z(i,j) + eta_z(i,j)).*dt);
                a_x_CPML(i,j) = (alpha_x_PML(i,j)/(beta_x(i,j)*(alpha_x_PML(i,j) + beta_x(i,j)*eta_x(i,j)))).*(beta_x(i,j) - 1);
                a_z_CPML(i,j) = (alpha_z_PML(i,j)/(beta_z(i,j)*(alpha_z_PML(i,j) + beta_z(i,j)*eta_z(i,j)))).*(beta_z(i,j) - 1);
    
            elseif i >= 1 && i <= PML_l && j > nx_PML - PML_l && j <= nx_PML        %   ���Ͻ�
                z_PML = dz*(PML_l - i);     x_PML = dx*(j - (nx_PML - PML_l) - 1);     
    
                eta_x(i,j) = eta_0*pi*f0*(1 - (x_PML/PML_x).^P_eta);
                eta_z(i,j) = eta_0*pi*f0*(1 - (z_PML/PML_z).^P_eta);
                beta_x(i,j) = 1 + (beta_0_x - 1).*(x_PML/PML_x).^P_beta;    % beta_0閸欐�???1-10
                beta_z(i,j) = 1 + (beta_0_z - 1).*(z_PML/PML_z).^P_beta;
    
                K_main_damping_x(i,j) = log(1./reflection_coefficient).*((n_power+1).*Vmin(i,j)./(2*PML_x));
                K_main_damping_z(i,j) = log(1./reflection_coefficient).*((n_power+1).*Vmin(i,j)./(2*PML_z));
    
                if attenuation_type == "exp"
                    alpha_x_PML_basic(i,j) = (x_PML/PML_x).^n_power;
                elseif attenuation_type == "trigonometric"
                    alpha_x_PML_basic(i,j) = 1 - cos(pi*(x_PML)/(2*PML_x));    %    余弦�?
                end
                alpha_x_PML_exp(i,j) = gamma_damping.*exp(-(delta_damping.*PML_x./x_PML));
                if attenuation_type == "exp"                
                    alpha_z_PML_basic(i,j) = (z_PML/PML_z).^n_power;
                elseif attenuation_type == "trigonometric"                    
                    alpha_z_PML_basic(i,j) = 1 - cos(pi*(z_PML)/(2*PML_z));
                end        
                alpha_z_PML_exp(i,j) = gamma_damping.*exp(-(delta_damping.*PML_z./z_PML));
                %   缂佸繐鍚?鐞涙澘鍣洪崙鑺ユ�? 
                alpha_x_PML(i,j) = K_main_damping_x(i,j).*alpha_x_PML_basic(i,j);
                alpha_z_PML(i,j) = K_main_damping_z(i,j).*alpha_z_PML_basic(i,j);
                %   娴兼ê瀵叉潻鍥╂畱鐞涙澘鍣洪崙鑺ユ�?
                alpha_x_PML_optimized(i,j) = K_main_damping_x(i,j).*(alpha_x_PML_basic(i,j) + alpha_x_PML_exp(i,j));
                alpha_z_PML_optimized(i,j) = K_main_damping_z(i,j).*(alpha_z_PML_basic(i,j) + alpha_z_PML_exp(i,j));
    
                % if media == "two layers"            %%%%    ���ɱ���߽�
                %     alpha_x_PML_optimized(i,j) = 0;
                %     alpha_z_PML_optimized(i,j) = 0;
                % end

                %%%%%% CPML %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                b_x_CPML(i,j) = exp(-(alpha_x_PML(i,j)./beta_x(i,j) + eta_x(i,j)).*dt);
                b_z_CPML(i,j) = exp(-(alpha_z_PML(i,j)./beta_z(i,j) + eta_z(i,j)).*dt);
                a_x_CPML(i,j) = (alpha_x_PML(i,j)/(beta_x(i,j)*(alpha_x_PML(i,j) + beta_x(i,j)*eta_x(i,j)))).*(beta_x(i,j) - 1);
                a_z_CPML(i,j) = (alpha_z_PML(i,j)/(beta_z(i,j)*(alpha_z_PML(i,j) + beta_z(i,j)*eta_z(i,j)))).*(beta_z(i,j) - 1);   
    
            elseif i > nz_PML - PML_l && i <= nz_PML && j >= 1 && j <= PML_l
                z_PML = dz*(i - (nz_PML - PML_l) - 1);      x_PML = dx*(PML_l - j);     
    
                eta_x(i,j) = eta_0*pi*f0*(1 - (x_PML/PML_x).^P_eta);
                eta_z(i,j) = eta_0*pi*f0*(1 - (z_PML/PML_z).^P_eta);
                beta_x(i,j) = 1 + (beta_0_x - 1).*(x_PML/PML_x).^P_beta;    % beta_0閸欐�???1-10
                beta_z(i,j) = 1 + (beta_0_z - 1).*(z_PML/PML_z).^P_beta;
    
                K_main_damping_x(i,j) = log(1./reflection_coefficient).*((n_power+1).*Vmin(i,j)./(2*PML_x));
                K_main_damping_z(i,j) = log(1./reflection_coefficient).*((n_power+1).*Vmin(i,j)./(2*PML_z));
    
                if attenuation_type == "exp"
                    alpha_x_PML_basic(i,j) = (x_PML/PML_x).^n_power;
                elseif attenuation_type == "trigonometric"
                    alpha_x_PML_basic(i,j) = 1 - cos(pi*(x_PML)/(2*PML_x));    %    余弦�?
                end
                alpha_x_PML_exp(i,j) = gamma_damping.*exp(-(delta_damping.*PML_x./x_PML));

                if attenuation_type == "exp"                
                    alpha_z_PML_basic(i,j) = (z_PML/PML_z).^n_power;
                elseif attenuation_type == "trigonometric"                    
                    alpha_z_PML_basic(i,j) = 1 - cos(pi*(z_PML)/(2*PML_z));
                end    

                alpha_z_PML_exp(i,j) = gamma_damping.*exp(-(delta_damping.*PML_z./z_PML));
                %   缂佸繐鍚?鐞涙澘鍣洪崙鑺ユ�? 
                alpha_x_PML(i,j) = K_main_damping_x(i,j).*alpha_x_PML_basic(i,j);
                alpha_z_PML(i,j) = K_main_damping_z(i,j).*alpha_z_PML_basic(i,j);
                %   娴兼ê瀵叉潻鍥╂畱鐞涙澘鍣洪崙鑺ユ�?
                alpha_x_PML_optimized(i,j) = K_main_damping_x(i,j).*(alpha_x_PML_basic(i,j) + alpha_x_PML_exp(i,j));
                alpha_z_PML_optimized(i,j) = K_main_damping_z(i,j).*(alpha_z_PML_basic(i,j) + alpha_z_PML_exp(i,j));
    
                %%%%%% CPML %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                b_x_CPML(i,j) = exp(-(alpha_x_PML(i,j)./beta_x(i,j) + eta_x(i,j)).*dt);
                b_z_CPML(i,j) = exp(-(alpha_z_PML(i,j)./beta_z(i,j) + eta_z(i,j)).*dt);
                a_x_CPML(i,j) = (alpha_x_PML(i,j)/(beta_x(i,j)*(alpha_x_PML(i,j) + beta_x(i,j)*eta_x(i,j)))).*(beta_x(i,j) - 1);
                a_z_CPML(i,j) = (alpha_z_PML(i,j)/(beta_z(i,j)*(alpha_z_PML(i,j) + beta_z(i,j)*eta_z(i,j)))).*(beta_z(i,j) - 1);  
    
            elseif i > nz_PML - PML_l && i <= nz_PML && j > nx_PML - PML_l && j <= nx_PML
    
                z_PML = dz*(i - (nz_PML - PML_l) - 1);      x_PML = dx*(j - (nx_PML - PML_l) - 1);     
    
                eta_x(i,j) = eta_0*pi*f0*(1 - (x_PML/PML_x).^P_eta);
                eta_z(i,j) = eta_0*pi*f0*(1 - (z_PML/PML_z).^P_eta);
                beta_x(i,j) = 1 + (beta_0_x - 1).*(x_PML/PML_x).^P_beta;    % beta_0閸欐�???1-10
                beta_z(i,j) = 1 + (beta_0_z - 1).*(z_PML/PML_z).^P_beta;
    
                K_main_damping_x(i,j) = log(1./reflection_coefficient).*((n_power+1).*Vmin(i,j)./(2*PML_x));
                K_main_damping_z(i,j) = log(1./reflection_coefficient).*((n_power+1).*Vmin(i,j)./(2*PML_z));
    
                if attenuation_type == "exp"
                    alpha_x_PML_basic(i,j) = (x_PML/PML_x).^n_power;
                elseif attenuation_type == "trigonometric"
                    alpha_x_PML_basic(i,j) = 1 - cos(pi*(x_PML)/(2*PML_x));    %    余弦�?
                end

                alpha_x_PML_exp(i,j) = gamma_damping.*exp(-(delta_damping.*PML_x./x_PML));
                if attenuation_type == "exp"                
                    alpha_z_PML_basic(i,j) = (z_PML/PML_z).^n_power;
                elseif attenuation_type == "trigonometric"                    
                    alpha_z_PML_basic(i,j) = 1 - cos(pi*(z_PML)/(2*PML_z));
                end              
                alpha_z_PML_exp(i,j) = gamma_damping.*exp(-(delta_damping.*PML_z./z_PML));
                %   缂佸繐鍚?鐞涙澘鍣洪崙鑺ユ�? 
                alpha_x_PML(i,j) = K_main_damping_x(i,j).*alpha_x_PML_basic(i,j);
                alpha_z_PML(i,j) = K_main_damping_z(i,j).*alpha_z_PML_basic(i,j);
                %   娴兼ê瀵叉潻鍥╂畱鐞涙澘鍣洪崙鑺ユ�?
                alpha_x_PML_optimized(i,j) = K_main_damping_x(i,j).*(alpha_x_PML_basic(i,j) + alpha_x_PML_exp(i,j));
                alpha_z_PML_optimized(i,j) = K_main_damping_z(i,j).*(alpha_z_PML_basic(i,j) + alpha_z_PML_exp(i,j));
    
                %%%%%% CPML %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                b_x_CPML(i,j) = exp(-(alpha_x_PML(i,j)./beta_x(i,j) + eta_x(i,j)).*dt);
                b_z_CPML(i,j) = exp(-(alpha_z_PML(i,j)./beta_z(i,j) + eta_z(i,j)).*dt);
                a_x_CPML(i,j) = (alpha_x_PML(i,j)/(beta_x(i,j)*(alpha_x_PML(i,j) + beta_x(i,j)*eta_x(i,j)))).*(beta_x(i,j) - 1);
                a_z_CPML(i,j) = (alpha_z_PML(i,j)/(beta_z(i,j)*(alpha_z_PML(i,j) + beta_z(i,j)*eta_z(i,j)))).*(beta_z(i,j) - 1);  
    
        %%%%%% 閸栧搫鐓?2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif i <= PML_l && j > PML_l && j < nx_PML - PML_l + 1                    %   �ϱ߽�

                z_PML = dz*(PML_l - i);     % X_PML = dX*(j - (nX_PML - PML_l) - 1);
                x_PML = 0; 

                eta_x(i,j) = 0;
                eta_z(i,j) = eta_0*pi*f0*(1 - (z_PML/PML_z).^P_eta);
%                 eta_X(i,j) = P_stable.*eta_Z(i,j);
                
                beta_x(i,j) = 1;    % beta_0閸欐�???1-10
                beta_z(i,j) = 1 + (beta_0_z - 1).*(z_PML/PML_z).^P_beta;
%                 beta_X(i,j) = P_stable.*beta_Z(i,j);                
    
    %             K_main_damping_X(i,j) = log(1./Reflection_coefficient).*((NPOWER+1).*Vmin./(2*PML_X));
                K_main_damping_z(i,j) = log(1./reflection_coefficient).*((n_power+1).*Vmin(i,j)./(2*PML_z));
 
                if attenuation_type == "exp"                
                    alpha_z_PML_basic(i,j) = (z_PML/PML_z).^n_power;
%                     alpha_z_PML_basic(i,j) = 0;

                elseif attenuation_type == "trigonometric"
                    alpha_z_PML_basic(i,j) = 1 - cos(pi*(z_PML)/(2*PML_z));  
%                     alpha_z_PML_basic(i,j) = 0;
                end
                alpha_z_PML_exp(i,j) = gamma_damping.*exp(-(delta_damping.*PML_z./z_PML));

%                 alpha_z_PML_exp(i,j) = 0;

%                 alpha_x_PML_basic(i,j) = (X_PML/PML_X).^P_eta;
%                 alpha_X_PML_exp(i,j) = gamma_damping.*exp(-(delta_damping.*PML_X./X_PML));
    
                %   缂佸繐鍚?鐞涙澘鍣洪崙鑺ユ�? 
                alpha_z_PML(i,j) = K_main_damping_z(i,j).*alpha_z_PML_basic(i,j);
    %             alpha_X_PML(i,j) = K_main_damping_X(i,j).*alpha_X_PML_basic(i,j);
                alpha_x_PML(i,j) = P_stable.*alpha_z_PML(i,j);

                %   娴兼ê瀵叉潻鍥╂畱鐞涙澘鍣洪崙鑺ユ�?
                alpha_z_PML_optimized(i,j) = K_main_damping_z(i,j).*(alpha_z_PML_basic(i,j) + alpha_z_PML_exp(i,j));
    %             alpha_X_PML_optimized(i,j) = K_main_damping_X(i,j).*(alpha_X_PML_basic(i,j) + alpha_X_PML_exp(i,j));
                alpha_x_PML_optimized(i,j) = P_stable.*alpha_z_PML_optimized(i,j);
    
                % if media == "two layers"            %%%%    ���ɱ���߽�
                %     alpha_x_PML_optimized(i,j) = 0;
                %     alpha_z_PML_optimized(i,j) = 0;
                % end

                %%%%%% CPML %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                b_x_CPML(i,j) = exp(-(alpha_x_PML(i,j)./beta_x(i,j) + eta_x(i,j)).*dt);
                b_z_CPML(i,j) = exp(-(alpha_z_PML(i,j)./beta_z(i,j) + eta_z(i,j)).*dt);
                a_x_CPML(i,j) = (alpha_x_PML(i,j)/(beta_x(i,j)*(alpha_x_PML(i,j) + beta_x(i,j)*eta_x(i,j)))).*(beta_x(i,j) - 1);
                a_z_CPML(i,j) = (alpha_z_PML(i,j)/(beta_z(i,j)*(alpha_z_PML(i,j) + beta_z(i,j)*eta_z(i,j)))).*(beta_z(i,j) - 1);              
    
            elseif i > nz_PML - PML_l && i <= nz_PML && j > PML_l && j <= nx_PML - PML_l       %    
%                 X_PML = 0;      
                
                z_PML = dz*(i - (nz_PML - PML_l) - 1);      x_PML = 0;
    
                eta_x(i,j) = 0;                
                eta_z(i,j) = eta_0*pi*f0*(1 - (z_PML/PML_z).^P_eta);
%                 eta_X(i,j) = P_stable.*eta_Z(i,j);
                
                
                beta_x(i,j) = 1;    % beta_0閸欐�???1-10
%                 beta_X(i,j) = 1 + (beta_0_X - 1).*(X_PML/PML_X).^P_beta;
                beta_z(i,j) = 1 + (beta_0_z - 1).*(z_PML/PML_z).^P_beta;
%                 beta_X(i,j) = P_stable.*beta_Z(i,j);
    
    %             K_main_damping_X(i,j) = log(1./Reflection_coefficient).*((NPOWER+1).*Vmin./(2*PML_X));
                K_main_damping_z(i,j) = log(1./reflection_coefficient).*((n_power+1).*Vmin(i,j)./(2*PML_z));
    
                if attenuation_type == "exp"                
                    alpha_z_PML_basic(i,j) = (z_PML/PML_z).^n_power;
                elseif attenuation_type == "trigonometric"
                    alpha_z_PML_basic(i,j) = 1 - cos(pi*(z_PML)/(2*PML_z));  
                end
                alpha_z_PML_exp(i,j) = gamma_damping.*exp(-(delta_damping.*PML_z./z_PML));
%                 alpha_X_PML_basic(i,j) = (X_PML/PML_X).^P_eta;
%                 alpha_X_PML_exp(i,j) = gamma_damping.*exp(-(delta_damping.*PML_X./X_PML));
    
                %   缂佸繐鍚?鐞涙澘鍣洪崙鑺ユ�? 
                alpha_z_PML(i,j) = K_main_damping_z(i,j).*alpha_z_PML_basic(i,j);
    %             alpha_X_PML(i,j) = K_main_damping_X(i,j).*alpha_X_PML_basic(i,j);
                alpha_x_PML(i,j) = P_stable.*alpha_z_PML(i,j);
    
                %   娴兼ê瀵叉潻鍥╂畱鐞涙澘鍣洪崙鑺ユ�?
                alpha_z_PML_optimized(i,j) = K_main_damping_z(i,j).*(alpha_z_PML_basic(i,j) + alpha_z_PML_exp(i,j));
    %             alpha_X_PML_optimized(i,j) = K_main_damping_X(i,j).*(alpha_X_PML_basic(i,j) + alpha_X_PML_exp(i,j));
                alpha_x_PML_optimized(i,j) = P_stable.*alpha_z_PML_optimized(i,j);
    
                %%%%%% CPML %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                b_x_CPML(i,j) = exp(-(alpha_x_PML(i,j)./beta_x(i,j) + eta_x(i,j)).*dt);
                b_z_CPML(i,j) = exp(-(alpha_z_PML(i,j)./beta_z(i,j) + eta_z(i,j)).*dt);
                a_x_CPML(i,j) = (alpha_x_PML(i,j)/(beta_x(i,j)*(alpha_x_PML(i,j) + beta_x(i,j)*eta_x(i,j)))).*(beta_x(i,j) - 1);
                a_z_CPML(i,j) = (alpha_z_PML(i,j)/(beta_z(i,j)*(alpha_z_PML(i,j) + beta_z(i,j)*eta_z(i,j)))).*(beta_z(i,j) - 1);     
    
        %%%%%%%%% 閸栧搫鐓?3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif i > PML_l && i <= nz_PML - PML_l && j <= PML_l      %    ����м�
    
                z_PML = 0;          x_PML = dx*(PML_l - j);      
    
                eta_x(i,j) = eta_0*pi*f0*(1 - (x_PML/PML_x).^P_eta);
                eta_z(i,j) = 0;
%                 eta_Z(i,j) = P_stable.*eta_X(i,j);

                beta_x(i,j) = 1 + (beta_0_x - 1).*(x_PML/PML_x).^P_beta;    % beta_0閸欐�???1-10
                beta_z(i,j) = 1;
%                 beta_Z(i,j) = P_stable.*beta_X(i,j);


    %             K_main_damping_X(i,j) = log(1./Reflection_coefficient).*((NPOWER+1).*Vmin./(2*PML_X));
                K_main_damping_x(i,j) = log(1./reflection_coefficient).*((n_power+1).*Vmin(i,j)./(2*PML_x));
    
                if attenuation_type == "exp"
                    alpha_x_PML_basic(i,j) = (x_PML/PML_x).^n_power;
                elseif attenuation_type == "trigonometric"
                    alpha_x_PML_basic(i,j) = 1 - cos(pi*(x_PML)/(2*PML_x));    %    余弦�?
                end
                alpha_x_PML_exp(i,j) = gamma_damping.*exp(-(delta_damping.*PML_x./x_PML));
    %             alpha_X_PML_basic(i,j) = (X_PML/PML_X).^P_eta;
    %             alpha_X_PML_exp(i,j) = gamma_damping.*exp(-(delta_damping.*PML_X./X_PML));
    
                %   缂佸繐鍚?鐞涙澘鍣洪崙鑺ユ�? 
                alpha_x_PML(i,j) = K_main_damping_x(i,j).*alpha_x_PML_basic(i,j);
    %             alpha_X_PML(i,j) = K_main_damping_X(i,j).*alpha_X_PML_basic(i,j);
                alpha_z_PML(i,j) = P_stable.*alpha_x_PML(i,j);
    
                %   娴兼ê瀵叉潻鍥╂畱鐞涙澘鍣洪崙鑺ユ�?
                alpha_x_PML_optimized(i,j) = K_main_damping_x(i,j).*(alpha_x_PML_basic(i,j) + alpha_x_PML_exp(i,j));
    %             alpha_X_PML_optimized(i,j) = K_main_damping_X(i,j).*(alpha_X_PML_basic(i,j) + alpha_X_PML_exp(i,j));
                alpha_z_PML_optimized(i,j) = P_stable.*alpha_x_PML_optimized(i,j);
    
                %%%%%% CPML %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                b_x_CPML(i,j) = exp(-(alpha_x_PML(i,j)./beta_x(i,j) + eta_x(i,j)).*dt);
                b_z_CPML(i,j) = exp(-(alpha_z_PML(i,j)./beta_z(i,j) + eta_z(i,j)).*dt);
                a_x_CPML(i,j) = (alpha_x_PML(i,j)/(beta_x(i,j)*(alpha_x_PML(i,j) + beta_x(i,j)*eta_x(i,j)))).*(beta_x(i,j) - 1);
                a_z_CPML(i,j) = (alpha_z_PML(i,j)/(beta_z(i,j)*(alpha_z_PML(i,j) + beta_z(i,j)*eta_z(i,j)))).*(beta_z(i,j) - 1);   
    
            elseif i > PML_l && i <= nz_PML - PML_l && j > nx_PML - PML_l && j <= nx_PML
    
                z_PML = 0;              x_PML = dx*(j - (nx_PML - PML_l) - 1);      
    
                eta_x(i,j) = eta_0*pi*f0*(1 - (x_PML/PML_x).^P_eta);
                eta_z(i,j) = 0;
%                 eta_Z(i,j) = P_stable.*eta_X(i,j);

                beta_x(i,j) = 1 + (beta_0_x - 1).*(x_PML/PML_x).^P_beta;    % beta_0閸欐�???1-10
%                 beta_Z(i,j) = P_stable.*beta_X(i,j);
                beta_z(i,j) = 1;
    
    %             K_main_damping_X(i,j) = log(1./Reflection_coefficient).*((NPOWER+1).*Vmin./(2*PML_X));
                K_main_damping_x(i,j) = log(1./reflection_coefficient).*((n_power+1).*Vmin(i,j)./(2*PML_x));
    
                if attenuation_type == "exp"
                    alpha_x_PML_basic(i,j) = (x_PML/PML_x).^n_power;
                elseif attenuation_type == "trigonometric"
                    alpha_x_PML_basic(i,j) = 1 - cos(pi*(x_PML)/(2*PML_x));    %    余弦�?
                end
                alpha_x_PML_exp(i,j) = gamma_damping.*exp(-(delta_damping.*PML_x./x_PML));
    %             alpha_X_PML_basic(i,j) = (X_PML/PML_X).^P_eta;
    %             alpha_X_PML_exp(i,j) = gamma_damping.*exp(-(delta_damping.*PML_X./X_PML));
    
                %   缂佸繐鍚?鐞涙澘鍣洪崙鑺ユ�? 
                alpha_x_PML(i,j) = K_main_damping_x(i,j).*alpha_x_PML_basic(i,j);
    %             alpha_X_PML(i,j) = K_main_damping_X(i,j).*alpha_X_PML_basic(i,j);
                alpha_z_PML(i,j) = P_stable.*alpha_x_PML(i,j);
    
                %   娴兼ê瀵叉潻鍥╂畱鐞涙澘鍣洪崙鑺ユ�?
                alpha_x_PML_optimized(i,j) = K_main_damping_x(i,j).*(alpha_x_PML_basic(i,j) + alpha_x_PML_exp(i,j));
    %             alpha_X_PML_optimized(i,j) = K_main_damping_X(i,j).*(alpha_X_PML_basic(i,j) + alpha_X_PML_exp(i,j));
                alpha_z_PML_optimized(i,j) = P_stable.*alpha_x_PML_optimized(i,j);
    
                %%%%%% CPML %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                b_x_CPML(i,j) = exp(-(alpha_x_PML(i,j)./beta_x(i,j) + eta_x(i,j)).*dt);
                b_z_CPML(i,j) = exp(-(alpha_z_PML(i,j)./beta_z(i,j) + eta_z(i,j)).*dt);
                a_x_CPML(i,j) = (alpha_x_PML(i,j)/(beta_x(i,j)*(alpha_x_PML(i,j) + beta_x(i,j)*eta_x(i,j)))).*(beta_x(i,j) - 1);
                a_z_CPML(i,j) = (alpha_z_PML(i,j)/(beta_z(i,j)*(alpha_z_PML(i,j) + beta_z(i,j)*eta_z(i,j)))).*(beta_z(i,j) - 1);               
            end
        end
    end

end