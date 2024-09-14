function [Vp_iso_N,Vs_N,Rho_N,K_N,Alpha_N,T_N,cV_N,Beta_N,Lambda_iso_N,Mu_N] = Add_PML_layers_CTE(nz,nx,PML_l,Vp_iso_N,Vs_N,Rho_N,K_N,Alpha_N,T_N,cV_N,Beta_N,Lambda_iso_N,Mu_N)

%   给预定义数组"镶边"，以包含PML边界层
%   Vp_iso_N
    tmp = zeros(nz + 2*PML_l,nx + 2*PML_l);
    tmp(PML_l + 1:nz + PML_l,PML_l + 1:nx + PML_l) = Vp_iso_N;
    tmp(1:PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*Vp_iso_N(1,:);
    tmp(nz + PML_l + 1:nz + 2*PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*Vp_iso_N(nz,:);
    tmp(:,1:PML_l) = tmp(:,PML_l + 1)*ones(1,PML_l);
    tmp(:,nx + PML_l + 1:nx + 2*PML_l) = tmp(:,nx + PML_l)*ones(1,PML_l);
    Vp_iso_N = tmp;

%   Vs_N
    tmp = zeros(nz + 2*PML_l,nx + 2*PML_l);
    tmp(PML_l + 1:nz + PML_l,PML_l + 1:nx + PML_l) = Vs_N;
    tmp(1:PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*Vs_N(1,:);
    tmp(nz + PML_l + 1:nz + 2*PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*Vs_N(nz,:);
    tmp(:,1:PML_l) = tmp(:,PML_l + 1)*ones(1,PML_l);
    tmp(:,nx + PML_l + 1:nx + 2*PML_l) = tmp(:,nx + PML_l)*ones(1,PML_l);
    Vs_N = tmp;

%   Rho_N
    tmp = zeros(nz + 2*PML_l,nx + 2*PML_l);
    tmp(PML_l + 1:nz + PML_l,PML_l + 1:nx + PML_l) = Rho_N;
    tmp(1:PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*Rho_N(1,:);
    tmp(nz + PML_l + 1:nz + 2*PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*Rho_N(nz,:);
    tmp(:,1:PML_l) = tmp(:,PML_l + 1)*ones(1,PML_l);
    tmp(:,nx + PML_l + 1:nx + 2*PML_l) = tmp(:,nx + PML_l)*ones(1,PML_l);
    Rho_N = tmp;

%   K_N
    tmp = zeros(nz + 2*PML_l,nx + 2*PML_l);
    tmp(PML_l + 1:nz + PML_l,PML_l + 1:nx + PML_l) = K_N;
    tmp(1:PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*K_N(1,:);
    tmp(nz + PML_l + 1:nz + 2*PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*K_N(nz,:);
    tmp(:,1:PML_l) = tmp(:,PML_l + 1)*ones(1,PML_l);
    tmp(:,nx + PML_l + 1:nx + 2*PML_l) = tmp(:,nx + PML_l)*ones(1,PML_l);
    K_N = tmp;

%   Alpha_N
    tmp = zeros(nz + 2*PML_l,nx + 2*PML_l);
    tmp(PML_l + 1:nz + PML_l,PML_l + 1:nx + PML_l) = Alpha_N;
    tmp(1:PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*Alpha_N(1,:);
    tmp(nz + PML_l + 1:nz + 2*PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*Alpha_N(nz,:);
    tmp(:,1:PML_l) = tmp(:,PML_l + 1)*ones(1,PML_l);
    tmp(:,nx + PML_l + 1:nx + 2*PML_l) = tmp(:,nx + PML_l)*ones(1,PML_l);
    Alpha_N = tmp;

%   T_N
    tmp = zeros(nz + 2*PML_l,nx + 2*PML_l);
    tmp(PML_l + 1:nz + PML_l,PML_l + 1:nx + PML_l) = T_N;
    tmp(1:PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*T_N(1,:);
    tmp(nz + PML_l + 1:nz + 2*PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*T_N(nz,:);
    tmp(:,1:PML_l) = tmp(:,PML_l + 1)*ones(1,PML_l);
    tmp(:,nx + PML_l + 1:nx + 2*PML_l) = tmp(:,nx + PML_l)*ones(1,PML_l);
    T_N = tmp;

%   cV_N
    tmp = zeros(nz + 2*PML_l,nx + 2*PML_l);
    tmp(PML_l + 1:nz + PML_l,PML_l + 1:nx + PML_l) = cV_N;
    tmp(1:PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*cV_N(1,:);
    tmp(nz + PML_l + 1:nz + 2*PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*cV_N(nz,:);
    tmp(:,1:PML_l) = tmp(:,PML_l + 1)*ones(1,PML_l);
    tmp(:,nx + PML_l + 1:nx + 2*PML_l) = tmp(:,nx + PML_l)*ones(1,PML_l);
    cV_N = tmp;

%   Beta_N
    tmp = zeros(nz + 2*PML_l,nx + 2*PML_l);
    tmp(PML_l + 1:nz + PML_l,PML_l + 1:nx + PML_l) = Beta_N;
    tmp(1:PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*Beta_N(1,:);
    tmp(nz + PML_l + 1:nz + 2*PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*Beta_N(nz,:);
    tmp(:,1:PML_l) = tmp(:,PML_l + 1)*ones(1,PML_l);
    tmp(:,nx + PML_l + 1:nx + 2*PML_l) = tmp(:,nx + PML_l)*ones(1,PML_l);
    Beta_N = tmp;

%   Lambda_iso_N
    tmp = zeros(nz + 2*PML_l,nx + 2*PML_l);
    tmp(PML_l + 1:nz + PML_l,PML_l + 1:nx + PML_l) = Lambda_iso_N;
    tmp(1:PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*Lambda_iso_N(1,:);
    tmp(nz + PML_l + 1:nz + 2*PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*Lambda_iso_N(nz,:);
    tmp(:,1:PML_l) = tmp(:,PML_l + 1)*ones(1,PML_l);
    tmp(:,nx + PML_l + 1:nx + 2*PML_l) = tmp(:,nx + PML_l)*ones(1,PML_l);
    Lambda_iso_N = tmp;

%   Mu_N
    tmp = zeros(nz + 2*PML_l,nx + 2*PML_l);
    tmp(PML_l + 1:nz + PML_l,PML_l + 1:nx + PML_l) = Mu_N;
    tmp(1:PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*Mu_N(1,:);
    tmp(nz + PML_l + 1:nz + 2*PML_l,PML_l + 1:nx + PML_l) = ones(PML_l,1)*Mu_N(nz,:);
    tmp(:,1:PML_l) = tmp(:,PML_l + 1)*ones(1,PML_l);
    tmp(:,nx + PML_l + 1:nx + 2*PML_l) = tmp(:,nx + PML_l)*ones(1,PML_l);
    Mu_N = tmp;


end