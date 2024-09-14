function [C11_iso_N,C12_iso_N,C13_iso_N,C14_iso_N,C15_iso_N,C16_iso_N,...
          C21_iso_N,C22_iso_N,C23_iso_N,C24_iso_N,C25_iso_N,C26_iso_N,...
          C31_iso_N,C32_iso_N,C33_iso_N,C34_iso_N,C35_iso_N,C36_iso_N,...
          C41_iso_N,C42_iso_N,C43_iso_N,C44_iso_N,C45_iso_N,C46_iso_N,...
          C51_iso_N,C52_iso_N,C53_iso_N,C54_iso_N,C55_iso_N,C56_iso_N,...
          C61_iso_N,C62_iso_N,C63_iso_N,C64_iso_N,C65_iso_N,C66_iso_N] = TEI_isothermal_isotropic_elastic_coefficients(nz_PML,nx_PML,PML_l,...
                                                                                                                       Vp_iso_N,Vs_N,Rho_N,...
                                                                                                                       K_N,Alpha_N,T_N,cV_N,Beta_N,...
                                                                                                                       Lambda_iso_N,Mu_N,A,B,C,...
                                                                                                                       thermoelasticity,initial_stress,P)



%   自然状态下的各向同性等温弹性系数
C11_iso_N = Lambda_iso_N + 2*Mu_N;          C12_iso_N = Lambda_iso_N;           C13_iso_N = C12_iso_N;                  C14_iso_N = zeros(nz_PML,nx_PML);    C15_iso_N = zeros(nz_PML,nx_PML);      C16_iso_N = zeros(nz_PML,nx_PML);
C21_iso_N = C12_iso_N;                                  C22_iso_N = C11_iso_N;                  C23_iso_N = Lambda_iso_N;          C24_iso_N = zeros(nz_PML,nx_PML);    C25_iso_N = zeros(nz_PML,nx_PML);       C26_iso_N = zeros(nz_PML,nx_PML);
C31_iso_N = C13_iso_N;                                  C32_iso_N = C23_iso_N;                  C33_iso_N = C11_iso_N;                  C34_iso_N = zeros(nz_PML,nx_PML);    C35_iso_N = zeros(nz_PML,nx_PML);      C36_iso_N = zeros(nz_PML,nx_PML);
C41_iso_N = C14_iso_N;                                  C42_iso_N = C24_iso_N;                  C43_iso_N = C34_iso_N;                  C44_iso_N = Mu_N;                               C45_iso_N = zeros(nz_PML,nx_PML);       C46_iso_N = zeros(nz_PML,nx_PML);
C51_iso_N = C15_iso_N;                                  C52_iso_N = C25_iso_N;                  C53_iso_N = C35_iso_N;                  C54_iso_N = C45_iso_N;                          C55_iso_N = C44_iso_N;                          C56_iso_N = zeros(nz_PML,nx_PML);
C61_iso_N = C16_iso_N;                                  C62_iso_N = C26_iso_N;                  C63_iso_N = C36_iso_N;                  C64_iso_N = C46_iso_N;                          C65_iso_N = C56_iso_N;                          C66_iso_N = C44_iso_N;



if initial_stress == 1
    
    %   预定义等温三阶弹性常数
    C123_iso_N = zeros(nz_PML,nx_PML);      C111_iso_N = zeros(nz_PML,nx_PML);     C333_iso_N = zeros(nz_PML,nx_PML);
    C113_iso_N = zeros(nz_PML,nx_PML);      C133_iso_N = zeros(nz_PML,nx_PML);     C331_iso_N = zeros(nz_PML,nx_PML); 
    C311_iso_N = zeros(nz_PML,nx_PML);      C131_iso_N = zeros(nz_PML,nx_PML);     C155_iso_N = zeros(nz_PML,nx_PML); 
    C551_iso_N = zeros(nz_PML,nx_PML);      C355_iso_N = zeros(nz_PML,nx_PML);     C553_iso_N = zeros(nz_PML,nx_PML); 


    C123_iso_N(1:nz_PML,1:nx_PML) = 2*C;
    C111_iso_N(1:nz_PML,1:nx_PML) = 2*A + 6*B + 2*C;
    C333_iso_N = C111_iso_N;
    C113_iso_N(1:nz_PML,1:nx_PML) = 2*B + 2*C;
    C133_iso_N = C113_iso_N;
    C331_iso_N = C113_iso_N;
    C311_iso_N = C113_iso_N;
    C131_iso_N = C113_iso_N;

    C155_iso_N(1:nz_PML,1:nx_PML) = B + A/2;
    C551_iso_N = C155_iso_N;
    C355_iso_N = C155_iso_N;
    C553_iso_N = C355_iso_N;

end



end


















