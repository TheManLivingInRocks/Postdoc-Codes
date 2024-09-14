function [OUTPUT1,OUTPUT2,OUTPUT3,OUTPUT4] = Spatial_Derivative_RSG_PSO(field_type,Field_Variable_1,Field_Variable_2,Field_Variable_3,Field_Variable_4,dx,dz,dfx,dfz)

    
    if field_type == "scalar"

            % dtheta_stage1_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",theta_stage1_x_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
            % dtheta_stage1_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",theta_stage1_z_MCFS_NPML_next,dx,dz,dfx,dfz,"RSG",1,1);
        dField_Variable_x_dx = Fourier_Pseudospectral_Operator("x",Field_Variable_1,dx,dz,dfx,dfz,"RSG",1,1);
        dField_Variable_z_dz = Fourier_Pseudospectral_Operator("z",Field_Variable_2,dx,dz,dfx,dfz,"RSG",1,1);

        OUTPUT1 = dField_Variable_x_dx;
        OUTPUT2 = dField_Variable_z_dz;
        OUTPUT3 = 0;
        OUTPUT4 = 0;


    elseif field_type == "vector"

            % dv_x_stage1_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",v_x_stage1_x_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            % dv_z_stage1_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",v_z_stage1_z_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            % dv_x_stage1_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",v_x_stage1_z_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
            % dv_z_stage1_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",v_z_stage1_x_MCFS_NPML_next_half,dx,dz,dfx,dfz,"RSG",1,-1);
        dField_Variable_x_x_dx = Fourier_Pseudospectral_Operator("x",Field_Variable_1,dx,dz,dfx,dfz,"RSG",1,-1);
        dField_Variable_z_z_dz = Fourier_Pseudospectral_Operator("z",Field_Variable_2,dx,dz,dfx,dfz,"RSG",1,-1);
        dField_Variable_x_z_dz = Fourier_Pseudospectral_Operator("z",Field_Variable_3,dx,dz,dfx,dfz,"RSG",1,-1);
        dField_Variable_z_x_dx = Fourier_Pseudospectral_Operator("x",Field_Variable_4,dx,dz,dfx,dfz,"RSG",1,-1);

        OUTPUT1 = dField_Variable_x_x_dx;
        OUTPUT2 = dField_Variable_z_z_dz;
        OUTPUT3 = dField_Variable_x_z_dz;
        OUTPUT4 = dField_Variable_z_x_dx;

    elseif field_type == "second-order tensor"

      % dT_PK1_xx_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",T_PK1_xx_x_MCFS_NPML_present,dx,dz,dfx,dfz,"RSG",1,1);
      % dT_PK1_zz_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",T_PK1_zz_z_MCFS_NPML_present,dx,dz,dfx,dfz,"RSG",1,1);
      % dT_PK1_xz_z_MCFS_NPML_dz = Fourier_Pseudospectral_Operator("z",T_PK1_xz_z_MCFS_NPML_present,dx,dz,dfx,dfz,"RSG",1,1);
      % dT_PK1_zx_x_MCFS_NPML_dx = Fourier_Pseudospectral_Operator("x",T_PK1_zx_x_MCFS_NPML_present,dx,dz,dfx,dfz,"RSG",1,1);

        dField_Variable_xx_x_dx = Fourier_Pseudospectral_Operator("x",Field_Variable_1,dx,dz,dfx,dfz,"RSG",1,1);
        dField_Variable_zz_z_dz = Fourier_Pseudospectral_Operator("z",Field_Variable_2,dx,dz,dfx,dfz,"RSG",1,1);
        dField_Variable_xz_z_dz = Fourier_Pseudospectral_Operator("z",Field_Variable_3,dx,dz,dfx,dfz,"RSG",1,1);
        dField_Variable_zx_x_dx = Fourier_Pseudospectral_Operator("x",Field_Variable_4,dx,dz,dfx,dfz,"RSG",1,1);

        OUTPUT1 = dField_Variable_xx_x_dx;
        OUTPUT2 = dField_Variable_zz_z_dz;
        OUTPUT3 = dField_Variable_xz_z_dz;
        OUTPUT4 = dField_Variable_zx_x_dx;

    end
end


