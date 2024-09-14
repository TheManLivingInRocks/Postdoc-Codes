function OUTPUT = Fourier_Pseudospectral_Operator(derivative,Input,dx,dz,dkx,dkz,grid_mode,order,direction)

%     [nZ_PML,nX_PML] = size(Input);

%     tmp_Input = zeros(nZ_PML,nX_PML);
%     tmp_dInput_dX = zeros(nZ_PML,nX_PML);
%     tmp_dInput_dZ = zeros(nZ_PML,nX_PML);
    if order == 1
        if derivative == "x"
    
    %         OUTPUT = Fourier_Differentiation_FFT_2D(Input,dfx,dfz,1,"x");
    
            OUTPUT = Fourier_Differentiation_FFT_2D(Input,dx,dz,dkx,dkz,grid_mode,1,"x",direction);
    
        elseif derivative == "z"
    
    %         OUTPUT = Fourier_Differentiation_FFT_2D(Input,dfx,dfz,1,"z");
    
            OUTPUT = Fourier_Differentiation_FFT_2D(Input,dx,dz,dkx,dkz,grid_mode,1,"z",direction);
    
        end
    elseif order == 2
        if derivative == "x"
    
    %         OUTPUT = Fourier_Differentiation_FFT_2D(Input,dfx,dfz,1,"x");
    
            OUTPUT = Fourier_Differentiation_FFT_2D(Input,dx,dz,dkx,dkz,grid_mode,2,"x",direction);
    
        elseif derivative == "z"
    
    %         OUTPUT = Fourier_Differentiation_FFT_2D(Input,dfx,dfz,1,"z");
    
            OUTPUT = Fourier_Differentiation_FFT_2D(Input,dx,dz,dkx,dkz,grid_mode,2,"z",direction);
    
        end        
    end


end





