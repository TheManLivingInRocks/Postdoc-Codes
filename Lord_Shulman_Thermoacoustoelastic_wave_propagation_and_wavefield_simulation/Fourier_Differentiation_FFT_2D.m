function Differentiation_m_th_function = Fourier_Differentiation_FFT_2D(f,dx,dz,dkx,dkz,grid_mode,m,derivative,direction)


%   基于Reddy and Weideman(2000)的DMsuite中的fourdifft.m函数。
%   f(x) 使用傅里叶微分过程。假设函数为 2pi 周期函数，
%   输入数据值 f 应对应于[0,2pi]上 N 个等距点的函数采样点。使用了FFT。

%   输入：
%   f :      f(x) 的采样点向量，在x = 2pi/N, 4pi/N, ... , (N-1)2pi/N
%   m :     所需的导数（非负整数）

%   输出：
%   Differentiation_m_th_function:       f 的第m导数

%  S.C. Reddy, J.A.C. Weideman 2000. Corrected for MATLAB R13 
%  by JACW, April 2003.

    PI  = 3.1415926;
%     wx  = 2*PI*dkx;
%     wz  = 2*PI*dkz;

%     x=2*pi*(0:N-1)'/N;                       % gridpoints
%     h=2*pi/N;                                % grid spacing
%     zi=sqrt(-1);
%     kk=(1:N-1)';
%     n1=floor((N-1)/2); n2=ceil((N-1)/2);



if grid_mode == "NG"

    if derivative == "x"
        %   计算x方向的导数
        %   输入 f(nz,nx)

%     f = f(:);       %   确保f是一个列向量
%     N = length(f);              %   采样点数量

        f = f';     %   确保x方向是列向量
        [nx,nz] = size(f);          %   x方向和z方向采样点数量

%     N1 = floor((N-1)/2);        %   设置波数，分为两部分
%     N2 = (-N/2)*rem(m+1,2)*ones(rem(N+1,2));
%     k = [(0:N1) N2 (-N1:-1)]';

        nx1 = floor((nx-1)/2);
        nx2 = -(nx/2)*rem(m+1,2)*ones(rem(nx+1,2));
        kx = (2*PI*dkx).*[(0:nx1) nx2 (-nx1:-1)]';

        Differentiation_m_th_function = ifft2(((1i*kx).^m).*fft2(f));      %   先变换到傅里叶空间，取导数，再返回物理空间

        Differentiation_m_th_function = Differentiation_m_th_function';


        if max(abs(imag(f))) == 0
            Differentiation_m_th_function = real(Differentiation_m_th_function);        %   读入数据，输出实导数
        end
    elseif derivative == "z"
        [nz,nx] = size(f);

        nz1 = floor((nz-1)/2);
        nz2 = -(nz/2)*rem(m+1,2)*ones(rem(nz+1,2));
        kz = (2*PI*dkz).*[(0:nz1) nz2 (-nz1:-1)]';

        Differentiation_m_th_function = ifft2(((1i*kz).^m).*fft2(f));      %   先变换到傅里叶空间，取导数，再返回物理空间

%         Differentiation_m_th_function = Differentiation_m_th_function';

        if max(abs(imag(f))) == 0
            Differentiation_m_th_function = real(Differentiation_m_th_function);        %   读入数据，输出实导数
        end
    end

elseif grid_mode == "SG"

    if derivative == "x"
        %   计算x方向的导数
        %   输入 f(nz,nx)

%     f = f(:);       %   确保f是一个列向量
%     N = length(f);              %   采样点数量

        f = f';     %   确保x方向是列向量
        [nx,nz] = size(f);          %   x方向和z方向采样点数量

%     N1 = floor((N-1)/2);        %   设置波数，分为两部分
%     N2 = (-N/2)*rem(m+1,2)*ones(rem(N+1,2));
%     k = [(0:N1) N2 (-N1:-1)]';

        nx1 = floor((nx-1)/2);
        nx2 = -(nx/2)*rem(m+1,2)*ones(rem(nx+1,2));
        kx = (2*PI*dkx).*[(0:nx1) nx2 (-nx1:-1)]';
        if direction == 1
            Differentiation_m_th_function = ifft2(((1i*kx).^m.*exp((1i*kx).*dx/2)).*fft2(f));      %   先变换到傅里叶空间，取导数，再返回物理空间
        elseif direction == -1
            Differentiation_m_th_function = ifft2(((1i*kx).^m.*exp(-(1i*kx).*dx/2)).*fft2(f));
        end

        Differentiation_m_th_function = Differentiation_m_th_function';
        if max(abs(imag(f))) == 0
            Differentiation_m_th_function = real(Differentiation_m_th_function);        %   读入数据，输出实导数
        end
    elseif derivative == "z"
        [nz,nx] = size(f);

        nz1 = floor((nz-1)/2);
        nz2 = -(nz/2)*rem(m+1,2)*ones(rem(nz+1,2));
        kz = (2*PI*dkz).*[(0:nz1) nz2 (-nz1:-1)]';

        if direction == 1
            Differentiation_m_th_function = ifft2(((1i*kz).^m.*exp((1i*kz).*dz/2)).*fft2(f));      %   先变换到傅里叶空间，取导数，再返回物理空间
        elseif direction == -1
            Differentiation_m_th_function = ifft2(((1i*kz).^m.*exp(-(1i*kz).*dz/2)).*fft2(f));      %   先变换到傅里叶空间，取导数，再返回物理空间
        end
%         Differentiation_m_th_function = Differentiation_m_th_function';

        if max(abs(imag(f))) == 0
            Differentiation_m_th_function = real(Differentiation_m_th_function);        %   读入数据，输出实导数
        end
    end

elseif grid_mode == "RSG"

    if derivative == "x"

        f = f';     %   确保x方向是列向量
        [nx,nz] = size(f);          %   x方向和z方向采样点数量

        nx1 = floor((nx-1)/2);
        nx2 = -(nx/2)*rem(m+1,2)*ones(rem(nx+1,2));
        kx = (2*PI*dkx).*[(0:nx1) nx2 (-nx1:-1)]';

        nz1 = floor((nz-1)/2);
        nz2 = -(nz/2)*rem(m+1,2)*ones(rem(nz+1,2));
        kz = (2*PI*dkz).*[(0:nz1) nz2 (-nz1:-1)];

        if direction == 1
            Differentiation_m_th_function = ifft2(((1i*kx).^m.*exp((1i*kx).*(dx/2)+(1i*kz).*(dz/2))).*fft2(f));      %   先变换到傅里叶空间，取导数，再返回物理空间
        elseif direction == -1
            Differentiation_m_th_function = ifft2(((1i*kx).^m.*exp(-(1i*kx).*(dx/2)-(1i*kz).*(dz/2))).*fft2(f));
        end
        Differentiation_m_th_function = Differentiation_m_th_function';

%         Differentiation_m_th_function = Differentiation_m_th_function';
        if max(abs(imag(f))) == 0
            Differentiation_m_th_function = real(Differentiation_m_th_function);        %   读入数据，输出实导数
        end
    elseif derivative == "z"
        [nz,nx] = size(f);

        nz1 = floor((nz-1)/2);
        nz2 = -(nz/2)*rem(m+1,2)*ones(rem(nz+1,2));
        kz = (2*PI*dkz).*[(0:nz1) nz2 (-nz1:-1)]';

        nx1 = floor((nx-1)/2);
        nx2 = -(nx/2)*rem(m+1,2)*ones(rem(nx+1,2));
        kx = (2*PI*dkx).*[(0:nx1) nx2 (-nx1:-1)];


        if direction == 1
            Differentiation_m_th_function = ifft2(((1i*kz).^m.*exp((1i*kx).*(dx/2)+(1i*kz).*(dz/2))).*fft2(f));      %   先变换到傅里叶空间，取导数，再返回物理空间
        elseif direction == -1
            Differentiation_m_th_function = ifft2(((1i*kz).^m.*exp(-(1i*kx).*(dx/2)-(1i*kz).*(dz/2))).*fft2(f));      %   先变换到傅里叶空间，取导数，再返回物理空间
        end
%         Differentiation_m_th_function = Differentiation_m_th_function';

        if max(abs(imag(f))) == 0
            Differentiation_m_th_function = real(Differentiation_m_th_function);        %   读入数据，输出实导数
        end


    end


end