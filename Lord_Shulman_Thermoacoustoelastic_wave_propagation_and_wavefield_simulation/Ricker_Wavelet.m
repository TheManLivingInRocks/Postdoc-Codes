function wavelet = Ricker_Wavelet(dt,fm,t0,nt)

    %   dt为采样率
    %   fm为峰值频率
    %   nt为采样点数
    t_min = -dt*round(nt/2);
    t_max = -t_min - dt;

    % t_wavelet = t_min:dt:t_max;

    t_wavelet = (-round(nt/2) + 1)*dt:dt:round(nt/2)*dt;

    wavelet = (1 - 2.*(t_wavelet - t0).^2*pi^2*fm^2).*exp(-(t_wavelet - t0).^2*pi^2*fm^2);

end