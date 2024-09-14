function OUTPUT = Cerjan_Absorbing_Boundary(a,nX_PML,nZ_PML,width)

%     a = 0;
%     a = 0.00016;
    % a = 0.016;

%     a = 0.08;


    OUTPUT = ones(nZ_PML,nX_PML);

    for jj=1:1:nZ_PML
        for ii=1:1:nX_PML
            CAB = 0;
            if ii <= width   
                % CAB = a*((width-1)-ii)^2;
                CAB = a*(width - ii)^2;
            end
            if ii >= nX_PML-width         
                % CAB = a*(ii-(nX_PML-(width-1)))^2;
                CAB = a*(ii - (nX_PML - width))^2;
            end
            OUTPUT(jj,ii) = (exp(-CAB))*OUTPUT(jj,ii);
        end
    end
      
    
    for jj=1:1:nZ_PML
        for ii=1:1:nX_PML
            CAB = 0;
            if jj <= width
                % CAB = a*((width-1)-jj)^2;
                CAB = a*(width - jj)^2;
            end
            if jj >= nZ_PML - width
                % CAB = a*(jj-(nZ_PML-(width-1)))^2;
                CAB = a*(jj - (nZ_PML - width))^2;
            end
            OUTPUT(jj,ii) = (exp(-CAB))*OUTPUT(jj,ii);
        end
    end
end