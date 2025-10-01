function Qs = calculate_volumetric_sed_flux_waves(H0,L0,rho_s,rho,g,D50,nu,slope,alpha0)

    w = calc_settling_velocity(rho_s,rho,g,D50,nu);
    R = rho_s/rho - 1;
    w_star = w/sqrt(g*D50);
    theta_0 = 0.1*((H0/D50)^2.3)*(sqrt(H0/L0))*exp(-6.1*w_star);
    

   for a = 1:numel(alpha0)
    if alpha0(a) < pi/2 && alpha0(a) > -pi/2
       
        alpha(a) = alpha0(a)./(pi/2);
        if alpha(a) < 0 
            alpha(a) = -alpha(a);
        end
        
        theta = (sin(2.*alpha0(a).*(1-(0.4.*(alpha(a)).*(1-(alpha(a)))))));        


        s = rho_s/rho;
        Qs(a) = theta*(H0*sqrt(slope)*sqrt(s-1)*g*D50);

    else

        Qs(a) = NaN;
   
    end
    % if neg_side
    %     Qs = -Qs;
    % end
   end


    

    function settling_velocity = calc_settling_velocity(rho_s,rho,g,D50,nu)
        R = (rho_s/rho) - 1;
        D_star = (R*g.*(D50.^3))./(nu^2);
        log_D_star = log10(D_star);  
        log_W_star = -3.76715 + ...
                         1.92944 .* log_D_star - ...
                         0.09815 .* (log_D_star.^2) - ...
                         0.00575 .* (log_D_star.^3) + ...
                         0.00056 .* (log_D_star.^4);
        W_star = 10.^(log_W_star);
        settling_velocity = ((W_star.*R.*g.*nu)).^(1/3);

       
    end

end