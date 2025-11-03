function [Qs_total,Qs_time] = calc_Qs_waves(x,y,wind_mag,wind_angle_deg,rho,g)
% angle in degrees (degrees CCW from east)
    
    make_plot = 0;
    Qs_total = NaN;
    Qs_time = NaN;

    %make_lake_wind_plot(x,y,wind_angle_deg(1))

    % Step 1: calculate shoreline angle
    window_size_ang = 50; 
    shoreline_angle_rad = calc_regional_shoreline_angle(x, y,window_size_ang);
    sang = rad2deg(shoreline_angle_rad);

    if make_plot
        figure('Name','shoreline angle');
        plot(x,y,'-k','LineWidth',2)
        hold on
        scatter(x,y,50,rad2deg(shoreline_angle_rad),'filled')
        colorbar
    end

     thetas = -90:90;
     Qsmax = zeros(size(sang));
 
     % Step 2: Calculate fetches for all points 
     fetch_km = calc_fetch(x,y,wind_angle_deg);
     fetch = fetch_km.*1000;

     if make_plot
        figure('Name','Max Fetch');
        plot(x,y,'-k','LineWidth',2)
        hold on
        scatter(x,y,50,max(fetch),'filled')
        colorbar
     end


     % add calculation of E pdf here !!!!

     
    
end

 


