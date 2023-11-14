
clc
clear
close all
% OBJECTIVE: Find periods of relative quiescence in Lake Superior Wind/Wave Data to compare to UMWM-Titan

makeplots = 1;
filen = 'BuoyData/Superior/45004h2019'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,varu,varht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,dir_t,g_t);

if ~isnan(qi)
    stru = "Avg |u| is " + avgu + " with std of " + varu + " m/s";
    strht = "Avg wave height is " + avght + " with std of " + varht + " m";
    
    
    if makeplots
        timeplot = 1:numel(umag);
        
        figure;
        boxplot([umag gust],'Notch','on','Labels',{'|u|','gust'},'Whisker',1)
        title(yr)
        
        figure;
        boxplot(dir,'Notch','on','Labels',{'direction'},'Whisker',1)
        title(yr)
        
        figure;
        boxplot([umag(qi) gust(qi)],'Notch','on','Labels',{'|u|','gust'},'Whisker',1)
        title(strcat('ROI',' ',yr))
        
        figure;
        boxplot(dir(qi),'Notch','on','Labels',{'direction'},'Whisker',1)
        title(strcat('ROI',' ',yr))
        
        figure;
        subplot(2,1,1)
        plot(timeplot,umag,'-k')
        grid on;
        hold on
        patch([timeplot(qi(1)) timeplot(qi(end)) timeplot(qi(end)) timeplot(qi(1))],[0 0 20 20],'red','FaceAlpha',0.5)
        ylim([0 20])
        xlabel('time (# measurement pts)')
        ylabel('u [m/s]')
        title(yr)
        subplot(2,1,2)
        yyaxis right
        h1 = plot(qi,umag(qi),'-k');
        hold on
        h2 = plot(qi,gust(qi),'--r');
        yyaxis left
        h3 = plot(qi,deg2rad(dir(qi)),'-g');
        ylim([0 2*pi])
        grid on;
        title('ROI')
        xlabel('time (# measurement pts)')
        
        
        figure;
        plot(qi,umag(qi),'ok','MarkerFaceColor','k')
        hold on
        plot(qi,waveht(qi),'-or','MarkerFaceColor','r')
        grid on;
        title({yr,stru,strht})
    else
        disp(stru);
        disp(strht);
    end

end

    
