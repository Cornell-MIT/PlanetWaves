
clc
clear
close all

makeplots = 0;
filen = "BuoyData/Superior/45004h2012"; % Lake Superior Buoy 45004 2022

data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
dir_t = 15;
g_t = 1.5;

[qi,umag,gust,dir,waveht,avgu,avght,varu,varht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,dir_t,g_t);

if ~isnan(qi)
    stru = "Avg |u| is " + avgu + " with std of " + varu + " m/s";
    strht = "Avg wave height is " + avght + " with std of " + varht + " m";
    
    
    if makeplots
        timeplot = 1:numel(umag);
        
        figure;
        boxplot([umag gust],'Notch','on','Labels',{'|u|','gust'},'Whisker',1)
        title('2022')
        
        figure;
        boxplot(dir,'Notch','on','Labels',{'direction'},'Whisker',1)
        title('2022')
        
        figure;
        boxplot([umag(qi) gust(qi)],'Notch','on','Labels',{'|u|','gust'},'Whisker',1)
        title('ROI 2022')
        
        figure;
        boxplot(dir(qi),'Notch','on','Labels',{'direction'},'Whisker',1)
        title('ROI 2022')
        
        figure;
        subplot(2,1,1)
        plot(timeplot,umag,'-k')
        grid on;
        hold on
        patch([timeplot(qi(1)) timeplot(qi(end)) timeplot(qi(end)) timeplot(qi(1))],[0 0 20 20],'red','FaceAlpha',0.5)
        ylim([0 20])
        xlabel('time (# measurement pts)')
        ylabel('u [m/s]')
        title('2022')
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
        title({stru,strht})
    else
        disp(stru);
        disp(strht);
    end

end

    
