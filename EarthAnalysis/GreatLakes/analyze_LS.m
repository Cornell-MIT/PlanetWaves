
clc
clear
close all
% OBJECTIVE: Find periods of relative quiescence in Lake Superior Wind/Wave Data to compare to UMWM-Titan

filen = 'BuoyData/Superior/45004h2022'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,avgdir,std_u,std_ht]   = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

RESULTS = [avgu' avgdir' avght'];
RESULTS = rmmissing(RESULTS);

dir_fetch = csvread('WindFetchLS.csv',1,0);
dir_fetch = dictionary(dir_fetch(:,1),dir_fetch(:,2));

avgfetch = dir_fetch(round(RESULTS(:,2)));


RESULTS = [RESULTS avgfetch];

RESULTS = table(round(RESULTS(:,1),1),round(RESULTS(:,2)),round(RESULTS(:,3),1),RESULTS(:,4),'VariableNames',{'WSPD','WDIR','SIGHT','FETCH'});

RESULTS = unique(RESULTS);
RESULTS = [table(2022.*ones(size(RESULTS(:,1))),'VariableNames',{'YEAR'}) RESULTS];

writetable(RESULTS,'LS_2022.csv')

u_roi = umag;
u_roi(logical(~qi)) = NaN;

figure;
plot(1:numel(umag),umag,'-ok')
hold on
plot(1:numel(umag),u_roi,'o','MarkerFaceColor','r')

pm_u = umag;
pm_u(isnan(pm_u)) = [];
[x,i] = sort(umag);


PM = (0.22.*(x).^2)./9.81;

cb = rand(43,3);
szs = (max(std_ht) + max(std_u))*20;
errorbaryes = 0;

% Create scatter plot with error bars
figure;
h1 = gscatter(avgu,avght,round(avgdir),[],[],szs);
grid on;


hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end
plot(x,PM,'--k','LineWidth',3)
grid on;
title(strcat('LakeSuperior:',buoy,'d1980-2022'))
xlabel('|u| [m/s]')
ylabel('sig H [m]')

%$$

filen = 'BuoyData/Superior/45004h2021'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,avgdir,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);


szs = (max(std_ht) + max(std_u))*100;
h2 = scatter(avgu, avght, szs,cb(2,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2020'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h3 = scatter(avgu, avght, szs,cb(3,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2019'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h4 = scatter(avgu, avght, szs,cb(4,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end
    
%%%

filen = 'BuoyData/Superior/45004h2018'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h5 = scatter(avgu, avght, szs,cb(5,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2017'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h6 = scatter(avgu, avght, szs,cb(6,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2016'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h7 = scatter(avgu, avght, szs,cb(7,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2015'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h8 = scatter(avgu, avght, szs,cb(8,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2014'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);


h9 = scatter(avgu, avght, szs,cb(9,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2013'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h10 = scatter(avgu, avght, szs,cb(10,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2012'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h11 = scatter(avgu, avght, szs,cb(11,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2011'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h12 = scatter(avgu, avght, szs,cb(12,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2010'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h13 = scatter(avgu, avght, szs,cb(13,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2009'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h14 = scatter(avgu, avght, szs,cb(14,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2008'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h15 = scatter(avgu, avght, szs,cb(15,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2007'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h16 = scatter(avgu, avght, szs,cb(16,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2006'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h17 = scatter(avgu, avght, szs,cb(17,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2005'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h18 = scatter(avgu, avght, szs,cb(18,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2004'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h19 = scatter(avgu, avght, szs,cb(19,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2003'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h20 = scatter(avgu, avght, szs,cb(20,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2002'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h21 = scatter(avgu, avght, szs,cb(21,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h2001'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h22 = scatter(avgu, avght, szs,cb(22,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

% filen = 'BuoyData/Superior/45004h2000'; % Lake Superior Buoy 45004 2022
% yr = filen(end-3:end); % year number for plots
% buoy = filen(end-9:end-5); % buoy number
% data_cadence = 6; % measurement every 10 minute [# pts/hr]
% window_size = 10; % window size in [hrs]
% u_t = 3; % maxmimum possible change in magnitude of wind
% dir_t = 15; % maximum possible change in direction over the window
% g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts
% 
% [qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);
% 
% 
% h23 = scatter(avgu, avght, szs,'filled',cb(23,:));
% hold on;
% errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'r');
% errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'b');

%%%

filen = 'BuoyData/Superior/45004h1999'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h24 = scatter(avgu, avght, szs,cb(24,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h1998'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h25 = scatter(avgu, avght, szs,cb(25,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h1997'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h26 = scatter(avgu, avght, szs,cb(26,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h1996'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h27 = scatter(avgu, avght, szs,cb(27,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h1995'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h28 = scatter(avgu, avght, szs,cb(28,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h1994'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h29 = scatter(avgu, avght, szs,cb(29,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h1993'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h30 = scatter(avgu, avght, szs,cb(30,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h1992'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h31 = scatter(avgu, avght, szs,cb(31,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

% filen = 'BuoyData/Superior/45004h1991'; % Lake Superior Buoy 45004 2022
% yr = filen(end-3:end); % year number for plots
% buoy = filen(end-9:end-5); % buoy number
% data_cadence = 6; % measurement every 10 minute [# pts/hr]
% window_size = 10; % window size in [hrs]
% u_t = 3; % maxmimum possible change in magnitude of wind
% dir_t = 15; % maximum possible change in direction over the window
% g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts
% 
% [qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);
% 
% 
% h32 = scatter(avgu, avght, szs,cb(32,:),'filled');
% hold on;
% errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'r');
% errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'b');

%%%

filen = 'BuoyData/Superior/45004h1990'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h33 = scatter(avgu, avght, szs,cb(33,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h1989'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h34 = scatter(avgu, avght, szs,cb(34,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h1988'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h35 = scatter(avgu, avght, szs,cb(35,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h1987'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h36 = scatter(avgu, avght, szs,cb(36,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h1986'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h37 = scatter(avgu, avght, szs,cb(37,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h1985'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h38 = scatter(avgu, avght, szs,cb(38,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h1984'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h39 = scatter(avgu, avght, szs,cb(39,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h1983'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h40 = scatter(avgu, avght, szs,cb(40,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h1982'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h41 = scatter(avgu, avght, szs,cb(41,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

%%%

filen = 'BuoyData/Superior/45004h1981'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h42 = scatter(avgu, avght, szs,cb(42,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end
%%%

filen = 'BuoyData/Superior/45004h1980'; % Lake Superior Buoy 45004 2022
yr = filen(end-3:end); % year number for plots
buoy = filen(end-9:end-5); % buoy number
data_cadence = 6; % measurement every 10 minute [# pts/hr]
window_size = 10; % window size in [hrs]
u_t = 3; % maxmimum possible change in magnitude of wind
dir_t = 15; % maximum possible change in direction over the window
g_t = 1.5; % maximum percent above the mean wind speed for concurrent gusts

[qi,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filen,data_cadence,window_size,u_t,dir_t,g_t);

szs = (max(std_ht) + max(std_u))*100;
h43 = scatter(avgu, avght, szs,cb(43,:),'filled');
hold on;
if errorbaryes 
    errorbar(avgu, avght, std_ht.*ones(size(avght)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
    errorbar(avgu, avght, std_u.*ones(size(avgu)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'k');
end

% h_in = [h1 h2 h3 h4 h24 h30 h43];
% legend(h_in,'2022' ,'2021', '2020', '2019', '1998', '1990', '1984')
