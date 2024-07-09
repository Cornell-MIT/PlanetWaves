clc
clear
close all

% MODELED VS OBSERVED AT LAKE SUPERIOR DEEPWATER BUOY

addpath(fullfile('..','planetwaves'))  
addpath(fullfile('..','planetwaves','pre_analysis'))
addpath(fullfile('..','data','Earth'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT OBSERVATIONS FROM BUOYS
loc = fullfile('..','data','Earth','WindFetchLS_45004.csv');
A = readtable(loc,'VariableNamingRule', 'preserve');
A.Properties.VariableNames = {'Dir_deg','Fetch_m'};
dir_to_fetch = containers.Map(A.Dir_deg,A.Fetch_m);

[u,h,angle] = table_quiet_times();

WIND_SPEED = [];
WIND_FETCH = [];
WIND_HEIGHT = [];

for i = 1:numel(u)
    
    if numel(u{i}) > 1
        WIND_SPEED = [WIND_SPEED; u{i}];
        WIND_HEIGHT = [WIND_HEIGHT; h{i}];
        for j = 1:numel(u{i})
            WIND_FETCH = [WIND_FETCH; dir_to_fetch(angle{i}(j))];
        end
    end

end


Events = table(WIND_SPEED,WIND_FETCH,WIND_HEIGHT);
Events.Properties.VariableNames = {'U','F','obs_H'};
Events = sortrows(Events,'obs_H');

% PLOT OBSERVATION 
mymap = [
    hex2rgb('1d55c7')
    hex2rgb('1d85f1')
    hex2rgb('5ca9ff')
    hex2rgb('9cd3ff')
    ];

figure('units', 'normalized', 'outerposition', [0 0 1 1]);
fax1 = axes;
scatter(fax1,Events.U,Events.obs_H,150,Events.F/1000,'filled');
cb = colorbar;
colormap(mymap)
ylabel(cb,'Fetch [km]','FontSize',16,'Rotation',270)
hold on;
grid on;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THEORY

avg_fetch = mean(Events.F);
std_fetch = std(Events.F);
p_fetch = avg_fetch + std_fetch;
m_fetch = avg_fetch - std_fetch;

% PIERSON-MOSKOWITZ
test_speeds = 0:20;
PM = (0.22.*(test_speeds).^2)./9.81;
% JONSWAP
JS = 4.*sqrt((1.67e-7).*((test_speeds.^2)./9.81).*avg_fetch);
JS_p = 4.*sqrt((1.67e-7).*((test_speeds.^2)./9.81).*p_fetch);
JS_m = 4.*sqrt((1.67e-7).*((test_speeds.^2)./9.81).*m_fetch);

plot(fax1,test_speeds, JS, '--k', 'LineWidth', 3);
plot(fax1,test_speeds,PM,':k','LineWidth',3)
fill(fax1,[test_speeds, fliplr(test_speeds)], [JS_p, fliplr(JS_m)], 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0);
fill(fax1,[test_speeds, fliplr(test_speeds)], [JS_m, fliplr(JS_p)], 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0);


xlabel('$|u|$ [m/s]','FontSize',25,'interpreter','latex')
ylabel('$H_{1/3}$ [m]','FontSize',25,'interpreter','latex')
title('Lake Superior: BUOY 45004')
grid on
box on;
xlim([0 20])
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
ylim([0 max(Events.obs_H)+1])







% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL
test_speeds = 1:20;
planet_to_run = 'Earth';
time_to_run = 1000;   
wind_direction = 0;      
grid_resolution = [20*1000 20*1000];
zDep = 273.5.*ones(12,12);
buoy_loc = [6 6];
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,0,zDep,buoy_loc);
Model.z_data = 3.6;
Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);   


directions = 0:(2*pi)/Model.Dirdim:2*pi;
directions(end) = [];

f_vec=(log(Model.max_freq)-log(Model.min_freq))/(Model.Fdim-1);             
f = exp(log(Model.min_freq)+(0:Model.Fdim-1)*f_vec);  

spectrogram.freq = f;
spectrogram.dir = directions;

dang =  360/Model.Dirdim;

figure('units', 'normalized', 'outerposition', [0 0 1 1]);
fax2 = axes;
for i = 1:numel(test_speeds)

    Wind.speed = test_speeds(i);

    [avgHsig, ~, E_spec, ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 
    save_avgHsig = avgHsig;
    save_avgHsig(avgHsig==0) = [];
    save_E_spec = E_spec(~cellfun('isempty',E_spec));
    if sum(avgHsig) ~= 0
        spectrogram.wave_height(i) = save_avgHsig(end);
        plot(fax2,1:numel(avgHsig),avgHsig,'-','DisplayName',num2str(Wind.speed))
        drawnow;
        hold on;
        yline(fax2,spectrogram.wave_height(i),'--k',num2str(Wind.speed))
        drawnow;
    else
        spectrogram.wave_height(i) = 0;
    end
    spectrogram.energy{i} = save_E_spec{end};
    spectrogram.wind(i) = test_speeds(i);
    save('EarthCompare.mat','spectrogram')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT MODEL
    % if spectrogram.wave_height(i) > 0
    %     plot(test_speeds(i), spectrogram.wave_height(i),'--sr','LineWidth',1,'MarkerFaceColor','#c7391d','MarkerSize',15,'MarkerEdgeColor','#c7391d')
    %     drawnow;
    % end

end


% add dashed line between points
plot(fax1,test_speeds, spectrogram.wave_height,'--sr','LineWidth',2,'MarkerFaceColor','#c7391d','MarkerSize',15,'MarkerEdgeColor','#c7391d','MarkerEdgeColor','#c7391d')

% figure('units', 'normalized', 'outerposition', [0 0 1 1]);
% cmap = linspecer(numel(test_speeds));
% set(gca,'ColorOrder',cmap)
% hold on;
% 
% for i = 3:11
%     plot(f,squeeze(sum(spectrogram.energy{1,i}(Model.long,Model.lat,:,:),4))*dang,'*-','DisplayName',num2str(spectrogram.wind(i)),'LineWidth',2)
%     hold on;
% end
% legend('show')
% grid on;
% xlabel('frequency [Hz]','FontSize',25,'interpreter','latex')
% ylabel('Omni-Directional E','FontSize',25,'interpreter','latex')
% set(gca,'Xscale','log')

%save('Earth_Compare_LakeSuperior_45004.mat','spectrogram')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rgb = hex2rgb(hex)
    hex = reshape(hex, [], 6);
    rgb = reshape(sscanf(hex.', '%2x'), [], 3) / 255;
end
