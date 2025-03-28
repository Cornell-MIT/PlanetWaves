clc
clear
close all

addpath('C:\Users\Owner\OneDrive\Documents\00_Main\Work\Github_Repos\umwm_titan\data\Earth')
addpath('C:\Users\Owner\OneDrive\Documents\00_Main\Work\Github_Repos\umwm_titan\validation')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT OBSERVATIONS FROM BUOYS
loc = fullfile('C:\Users\Owner\OneDrive\Documents\00_Main\Work\Github_Repos\umwm_titan\data\Earth\WindFetchLS_45004.csv');
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
p1 = scatter(fax1,Events.U,Events.obs_H,150,Events.F/1000,'filled');
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
% JONSWAP
JS = 4.*sqrt((1.67e-7).*((test_speeds.^2)./9.81).*avg_fetch);
JS_p = 4.*sqrt((1.67e-7).*((test_speeds.^2)./9.81).*p_fetch);
JS_m = 4.*sqrt((1.67e-7).*((test_speeds.^2)./9.81).*m_fetch);

JS_RMS = calc_rms(test_speeds,JS,Events.U,Events.obs_H);
fprintf('RMS of JONSWAP Curve: %f\n',JS_RMS)

p2 = plot(fax1,test_speeds, JS, '--','Color',[0.7 0.7 0.7], 'LineWidth', 3);
p4 = fill(fax1,[test_speeds, fliplr(test_speeds)], [JS_p, fliplr(JS_m)], 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0);
p5 = fill(fax1,[test_speeds, fliplr(test_speeds)], [JS_m, fliplr(JS_p)], 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0);


xlabel('$|u|$ [m/s]','FontSize',25,'interpreter','latex')
ylabel('$H_{1/3}$ [m]','FontSize',25,'interpreter','latex')
%title('Lake Superior: BUOY 45004')
grid on
box on;
xlim([0 20])
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
ylim([0 max(Events.obs_H)+1])

wavetable = readtable('umwm_wind_waveheights.xlsx');

PM = (0.22.*(wavetable.u).^2)./9.81;

PM_RMS = calc_rms(test_speeds,PM,Events.U,Events.obs_H);
fprintf('RMS of Pierson-Moskowtiz Curve: %f\n',PM_RMS)

p6 = plot(wavetable.u,wavetable.umwm,'-ok','MarkerFaceColor','k','LineWidth',5,'DisplayName','UMWM');
p7 = plot(wavetable.u,wavetable.planetwaves,'-or','MarkerFaceColor','r','LineWidth',5,'DisplayName','PlanetWaves');
p8 = plot(wavetable.u,PM,':','Color',[0.7 0.7 0.7],'LineWidth',7,'DisplayName','Pierson-Moskowitz');

% legend([p1 p8 p2 p6 p7 p8],'Observation','Pierson-Moskowitz','JONSWAP','UMWM','PlanetWaves','Location','best')
%exportgraphics(gcf, 'LSCompare.png', 'ContentType', 'vector');

figure('units', 'normalized', 'outerposition', [0 0 1 1]);
plot(wavetable.umwm,wavetable.umwm,'-k','LineWidth',4)
hold on;
plot(wavetable.umwm,wavetable.planetwaves,'o','MarkerFaceColor',[219, 84, 97]./255,'MarkerSize',10)
grid on;
xlabel('UMWM H_{1/3} [m]')
ylabel('PlanetWaves H_{1/3} [m]')

x = wavetable.umwm;
y = wavetable.planetwaves;

SSR = sum((y - x).^2);
SST = sum((y - mean(y)).^2);
R_sqr = 1 - (SSR/SST);

title(['R^2 = ' num2str(R_sqr)])

UMWM_RMS = calc_rms(wavetable.u,wavetable.umwm,Events.U,Events.obs_H);
fprintf('RMS of UMWM Curve: %f\n',UMWM_RMS)
PlanetWave_RMS = calc_rms(wavetable.u,wavetable.planetwaves,Events.U,Events.obs_H);
fprintf('RMS of PlanetWaves Curve: %f\n',PlanetWave_RMS)

%exportgraphics(gcf, 'LSCompare_Rsqr.pdf', 'ContentType', 'vector');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rgb = hex2rgb(hex)
    hex = reshape(hex, [], 6);
    rgb = reshape(sscanf(hex.', '%2x'), [], 3) / 255;
end


function RMS = calc_rms(x_model,y_model,x_data,y_data)

    y_interp = interp1(x_model,y_model,x_data,'linear');
    RMS = sqrt(mean((y_data - y_interp).^2));

end
