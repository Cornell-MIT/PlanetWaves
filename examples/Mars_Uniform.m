clc
clear
close all

% MARS JEZERO CRATER

addpath(fullfile('..','planetwaves'))  

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL
test_speeds = [9:20];
planet_to_run = 'Mars';
time_to_run = 5;   
wind_direction = 0;      
grid_resolution = [7.5*1000 7.5*1000];
zDep = 10.*ones(12,12);
buoy_loc = [6 6];
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,0,zDep,buoy_loc);
Model.z_data = 10;
Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);   

figure;
hold on;
% xlabel('$|u|$ [m/s]','FontSize',25,'interpreter','latex')
% ylabel('$H_{1/3}$ [m]','FontSize',25,'interpreter','latex')
% title('Mars: Jezero Crater Lake')
% grid on
% box on;
% xlim([0 20])
% set(gca,'FontSize',16)
% set(gca,'FontWeight','bold')

save_time = NaN(numel(test_speeds),time_to_run);
for i = 1:numel(test_speeds)

    Wind.speed = test_speeds(i);

    [avgHsig, ~, E_spec, ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 
    save_time(i,:) = avgHsig;
    save_avgHsig = avgHsig;
    save_avgHsig(avgHsig==0) = [];
    save_E_spec = E_spec(~cellfun('isempty',E_spec));
    if sum(avgHsig) ~= 0
        spectrogram.wave_height(i) = save_avgHsig(end);
        plot(1:numel(avgHsig),avgHsig,'-','DisplayName',num2str(Wind.speed))
        drawnow;
        hold on;
        yline(spectrogram.wave_height(i),'--k',num2str(Wind.speed))
        drawnow;
    else
        spectrogram.wave_height(i) = 0;
    end
    spectrogram.energy{i} = save_E_spec{end};
    spectrogram.wind(i) = test_speeds(i);

    %save('test_Mars_loop2.mat','spectrogram','save_time')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT MODEL
    % if spectrogram.wave_height(i) > 0
    %     plot(test_speeds(i), spectrogram.wave_height(i),'--s','LineWidth',1,'MarkerFaceColor','#FF6F59','MarkerSize',15,'MarkerEdgeColor','#FF6F59','Color','#FF6F59')
    %     drawnow;
    % end

end
legend('show')


% % add dashed line between points
% plot(test_speeds, spectrogram.wave_height,'--s','LineWidth',2,'MarkerFaceColor','#FF6F59','MarkerSize',15,'MarkerEdgeColor','#FF6F59','MarkerEdgeColor','#FF6F59','Color','#FF6F59')
% 
% Planet.surface_press = 4*50000; 
% 
% for i = 1:numel(test_speeds)
% 
%     Wind.speed = test_speeds(i);
% 
%     [avgHsig, ~, E_spec, ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 
%     save_avgHsig = avgHsig;
%     save_avgHsig(avgHsig==0) = [];
%     save_E_spec = E_spec(~cellfun('isempty',E_spec));
%     if sum(avgHsig) ~= 0
%         spectrogram.wave_height(i) = save_avgHsig(end);
%     else
%         spectrogram.wave_height(i) = 0;
%     end
%     spectrogram.energy{i} = save_E_spec{end};
%     spectrogram.wind(i) = test_speeds(i);
%     %save('Mars2_loop.mat','spectrogram')
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % PLOT MODEL
%     % if spectrogram.wave_height(i) > 0
%     %     plot(test_speeds(i), spectrogram.wave_height(i),'--s','LineWidth',1,'MarkerFaceColor','#2B3A67','MarkerSize',15,'MarkerEdgeColor','#2B3A67','Color','#2B3A67')
%     %     drawnow;
%     % end
% 
% end
% 
% 
% % add dashed line between points
% % plot(test_speeds, spectrogram.wave_height,'--s','LineWidth',2,'MarkerFaceColor','#2B3A67','MarkerSize',15,'MarkerEdgeColor','#2B3A67','MarkerEdgeColor','#2B3A67','Color','#2B3A67')
% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rgb = hex2rgb(hex)
    hex = reshape(hex, [], 6);
    rgb = reshape(sscanf(hex.', '%2x'), [], 3) / 255;
end