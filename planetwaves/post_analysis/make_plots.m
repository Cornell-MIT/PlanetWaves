function make_plots(Planet,Model,Wind,test_speeds,myHsig, htgrid,energy,wn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKES PLOTS FOR WAVEFIELD RESULTS OF:
% (1) Signifigant wave height growth with time (good to check this for numerical ringing)
% (2) Wavefield height in the basin with the wave height at the buoy -- gif
% (3) Wave height vs wind speed at the buoy
% (4) Wavenumber frequency spectrum with average of linear fit to high frequency slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    for i = 1:numel(test_speeds)
        plot(myHsig{i}, '-', 'LineWidth', 3, 'DisplayName', num2str(test_speeds(i)))
        hold on
    end
    grid on;
    legend('show', 'Location', 'northwest','interpreter','latex');
    title(['Waves on',' ',Planet.name],'interpreter','latex');
    xlabel('model time step [$\Delta$ t]','interpreter','latex')
    ylabel('significant wave height [m]','interpreter','latex')

    for k = 1:numel(test_speeds)
        if ~isempty(htgrid{k}{end})
            buoy_waves(k) = htgrid{k}{end}(Model.long,Model.lat);
        end
    end

    frame = 1;

    figure('units','normalized','outerposition',[0 0 1 1],'Color',[1 1 1])
    for speed = 1:numel(test_speeds)
        if ~isempty(htgrid{speed}{end})
            subplot(1,3,3)
            plot(test_speeds(buoy_waves>0),buoy_waves(buoy_waves>0),'-sb','LineWidth',1,'MarkerFaceColor','b')
            hold on
            plot(test_speeds(speed),htgrid{speed}{end}(Model.long,Model.lat),'-sr','LineWidth',1,'MarkerFaceColor','r')
            xlim([0 max(test_speeds)+0.5])
        
            hold off;
            xlabel('u [m/s]','interpreter','latex')
            ylabel('Hsig [m]','interpreter','latex')
            grid on;
        
            plot_grid = htgrid{speed}{end}';
            plot_grid(isnan(plot_grid)) = 0;
            plot_alpha_data = ones(size(plot_grid));
            plot_alpha_data(plot_grid==0) = 0;
        
            ax1 = subplot(1,3,[1,2]);
            h1 = imagesc(plot_grid);
            colormap linspecer
            xlabel('longitude [km]')
            ylabel('latitude [km]')
            title(sprintf('u = %i m/s',test_speeds(speed)))
            c1 = colorbar;
            c1.Label.String = 'Hsig [m]';
            clim([0 ceil(max(buoy_waves))])
            set(h1, 'AlphaData', plot_alpha_data);
            hold on;
        
            %contour(Model.bathy_map,linspace(min(min(Model.bathy_map)),max(max(Model.bathy_map)),10),':k','LineWidth',2)
            contour(plot_grid,[0:max(buoy_waves)/10:max(buoy_waves)],'-k','LineWidth',2)
        
            grid on
            new_xtick = get(gca, 'XTick')*(Model.gridX)/1000;
            new_ytick = get(gca, 'YTick')*(Model.gridY)/1000;
            set(gca, 'XTick',  get(gca, 'XTick'), 'XTickLabel', arrayfun(@(x) sprintf('%d', x), new_xtick, 'UniformOutput', false));
            set(gca, 'YTick',  get(gca, 'YTick'), 'YTickLabel', arrayfun(@(y) sprintf('%d', y), new_ytick, 'UniformOutput', false));
        
        
            
            plot(Model.long,Model.lat,'pentagram','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',20)
        
            [wx,wy] = pol2cart(Wind.dir,1);
            size_plot = size(htgrid{speed}{end});
        
            quiver(Model.long,size_plot(1).*0.85,(speed/5)*wx, (speed/5)*wy, 'k', 'MaxHeadSize', 1);
            % for i = 1:Model.LonDim
            %     for j = 1:Model.LatDim
            %         if Model.bathy_map(j,i) <=0
            %             quiver(i, j, (speed/5)*wx, (speed/5)*wy, 'k', 'MaxHeadSize', 1);
            %         end
            %     end
            % end
            set(ax1,'Ydir','reverse')
            %set(ax1,'Xdir','reverse')
            drawnow
            if frame == 1
                gif([Planet.name,'_waves.gif'],'DelayTime',1,'overwrite',true)
                frame = 0;
            else
                gif
            end
        end
    end
    
    exportgraphics(gcf, [Planet.name,'_waves.pdf'], 'ContentType', 'vector');

    figure('units','normalized','outerposition',[0 0 1 1])
    plot(test_speeds,buoy_waves,'-sb','LineWidth',1,'MarkerFaceColor','b')
    ylim([0 max(buoy_waves)+1])
    xlabel('u [m/s]','interpreter','latex')
    ylabel('Hsig [m]','interpreter','latex')
    grid on;
    
    
    

    dang =  360/Model.Dirdim; 
    
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    
    
    color1 = [229, 211, 82] / 255;
    color2 = [172, 57, 49] / 255;
    
    
    num_colors = numel(test_speeds);
    
    cmap = [linspace(color1(1), color2(1), num_colors)', ...
                       linspace(color1(2), color2(2), num_colors)', ...
                       linspace(color1(3), color2(3), num_colors)'];
    
    
    set(gca,'ColorOrder',colormap(cmap))
    hold on;
    
    for i = 1:length(test_speeds)
        
        if ~isempty(htgrid{i}{end})
            k_wave = wn{i};
            e_wave = energy{i}*dang;
            plot(k_wave,e_wave,'-','DisplayName',num2str(test_speeds(i)),'LineWidth',5)
            hold on;
        
            [~,max_ind] = max(e_wave);
            inds = k_wave >= k_wave(max_ind);
        
            X = k_wave(inds); 
            logX = log(X);
            Y = e_wave(inds); 
            logY = log(Y);
          
            validInds = ~isinf(logY);
            logX = logX(validInds);
            X = X(validInds);
            Y = Y(validInds);
            logY = logY(validInds);
        
            coeffs = polyfit(logX, logY, 1);
            y_int(i) = coeffs(2);
            n = abs(coeffs(1));
            fittedLogY = polyval(coeffs, logX);
        
            residuals = logY - fittedLogY;
            stdResiduals = std(residuals,"omitmissing");
            Sxx = sum((logX - mean(logX)).^2);
            n_std = stdResiduals / sqrt(Sxx);
        
            n_loop(i) = n;
            n_std_loop(i) = n_std(1);
        else
            n_loop(i) = NaN;
            n_std_loop(i) = NaN;
            y_int(i) = NaN;
        end
    end
    
    
    
    meanY = X.^(-mean(n_loop,"omitmissing"))*exp(mean(y_int,"omitmissing"));
    plot(X, meanY, '--k', 'LineWidth', 3, 'DisplayName', 'fit')
    
    legend('show')
    grid on;
    xlabel('k [m-1]', 'FontSize', 25, 'interpreter', 'latex')
    ylabel('E(k)', 'FontSize', 25, 'interpreter', 'latex')
    set(gca, 'Xscale', 'log')
    set(gca, 'Yscale', 'log')
    ylim([1e-6 1e5])
    set(gca, 'FontSize', 16)
    set(gca, 'FontWeight', 'bold')
    box on;
    
    
    title(sprintf('High frequency tail: f^n where n is: %f ± %f\n', mean(n_loop,"omitmissing"), mean(n_std_loop,"omitmissing")))


end