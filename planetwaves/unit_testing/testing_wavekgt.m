clc
clear
close all


% load('.\EarthAnalysis\GreatLakes\LakeData\LakeSuperior_cleaned.mat')
% LS = squeeze(LS);
% LS_orig = -LS;
% resizeFactor = 0.002;
% % from find_fetch.py
% blon = 1729;
% blat = 6618;
% gridcellsizeX = 4542.948547909539*(1/resizeFactor);
% gridcellsizeY = 92.66280063299297*(1/resizeFactor);
% pos =  [blon, blat];
% pos_orig = pos;
% pos = round(pos * resizeFactor);
% LS = imresize(LS, resizeFactor, "bilinear");
% size_lake = size(LS);
% alphaData = ones(size_lake);
% alphaData(LS==0) = 0;
% D = -LS;
% LS = round(LS);

o = 35;
min_freq = 0.05;                                               % minimum frequency to model
max_freq = 35; 
rho_liquid = 997;
surface_tension = 0.072;  
g = 9.81;

dlnf=(log(max_freq)-log(min_freq))/(o-1);                      % frequency step size for log normal distribution
freq = exp(log(min_freq)+(0:o-1)*dlnf);                           % frequencies for spectrum
D = 0.0001:0.001:300;

for i = 1:numel(freq)
    f = freq(i);

    
    
    T = 1./f;
    
    shallow_L = T.*sqrt(g*D(1));
    shallow_wn = (2*pi)./shallow_L;
    
    deep_L = (g/(2*pi)).*(T.^2);
    deep_wn = (2*pi)./deep_L;
    
    eckart_wn = deep_L.*sqrt(tanh((4*(pi^2).*D)./(g*(T^2))));
    eckart_wn = (2*pi)./eckart_wn;
    
    for x = 1:numel(D)
        wn(:,x) = wavekgt(f,D(x),g,surface_tension,rho_liquid); 
    end
    
   
    % figure;
    if i == 1
        figure;
        hold on
        plot(D,wn(1,:),'Color', 'b', 'LineWidth', 1.5);
        hold on;    
        plot(D,eckart_wn,'-g','LineWidth', 1.5)
        yline(shallow_wn,'--k','shallow')
        yline(deep_wn,'--k','deep')
        hold off
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')
        title(sprintf('freq: %0.2f Hz',f))
        ylabel('wavenumber (k)')
        xlabel('depth')
        legend('estimated k','eckart k','location','best')
        ylim padded
        gif('wavenumber_estimate.gif','DelayTime',0.25)
        hold off;
    else
        plot(D,wn(1,:),'Color', 'b', 'LineWidth', 1.5);
        hold on;    
        plot(D,eckart_wn,'-g','LineWidth', 1.5)
        yline(shallow_wn,'--k','shallow')
        yline(deep_wn,'--k','deep')
        hold off
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')
        title(sprintf('freq: %0.2f Hz',f))
        ylabel('wavenumber (k)')
        xlabel('depth')
        legend('estimated k','eckart k','location','best')
        ylim padded
        gif
        hold off
    end

    
    error(i,:) = eckart_wn./wn(1,:);
  
    
end

figure;
contourf(D,freq,error)
colormap(hot)
title('eckart k/estimated k')
colorbar
xlabel('depth')
ylabel('frequency')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')