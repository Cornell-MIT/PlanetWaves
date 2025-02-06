clc
clear
close all

load('ontario_lacus_bathtub_0.002.mat','x','y')
[x,y] = smooth_path(x,y,10);

wind_direction = {'West -> East (GCM)', 'South -> North','East -> West ', 'North -> South'};
clr = {'blue','red','black','green'};

load('zero.mat','psi','sinu')
x1 = psi;
y1 = 1 - sinu;
clearvars psi sinu

load('pi2.mat','psi','sinu')
x2 = psi;
y2 = 1 - sinu;
clearvars psi sinu

load('pi.mat','psi','sinu')
x3 = psi;
y3 = 1 - sinu;

clearvars psi sinu
load('3pi2.mat','psi','sinu')
x4 = psi;
y4 = 1 - sinu;
clearvars psi sinu

psi = [x1; x2; x3; x4];
sinu = [y1; y2; y3; y4];

 
figure;
scatter(x,y,100,sinu(2,:),"filled")
colorbar
title('sinuuosity')

figure;
for i = 1:numel(wind_direction)
    diff_sinu(i,:) = psi(i,:)./sinu(i,:);
    scatter(psi(i,:),sinu(i,:),50,clr{i},'filled','DisplayName',wind_direction{i})
    hold on;
end
xlabel('diffusivity')
ylabel('sinuuosity')
legend('show')

figure;
for i = 1:numel(wind_direction)
    histogram(real(psi(i,:)),[-1:0.1:1],'DisplayName',wind_direction{i})
    hold on
end
legend('show')
xlabel('diffusivity')

figure;
hold on
for i = 1:numel(wind_direction)
    diff = psi(i,:);
    sinuu = sinu(i,:);
    pos = diff>0;
    median_group = median(sinuu(pos),'omitmissing');
    std_group = std(sinuu(pos),'omitmissing');
    xrange = linspace(-10, 10,10);
    scatter(diff(pos), sinuu(pos),10,clr{i},'filled','DisplayName',wind_direction{i})
    yline(median_group,'-','Color',clr{i},'LineWidth',5)
    fill([xrange, fliplr(xrange)], [(median_group + std_group) * ones(size(xrange)), ...
                      (median_group - std_group) * ones(size(xrange))], ...
     clr{i}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

end
legend('show')
xlabel('diffusivity')
ylabel('sinuuosity')
title('Smoothing')
xlim([min(min(diff))-0.5 max(max(diff))+0.5])

figure;
hold on
for i = 1:numel(wind_direction)
    diff = psi(i,:);
    sinuu = sinu(i,:);
    neg = diff<0;
    median_group = median(sinuu(neg),'omitmissing');
    std_group = std(sinuu(neg),'omitmissing');
    xrange = linspace(-10, 10,10);
    scatter(diff(neg), sinuu(neg),10,clr{i},'filled','DisplayName',wind_direction{i})
    yline(median_group,'-','Color',clr{i},'LineWidth',5)
    fill([xrange, fliplr(xrange)], [(median_group + std_group) * ones(size(xrange)), ...
                      (median_group - std_group) * ones(size(xrange))], ...
     clr{i}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

end
legend('show')
xlabel('diffusivity')
ylabel('sinuuosity')
title('Roughening')
xlim([min(min(diff))-0.5 max(max(diff))+0.5])


plot_psi = psi(4,:);
max_color = max(plot_psi,[],'omitmissing');
min_color = min(plot_psi,[],'omitmissing');
plot_psi(plot_psi==0) = -999;

figure;
scatter(x(plot_psi~=-999),y(plot_psi~=-999),100,plot_psi(plot_psi~=-999),"filled")
colorbar
clim([-0.5 0.5])
hold on
scatter(x(plot_psi==-999),y(plot_psi==-999),100,plot_psi(plot_psi==-999))
title('Diff')