function plot_against_lake_superior_45004(model_winds,model_waves)
% This code will extract quiet periods in terms of variability in wind speed, direction, and gust from Lake Superior data
% at a deepwater buoy. It will produce plots of the sig wave height vs wind speed with direction information against the 
% first order estimate by Pierson-Moskowitz. It will then extract fetch information for each time there is a wind speed
% wave height pair and plot that versus the JONSWAP prediction given the average +/- std of the fetches. It will also produce
% a box plot of the fetches.

years = 2002:2022;
data = cell(1, length(years));

cwd = pwd;
cd('../data/Earth/GreatLakes/LakeSuperior/BuoyData/45004');

for i = 1:length(years)
    data{i} = analyze_buoy_data(['45004h' num2str(years(i))]);
end



num_colors = 360;
hues = linspace(0, 1, num_colors)';
saturation = ones(num_colors, 1);
value = ones(num_colors, 1);
hsv_colors = [hues, saturation, value];
cyclic_colormap = hsv2rgb(hsv_colors);


figure('units', 'normalized', 'outerposition', [0 0 1 1]);

for i = 1:length(years)

    if ~isempty(data{i})
        u = data{i}.WSPD;
        h = data{i}.SIGHT;
        angle = data{i}.WDIR;

        scatter(u,h,[],angle,'filled')
        colormap(cyclic_colormap);
        cb = colorbar;
        ylabel(cb,'Direction (deg)','FontSize',16,'Rotation',270)
        hold on;
    end
end


x = 0:20;
PM = (0.22.*(x).^2)./9.81;
grid on;
plot(x, PM, '--k', 'LineWidth', 3)
xlabel('wind speed')
ylabel('sig wave height')
title('Pierson-Moskwotiz Lake Superior [2002-2022]')


%%



max_num_fetch = max(cellfun(@numel, data));
myfetch = cell(length(years),max_num_fetch);

mymap = [
    hex2rgb('1d55c7')
    hex2rgb('1d85f1')
    hex2rgb('5ca9ff')
    hex2rgb('9cd3ff')
    ];

figure('units', 'normalized', 'outerposition', [0 0 1 1]);
for i = 1:length(years)

    pad_ends = NaN(1,max_num_fetch);

    if ~isempty(data{i})
        u = data{i}.WSPD;
        h = data{i}.SIGHT;
        fetch = data{i}.FETCH;
    
        pad_ends(1:numel(fetch)) = fetch;
       
        scatter(u,h,150,fetch./1000,'filled');
        cb = colorbar;
        colormap(mymap)
        ylabel(cb,'Fetch [km]','FontSize',16,'Rotation',270)
        hold on;
    end

    myfetch{i} = pad_ends;
end

myfetch = cell2mat(myfetch);

myfetch = reshape(myfetch.', 1, []);
myfetch = myfetch(~isnan(myfetch));

avg_fetch = mean(myfetch);
std_fetch = std(myfetch);
p_fetch = avg_fetch + std_fetch;
m_fetch = avg_fetch - std_fetch;

JS = 4.*sqrt((1.67e-7).*((x.^2)./9.81).*avg_fetch);
JS_p = 4.*sqrt((1.67e-7).*((x.^2)./9.81).*p_fetch);
JS_m = 4.*sqrt((1.67e-7).*((x.^2)./9.81).*m_fetch);

plot(x, JS, '--k', 'LineWidth', 3);
fill([x, fliplr(x)], [JS_p, fliplr(JS_m)], 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0);
fill([x, fliplr(x)], [JS_m, fliplr(JS_p)], 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0);
%load('waveheights.mat') % from model run
%load('ht_sig_04162024.mat','test_speeds')
plot(model_winds, model_waves,'--sr','LineWidth',1,'MarkerFaceColor','#c7391d','MarkerSize',15,'MarkerEdgeColor','#c7391d')
xlabel('$|u|$ [m/s]','FontSize',25,'interpreter','latex')
ylabel('$H_{1/3}$ [m]','FontSize',25,'interpreter','latex')
%title('Lake Superior: BUOY 45004')
grid on
box on;
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')


% figure;
% boxplot(myfetch)
% ylabel('fetch [km]')
% title('lake superior 2002-2022')
% grid on

cd(cwd)
function rgb = hex2rgb(hex)
    hex = reshape(hex, [], 6);
    rgb = reshape(sscanf(hex.', '%2x'), [], 3) / 255;
end

end