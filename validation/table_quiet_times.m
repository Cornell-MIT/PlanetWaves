function [u,h,angle] = table_quiet_times()

years = 2002:2022;
data = cell(1, length(years));
u = data;
h = data;
angle = data;

cwd = pwd;
cd('..');
cd('data/Earth/GreatLakes/LakeSuperior/45004_Buoy');

for i = 1:length(years)
    data{i} = analyze_buoy_data(['45004h' num2str(years(i))]);
end


for i = 1:length(years)

    if ~isempty(data{i})
        u{i} = data{i}.WSPD;
        h{i} = data{i}.SIGHT;
        angle{i} = data{i}.WDIR;
    else
        u{i} = NaN;
        h{i} = NaN;
        angle{i} = NaN;
    end
end




cd(cwd)

end
