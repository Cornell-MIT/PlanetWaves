function [x,y,Xmesh, Ymesh, zDep] = make_fake_shoreline(radius,num_points,noise_amplitude)

addpath('C:\Users\Owner\OneDrive\Desktop\Main\Work\Github_Repos\umwm_titan\data\Titan\TitanLakes\Bathymetries\bathtub_bathy')
theta = linspace(0, 2*pi, num_points); % Angle values

% basic circle
x = radius * cos(theta);
y = radius * sin(theta);

% Add irregularity to the left side
for i = 1:length(x)
    noise_amp = noise_amplitude * (rand - 0.5)*((abs(x(i)-max(x)))/(min(x) - max(x)));
    x(i) = x(i) + noise_amp * cos(theta(i));
    y(i) = y(i) + noise_amp * sin(theta(i));

end

figure;
plot(x,y)
[Xmesh,Ymesh,zDep] = make_bathtub_lake(0.002,[x' y']);
figure;
surf(Xmesh, Ymesh, zDep,'EdgeColor','none');
end