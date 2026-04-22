clc
clear all
%close all


%% COLD GRAVEL

R = 950./670. - 1; %(sg-1) for Titan

g = 1.35; %m/s^2
%alpha_s = 0.1*R/((1.+R)^0.77);
%alpha_b = 20./((1.+R)*sqrt(R));
%alpha_h = 0.18*(1.+R)^0.77;
%alpha_y = 0.007/(1.+R); %eq 25
alpha_s = 0.11*R/((1.+R)^0.73);
alpha_b = 18./((1.+R)*sqrt(R));
alpha_h = 0.22*(1.+R)^0.73;
alpha_y = 0.01/(1.+R);
%n_b = 0.05;
%n_s = -0.31;
%n_h = -0.006;
n_b = 0.06;
n_s = -0.33;
n_h = -0.02;
n_y = n_b + 1.5*(n_s + n_h) + 1.; %eq 21

set(gcf, 'visible', 'on')

B = 175:25:700;
D50_gravel = zeros(size(B));
Q_gravel = zeros(size(B));
Qs_gravel = zeros(size(B));

j = 0;
%for S = 0.0004:0.0009:0.002 %slope
for S = [0.0004,0.0013,0.002]
    j = j + 1;
    i = 0;
    for B = 175:25:700 %channel width
        i = i + 1;
        D50_gravel(j,i) = (B/alpha_b)*(alpha_s/S)^((n_b+0.4)/n_s);
        Q_gravel(j,i) = ((B/alpha_b)*(D50_gravel(j,i)^(2.5*n_b))*g^(0.2+0.5*n_b))^(1./(n_b+0.4));
        Qs_gravel(j,i) = alpha_y*(Q_gravel(j,i)^n_y)*(g^((1.-n_y)/2.))*(D50_gravel(j,i)^(2.5*(1.-n_y)));
    end
end

B = 175:25:700;

figure
box on
plot(B,D50_gravel*100,linewidth=1)
rectangle('Position',[175 1.6 700-175 61-1.6], 'EdgeColor','r')
xlim([min(B),max(B)])
xlabel('channel width (m)')
ylabel('the dimensional bed grain size (cm)')
title('D_{50} in Saraswati Flumen [cold,gravel]')
legend('S=0.0004','S=0.0013','S=0.0020','Location','northwest')
%saveas(gcf,'D50_Vid_Flumina_cold_gravel.png')

figure
plot(B,Q_gravel,linewidth=1)
rectangle('Position',[175 35 700-175 2300-35], 'EdgeColor','r')
xlim([min(B),max(B)])
xlabel('channel width (m)')
ylabel('flow discharge (m^3/s)')
title('Q in Saraswati Flumen [cold,gravel]')
legend('S=0.0004','S=0.0013','S=0.0020','Location','northwest')

figure
plot(B,Qs_gravel,linewidth=1)
rectangle('Position',[175 0.00065 700-175 0.36-0.00065], 'EdgeColor','r')
xlim([min(B),max(B)])
xlabel('channel width (m)')
ylabel('sediment flow discharge (m^3/s)')
title('Qs in Saraswati Flumen [cold,gravel]')
legend('S=0.0004','S=0.0013','S=0.0020','Location','northwest')

qw = Q_gravel;
nyu = diffusivity(Q_gravel,B);

figure
plot(B,nyu,linewidth=3)
xlim([min(B),max(B)])
ax=gca;
ax.LineWidth=2;
ax.FontSize = 20;
xlabel('channel width (m)')
ylabel('diffusivity (m^2/year)')
title('Saraswati Flumen [cold, bedload-dominated river]')
legend('S=0.0004','S=0.0013','S=0.0020','Location','northwest')

outfile=['diffusivity_SF.pdf'];
%gcf.Units = 'centimeters';
%gcf.Position = [5,5,15,10];
set(gcf, 'PaperOrientation', 'landscape');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COLD SAND
R = 950./670. - 1; %(sg-1) for Titan

g = 1.35; %m/s^2
%alpha_s = 0.02; %0.01*R;
%alpha_b = 0.7; %0.93/sqrt(R);
%alpha_h = 3.7;
alpha_s = 0.01*R;
alpha_b = 0.95/sqrt(R);
alpha_h = 3.8;
n_b = 0.11;
n_s = -0.17;
n_h = -0.06;
n_y = 1.1; %eq 31

m_s = -0.05;
m_b = 0.10;
v = 3.*10.^(-7); %m^2/s

B = 175:25:700;
D50_sand = zeros(size(B));
Q_sand = zeros(size(B));
Qs_sand = zeros(size(B));
Rep = zeros(size(B));
alpha_y = zeros(size(B));

j = 0;
for S = [0.0004,0.0013,0.002] %slope
    j = j + 1;
    i = 0;
    for B = 175:25:700 %channel width
        i = i + 1;
        %D50_sand(j,i) = ((alpha_b/B)*(((S/alpha_s)^((n_b+0.4)/n_s)))*(R*g/(v^2.))^(-0.5*n_b*(m_s/n_s)-0.2*(m_s/n_s)+0.5*m_b))^(1./(1.5*n_b*(m_s/n_s)+0.6*(m_s/n_s)-1.5*m_b-1.));
        D50_sand(j,i) = 0.00035;
        Q_sand(j,i) = ((B/alpha_b)*(g^(0.5*(n_b-m_b)+0.2))*(D50_sand(j,i)^(2.5*n_b-1.5*m_b))*(R^(-0.5*m_b))*(v^m_b))^(1./(n_b+0.4));
        
        Rep(j,i) = (D50_sand(j,i)*sqrt(R*g*D50_sand(j,i)))/v;
        alpha_y(j,i) = 3.6*(10.^(-5))*(Rep(j,i)^(-0.12)); %Eq33 in SI

        Qs_sand(j,i) = alpha_y(j,i)*(Q_sand(j,i)^n_y)*(g^((1.-n_y)/2.))*(D50_sand(j,i)^(2.5*(1.-n_y)));    
    end
end

B = 175:25:700;
nyu = diffusivity(Q_sand,B);

figure
box on
plot(B,D50_sand,linewidth=1)
xlim([min(B),max(B)])
xlabel('channel width (m)')
ylabel('the dimensional bed grain size (cm)')
title('D_{50} in Saraswati Flumen [cold,sand]')
legend('S=0.0004','S=0.0013','S=0.002','Location','northwest')
%saveas(gcf,'D50_Vid_Flumina_cold_sand.png')

plot(B,Q_sand,linewidth=1)
xlim([min(B),max(B)])
rectangle('Position',[175 680 700-175 45000-680], 'EdgeColor','r')
%ylim([0,2000])
xlabel('channel width (m)')
ylabel('flow discharge (m^3/s)')
title('Q in Saraswati Flumen [cold,sand]')
legend('S=0.0004','S=0.0013','S=0.002','Location','northwest')
%saveas(gcf,'Q_Vid_Flumina_cold_sand.png')

plot(B,Qs_sand,linewidth=1)
rectangle('Position',[175 1.2 700-175 9.4-1.2], 'EdgeColor','r')
xlim([min(B),max(B)])
xlabel('channel width (m)')
ylabel('sediment flow discharge (m^3/s)')
title('Qs in Saraswati Flumen [cold,sand]')
legend('S=0.0004','S=0.0013','S=0.002','Location','northwest')
%saveas(gcf,'Qs_Vid_Flumina_cold_sand.png')

figure
plot(B,nyu,linewidth=3);
ax=gca;
ax.LineWidth=2;
ax.FontSize = 20;
xlim([min(B),max(B)])
xlabel('channel width (m)')
ylabel('diffusivity (m^2/year)')
title('Saraswati Flumen [cold,suspended-load-dominated river]')
legend('S=0.0004','S=0.0013','S=0.002','Location','northwest')

outfile=['diffusivity_SF_susp.pdf'];
%gcf.Units = 'centimeters';
%gcf.Position = [5,5,15,10];
set(gcf, 'PaperOrientation', 'landscape');
exportgraphics(gcf, outfile, 'ContentType','vector')

function nyu = diffusivity(qw,B)
    Cf = 10.^(-2);
    phy0 = 0.5;
    %sg = 2.20; %%water ice (sg=2.20), organics (sg=3.33)
    %sg = 0.95/0.54; %warm (Birch+2023)
    sg = 0.95/0.67; %cold (Birch+2023)
    k = 1; %% meandering (k=1), braded (k=0.11)

    nyu = qw*(8.*k*sqrt(Cf))/((1.-phy0)*(sg-1.));
    nyu = nyu*60.*60.*24.*365./B;
end