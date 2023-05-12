clc
clear
close all

% contour plots of wavenumber vs angle 
addpath 'C:\Users\schne\OneDrive\Desktop\Main\Work\MIT\Github_Repos\umwm_titan\FetchLaws\Titan'
p = 64;
o = 25;
dr = pi/180;                                                               % conversion from degrees to radians
dthd = 360/(p);
th = rad2deg(([0:p-1]-p/2+0.5)*dthd*dr);
f1 = 0.05;                                                                 % minimum frequency
f2 = 35;                                                                   % maximum frequency
% create frequency limits for spectrum
dlnf=(log(f2)-log(f1))/(o-1);                                              % frequency step size for log normal distribution
f = exp(log(f1)+[0:o-1]*dlnf);                                             % frequencies for spectrum
dom = 2*pi*dlnf.*f;                                                        % angular frequency (w = 2pi*f)
idx = 1;

[rr,tt] = meshgrid(th,dom);
for ii = 1:100
    load(['Titan_1_' num2str(ii) '.mat'])
    contourf(th,dom,squeeze(E(5,5,:,:)))
    colorbar
    set(gca, 'YScale', 'log')
    xlabel('direction [deg]')
    ylabel('angular frequency [rad/s]')
    title(['E Time Step: ' num2str(ii)])
    ylim([dom(1) dom(end)])
    xlim([th(1) th(end)])
%     xx = rr.*cos(tt);
%     yy = rr.*sin(tt);
%     hh = polar(xx,yy);
%     hold on
%     contour(xx,yy,squeeze(E(5,5,:,:)))
%     set(hh,'Visible','off')
%     axis off
%     axis image
    drawnow
%     frame = getframe(gcf);
%     im{idx} = frame2im(frame);
%     idx = idx + 1;
end
% filename = 'E.gif';
% idx_end = idx;
% for idx = 1:idx_end-1
%     [A,map] = rgb2ind(im{idx},256);
% if idx == 1
%     imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
% else
%     imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
% end
%     pause(0.1)
% end

