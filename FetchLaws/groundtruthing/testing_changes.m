clc
clear
close all

% checking that new edits and organization haven't changed the original model

% --- %

% CCg %
% load('04_14_24_run_from_old.mat','Ccg')
% CCg_original = Ccg;
% load('04_14_24_run_from_main.mat','Ccg')
% 
% local_max = 0;
% global_max = 0;
% 
% 
% for k = 1:25
%     for m = 1:288
%     %image(CCg_original(:,:,j,k) - Ccg(:,:,j,k))
%     local_max = max(max(CCg_original(:,:,k,m) - Ccg(:,:,k,m)));
%         if local_max > global_max
%            global_max = local_max;
%         end
%     %title(['k ',num2str(k),' m ',num2str(m)])
%     %pause(0.1)
%     end
% end



% Cg %
% load('04_14_24_run_from_old.mat','Cg')
% Cg_original = Cg;
% load('04_14_24_run_from_main.mat','Cg')
% 
% local_max = 0;
% global_max = 0;
% 
% 
% for k = 1:25
%     for m = 1:288
%     local_max = max(max(Cg_original(:,:,k,m) - Cg(:,:,k,m)));
%         if local_max > global_max
%            global_max = local_max;
%         end
%     end
% end

%D
% load('04_14_24_run_from_old.mat','D')
% D_original = D;
% load('04_14_24_run_from_main.mat','D')
% 
% local_max = 0;
% global_max = 0;
% 
% 
% for k = 1:25
%     for m = 1:288
%     local_max = max(max(D_original(:,:,k,m) - D(:,:,k,m)));
%         if local_max > global_max
%            global_max = local_max;
%         end
%     end
% end



disp('done')