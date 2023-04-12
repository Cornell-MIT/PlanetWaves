clc
clear
close all

%% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Author: Una Schneck (schneck.una@gmail.com)
% 
% This code performs unit testing for the ref variables in the umwm model, valid for:
%
%                   windspeeds = 0.4:1:3.3
%                   m = 31                                                                    
%                   n = 15                                                                  
%                   rho_liquid = 465
%                   nu_liquid = 0.0031/1e4
%                   bathy_map = 100.*ones(m,n)
%
% If you want to unit test for different parameters, re-run UMWM_Titan_FetchLaws and put resulting .mat files in umwm_old folder under unit_test folder
%
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% External function requirements:
%   (1) none
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Code produces:
%   (1) List of passed/failed tests
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Aprox time to run: < 1 minute 
%---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ===========================================================================================================================================================================================================================
% ===========================================================================================================================================================================================================================

% location of .mat files for testing
addpath('C:\Users\schne\OneDrive\Desktop\Main\Work\MIT\Github_Repos\umwm_titan\FetchLaws\unit_test\umwm_mat\umwm_old')
addpath('C:\Users\schne\OneDrive\Desktop\Main\Work\MIT\Github_Repos\umwm_titan\FetchLaws\unit_test\umwm_mat\umwm_new')
 
% ===========================================================================================================================================================================================================================
% -- check that grid size is same 
% ===========================================================================================================================================================================================================================

% m: number of gridpoints in x-direction 
load('Old_1_ref.mat','m'); om = m;
load('New_1_ref.mat','m'); nm = m;
% n: number of gridpoints in y-direction
load('New_1_ref.mat','n'); on = n;
load('New_1_ref.mat','n'); nn = n;

if om ~= nm 
    error('m-dimension for original model has changed. Need to re-run original model for unit testing')
    
elseif on ~= nn
    error('n-dimension for original model has changed. Need to re-run original model for unit testing')
   
end

[xplot,yplot] = meshgrid(1:m,1:n);

passTest = 0; % counts when a test passes
numTest = 0; % counts when a test is attempted

% ===========================================================================================================================================================================================================================
% -- checking D is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_ref.mat','D'); oD = D;
load('New_1_ref.mat','D'); nD = D;

clearvars D

% for o = 1:25
%     for p = 1:287 
%     subplot(1,2,1)
%     surf(xplot',yplot',nD(:,:,o,p),'EdgeColor','k')
%     colormap(jet(256));
%     title('Original')
%     subplot(1,2,2)
%     surf(xplot',yplot',nD(:,:,o,p),'EdgeColor','k')
%     title('New')
%     sgtitle("o "+o+"p "+p)
%     colormap(jet(256));
%     drawnow;
%     pause(0.1)
%     end
% end

problem = 0;
for o = 1:25
    for p = 1:287
        A = oD(:,:,o,p);
        B = nD(:,:,o,p);
        Index{o,p} = find(A~=B);
        if ~isempty(Index{o,p})
            problem = problem + 1;
        end
        %surf(xplot',yplot',abs(nD(:,:,o,p)-oD(:,:,o,p)))
        %title("o "+o+"p "+p)
        %drawnow;
        %pause(0.1)
    end
end

numTest = numTest + 1;

if ~problem
    disp('Depth Test: Pass')
    passTest = passTest + 1;
else
    disp('Depth Test: FAIL < ----------------------------------------')
end

clearvars -except xplot yplot m n passTest numTest

% ===========================================================================================================================================================================================================================
% -- checking dwn is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_ref.mat','dwn'); odwn = dwn;
load('New_1_ref.mat','dwn'); ndwn = dwn;

clearvars dwn

problem = 0;
for o = 1:25
    for p = 1:287
        A = odwn(:,:,o,p);
        B = ndwn(:,:,o,p);
        Index{o,p} = find(A~=B);
        if ~isempty(Index{o,p})
            problem = problem + 1;
        end
        %surf(xplot',yplot',abs(ndwn(:,:,o,p)-odwn(:,:,o,p)))
        %title("o "+o+"p "+p)
        %drawnow;
        %pause(0.1)
    end
end

numTest = numTest + 1;

if ~problem
    disp('dwn Test: Pass')
    passTest = passTest + 1;
else
    disp('dwn Test: FAIL < ----------------------------------------')
end

clearvars -except xplot yplot m n passTest numTest

% ===========================================================================================================================================================================================================================
% --- checking if wn is the same  
% ===========================================================================================================================================================================================================================

load('Old_1_ref.mat','wn'); own = wn;
load('New_1_ref.mat','wn'); nwn = wn;

clearvars wn

problem = 0;
for o = 1:25
    for p = 1:287
        A = own(:,:,o,p);
        B = nwn(:,:,o,p);
        Index{o,p} = find(A~=B);
        if ~isempty(Index{o,p})
            problem = problem + 1;
        end
%         surf(xplot',yplot',abs(nwn(:,:,o,p)-own(:,:,o,p)))
%         title("o "+o+"p "+p)
%         drawnow;
%         pause(0.1)
    end
end

numTest = numTest + 1;

if ~problem
    disp('wn Test: Pass')
    passTest = passTest + 1;
else
    disp('wn Test: FAIL < ----------------------------------------')
end

clearvars -except xplot yplot m n passTest numTest

% ===========================================================================================================================================================================================================================
% -- checking if Cg is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_ref.mat','Cg'); ocg = Cg;
load('New_1_ref.mat','Cg'); ncg = Cg;

clearvars Cg

problem = 0;
for o = 1:25
    for p = 1:287
        A = ocg(:,:,o,p);
        B = ncg(:,:,o,p);
        Index{o,p} = find(A~=B);
        if ~isempty(Index{o,p})
            problem = problem + 1;
        end
    end
end

numTest = numTest + 1;

if ~problem
    disp('Cg Test: Pass')
    passTest = passTest + 1;
else
    disp('Cg Test: FAIL < ----------------------------------------')
end

clearvars -except xplot yplot m n passTest numTest

% ===========================================================================================================================================================================================================================
% -- checking if c is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_ref.mat','c'); oc = c;
load('New_1_ref.mat','c'); nc = c;

clearvars Cg

problem = 0;
for o = 1:25
    for p = 1:287
        A = oc(:,:,o,p);
        B = nc(:,:,o,p);
        Index{o,p} = find(A~=B);
        if ~isempty(Index{o,p})
            problem = problem + 1;
        end

    end
end

numTest = numTest + 1;

if ~problem
    disp('c Test: Pass')
    passTest = passTest + 1;
else
    disp('c Test: FAIL < ----------------------------------------')
end

clearvars -except xplot yplot m n passTest numTest

% ===========================================================================================================================================================================================================================
% -- checking if delx is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_ref.mat','delx'); odelx = delx;
load('New_1_ref.mat','delx'); ndelx = delx;

clearvars delx

problem = 0;
for o = 1:25
    for p = 1:287
        A = odelx(:,:,o,p);
        B = ndelx(:,:,o,p);
        Index{o,p} = find(A~=B);
        if ~isempty(Index{o,p})
            problem = problem + 1;
        end

    end
end

numTest = numTest + 1;

if ~problem
    disp('delx Test: Pass')
    passTest = passTest + 1;
else
    disp('delx Test: FAIL < ----------------------------------------')
end

clearvars -except xplot yplot m n passTest numTest

% ===========================================================================================================================================================================================================================
% -- checking if dely is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_ref.mat','dely'); odely = dely;
load('New_1_ref.mat','dely'); ndely = dely;

clearvars dely

problem = 0;
for o = 1:25
    for p = 1:287
        A = odely(:,:,o,p);
        B = ndely(:,:,o,p);
        Index{o,p} = find(A~=B);
        if ~isempty(Index{o,p})
            problem = problem + 1;
        end

    end
end

numTest = numTest + 1;

if ~problem
    disp('dely Test: Pass')
    passTest = passTest + 1;
else
    disp('dely Test: FAIL < ----------------------------------------')
end

clearvars -except xplot yplot m n passTest numTest

% ===========================================================================================================================================================================================================================
% -- check if f is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_ref.mat','f'); of = f;
load('New_1_ref.mat','f'); nf = f;

clearvars f

problem = 0;
for o = 1:25
    for p = 1:287
        A = of(:,:,o,p);
        B = nf(:,:,o,p);
        Index{o,p} = find(A~=B);
        if ~isempty(Index{o,p})
            problem = problem + 1;
        end

    end
end

numTest = numTest + 1;

if ~problem
    disp('f Test: Pass')
    passTest = passTest + 1;
else
    disp('f Test: FAIL < ----------------------------------------')
end

clearvars -except xplot yplot m n passTest numTest

% ===========================================================================================================================================================================================================================
% -- check if Uei is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_ref.mat','Uei'); oUei = Uei;
load('New_1_ref.mat','Uei'); nUei = Uei;

clearvars Uei

problem = 0;
for o = 1:25
    for p = 1:287
        A = oUei(:,:,o,p);
        B = nUei(:,:,o,p);
        Index{o,p} = find(A~=B);
        if ~isempty(Index{o,p})
            problem = problem + 1;
        end

    end
end

numTest = numTest + 1;

if ~problem
    disp('Uei Test: Pass')
    passTest = passTest + 1;
else
    disp('Uei Test: FAIL < ----------------------------------------')
end

clearvars -except xplot yplot m n passTest numTest

% ===========================================================================================================================================================================================================================
% -- check if Uer is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_ref.mat','Uer'); oUer = Uer;
load('New_1_ref.mat','Uer'); nUer = Uer;

clearvars dely

problem = 0;
for o = 1:25
    for p = 1:287
        A = oUer(:,:,o,p);
        B = nUer(:,:,o,p);
        Index{o,p} = find(A~=B);
        if ~isempty(Index{o,p})
            problem = problem + 1;
        end
    end
end

numTest = numTest + 1;

if ~problem
    disp('Uer Test: Pass')
    passTest = passTest + 1;
else
    disp('Uer Test: FAIL < ----------------------------------------')
end

clearvars -except xplot yplot m n passTest numTest

% ===========================================================================================================================================================================================================================
% -- check if waveang is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_ref.mat','waveang'); owaveang = waveang;
load('New_1_ref.mat','waveang'); nwaveang = waveang; 

clearvars waveang

problem = 0;
for o = 1:25
    for p = 1:287
        A = owaveang(:,:,o,p);
        B = nwaveang(:,:,o,p);
        Index{o,p} = find(A~=B);
        if ~isempty(Index{o,p})
            problem = problem + 1;
        end

    end
end

numTest = numTest + 1;


if ~problem
    disp('waveang Test: Pass')
    passTest = passTest + 1;
else
    disp('waveang Test: FAIL < ----------------------------------------')
end

clearvars -except xplot yplot m n passTest numTest
% ===========================================================================================================================================================================================================================
% ===========================================================================================================================================================================================================================

% report pass rate
passRate = (passTest/numTest)*100;
fprintf('%.2f percent of tests passed\n',passRate)
disp('testing_changes_ref.m completed')

% ===========================================================================================================================================================================================================================
% ===========================================================================================================================================================================================================================