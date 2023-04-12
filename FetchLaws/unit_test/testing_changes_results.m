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

% location of .mat files 
addpath('C:\Users\schne\OneDrive\Desktop\Main\Work\MIT\Github_Repos\umwm_titan\FetchLaws\unit_test\umwm_mat\umwm_old') % original model
addpath('C:\Users\schne\OneDrive\Desktop\Main\Work\MIT\Github_Repos\umwm_titan\FetchLaws\unit_test\umwm_mat\umwm_new') % new model with changes

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
% -- checking that E is the same
% ===========================================================================================================================================================================================================================
load('Old_1_1.mat','E'); oE = E;
load('New_1_1.mat','E'); nE = E;
clearvars E

problem = 0;
for o = 1:25
    for p = 1:287
        A = oE(:,:,o,p);
        B = nE(:,:,o,p);
        Index{o,p} = find(A~=B);
        if ~isempty(Index{o,p})
            problem = problem + 1;
        end
    end
end


if ~problem
    disp('E Test 1_1: Pass')
    passTest = passTest + 1;
else
    disp('E Test 1_1: FAIL < ----------------------------------------')
end

numTest = numTest + 1;

clearvars oE nE

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_1_2.mat','E'); oE = E;
load('New_1_2.mat','E'); nE = E;

clearvars E

problem = 0;
for o = 1:25
    for p = 1:287
        A = oE(:,:,o,p);
        B = nE(:,:,o,p);
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


if ~problem
    disp('E Test 1_2: Pass')
    passTest = passTest + 1;
else
    disp('E Test 1_2: FAIL < ----------------------------------------')
end

numTest = numTest + 1;

clearvars oE nE

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_1.mat','E'); oE = E;
load('New_2_1.mat','E'); nE = E;

clearvars E

problem = 0;
for o = 1:25
    for p = 1:287
        A = oE(:,:,o,p);
        B = nE(:,:,o,p);
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


if ~problem
    disp('E Test 2_1: Pass')
    passTest = passTest + 1;
else
    disp('E Test 2_1: FAIL < ----------------------------------------')
end

numTest = numTest + 1;

clearvars oE nE

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_2.mat','E'); oE = E;
load('New_2_2.mat','E'); nE = E;

clearvars E

problem = 0;
for o = 1:25
    for p = 1:287
        A = oE(:,:,o,p);
        B = nE(:,:,o,p);
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

if ~problem
    disp('E Test 2_2: Pass')
    passTest = passTest + 1;
else
    disp('E Test 2_2: FAIL < ----------------------------------------')
end

numTest = numTest + 1;

clearvars oE nE

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_1.mat','E'); oE = E;
load('New_3_1.mat','E'); nE = E;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oE(:,:,o,p);
        B = nE(:,:,o,p);
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



if ~problem
    disp('E Test 3_1: Pass')
    passTest = passTest + 1;
else
    disp('E Test 3_1: FAIL < ----------------------------------------')
end


numTest = numTest + 1;

clearvars oE nE
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_2.mat','E'); oE = E;
load('New_3_2.mat','E'); nE = E;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oE(:,:,o,p);
        B = nE(:,:,o,p);
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



if ~problem
    disp('E Test 3_2: Pass')
    passTest = passTest + 1;
else
    disp('E Test 3_2: FAIL < ----------------------------------------')
end

numTest = numTest + 1;

clearvars -except xplot yplot m n passTest numTest

% ===========================================================================================================================================================================================================================
% -- checking that ht is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_1.mat','ht'); oht = ht;
load('New_1_1.mat','ht'); nht = ht;

clearvars ht

problem = 0 ;

if nansum(nansum(oht - nht)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('ht test 1_1: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oht - nht)))
else
    disp('ht test 1_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars oht nht problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_1_2.mat','ht'); oht = ht;
load('New_1_2.mat','ht'); nht = ht;

clearvars ht

problem = 0 ;

if nansum(nansum(oht - nht)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('ht test 1_2: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n ',nansum(nansum(oht - nht)))
else
    disp('ht test 1_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars oht nht problem


% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_2_1.mat','ht'); oht = ht;
load('New_2_1.mat','ht'); nht = ht;

clearvars ht

problem = 0 ;

if nansum(nansum(oht - nht)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('ht test 2_1: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oht - nht)))
else
    disp('Cd test 2_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars oht nht problem


% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_2_2.mat','ht'); oht = ht;
load('New_2_2.mat','ht'); nht = ht;

clearvars ht

problem = 0 ;

if nansum(nansum(oht - nht)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('ht test 2_2: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n ',nansum(nansum(oht - nht)))
else
    disp('ht test 2_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;
clearvars oht nht problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_3_1.mat','ht'); oht = ht;
load('New_3_1.mat','ht'); nht = ht;

clearvars ht

problem = 0 ;

if nansum(nansum(oht - nht)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('ht test 3_1: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n ',nansum(nansum(oht - nht)))
else
    disp('ht test 3_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;
clearvars oht nht problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_3_2.mat','ht'); oht = ht;
load('New_3_2.mat','ht'); nht = ht;

clearvars ht

problem = 0 ;

if nansum(nansum(oht - nht)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('ht test 3_2: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n ',nansum(nansum(oht - nht)))
else
    disp('ht test 3_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

disp('ht tests done')

clearvars -except xplot yplot m n passTest numTest

% ===========================================================================================================================================================================================================================
% -- checking that freqs is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_1.mat','freqs'); ofreqs = freqs;
load('New_1_1.mat','freqs'); nfreqs = freqs;

clearvars freqs

if ~any(nfreqs == ofreqs)
    disp('freqs test 1_1: FAIL < ----------------------------------------')
else
    disp('freqs test 1_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars ofreqs nfreqs
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_1_2.mat','freqs'); ofreqs = freqs;
load('New_1_2.mat','freqs'); nfreqs = freqs;

clearvars freqs

if ~any(nfreqs == ofreqs)
    disp('freqs test 1_2: FAIL < ----------------------------------------')
else
    disp('freqs test 1_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars ofreqs nfreqs
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_1.mat','freqs'); ofreqs = freqs;
load('New_2_1.mat','freqs'); nfreqs = freqs;

clearvars freqs

if ~any(nfreqs == ofreqs)
    disp('freqs test 2_1: FAIL < ----------------------------------------')
else
    disp('freqs test 2_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;
clearvars ofreqs nfreqs
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_2.mat','freqs'); ofreqs = freqs;
load('New_2_2.mat','freqs'); nfreqs = freqs;

clearvars freqs

if ~any(nfreqs == ofreqs)
    disp('freqs test 2_2: FAIL < ----------------------------------------')
else
    disp('freqs test 2_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;
clearvars ofreqs nfreqs
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_1.mat','freqs'); ofreqs = freqs;
load('New_3_1.mat','freqs'); nfreqs = freqs;

clearvars freqs

if ~any(nfreqs == ofreqs)
    disp('freqs test 3_1: FAIL < ----------------------------------------')
else
    disp('freqs test 3_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars ofreqs nfreqs
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_2.mat','freqs'); ofreqs = freqs;
load('New_3_2.mat','freqs'); nfreqs = freqs;

clearvars freqs

if ~any(nfreqs == ofreqs)
    disp('freqs test 3_2: FAIL < ----------------------------------------')
else
    disp('freqs test 3_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

disp('freqs test done')
clearvars -except xplot yplot m n passTest numTest

% ===========================================================================================================================================================================================================================
% -- checking that Cd is the same -- %
% ===========================================================================================================================================================================================================================

load('Old_1_1.mat','Cd'); oCd = Cd;
load('New_1_1.mat','Cd'); nCd = Cd;

clearvars Cd

problem = 0 ;

if nansum(nansum(oCd - nCd)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cd test 1_1: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCd - nCd)))
else
    disp('Cd test 1_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;
clearvars oCd nCd problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_1_2.mat','Cd'); oCd = Cd;
load('New_1_2.mat','Cd'); nCd = Cd;

clearvars Cd

problem = 0 ;

if nansum(nansum(oCd - nCd)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cd test 1_2: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCd - nCd)))
else
    disp('Cd test 1_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;
clearvars oCd nCd problem


% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_2_1.mat','Cd'); oCd = Cd;
load('New_2_1.mat','Cd'); nCd = Cd;

clearvars Cd

problem = 0 ;

if nansum(nansum(oCd - nCd)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cd test 2_1: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCd - nCd)))
else
    disp('Cd test 2_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;
clearvars oCd nCd problem


% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_2_2.mat','Cd'); oCd = Cd;
load('New_2_2.mat','Cd'); nCd = Cd;

clearvars Cd

problem = 0 ;

if nansum(nansum(oCd - nCd)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cd test 2_2: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCd - nCd)))
else
    disp('Cd test 2_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;
clearvars oCd nCd problem
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_3_1.mat','Cd'); oCd = Cd;
load('New_3_1.mat','Cd'); nCd = Cd;

clearvars Cd

problem = 0 ;

if nansum(nansum(oCd - nCd)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cd test 3_1: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCd - nCd)))
else
    disp('Cd test 3_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;
clearvars oCd nCd problem


% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_3_2.mat','Cd'); oCd = Cd;
load('New_3_2.mat','Cd'); nCd = Cd;

clearvars Cd

problem = 0 ;

if nansum(nansum(oCd - nCd)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cd test 3_2: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCd - nCd)))
else
    disp('Cd test 3_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars -except xplot yplot m n passTest numTest


disp('Cd tests done')

% ===========================================================================================================================================================================================================================
% -- checking that Cdf is the same
% ===========================================================================================================================================================================================================================

load('Old_1_1.mat','Cdf'); oCdf = Cdf;
load('New_1_1.mat','Cdf'); nCdf = Cdf;

clearvars Cdf

problem = 0 ;

if nansum(nansum(oCdf - nCdf)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cdf test 1_1: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCdf - nCdf)))
else
    disp('Cdf test 1_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;
clearvars oCdf nCdf problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_1_2.mat','Cdf'); oCdf = Cdf;
load('New_1_2.mat','Cdf'); nCdf = Cdf;

clearvars Cdf

problem = 0 ;

if nansum(nansum(oCdf - nCdf)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cdf test 1_2: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCdf - nCdf)))
else
    disp('Cdf test 1_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;
clearvars oCdf nCdf problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_1.mat','Cdf'); oCdf = Cdf;
load('New_2_1.mat','Cdf'); nCdf = Cdf;

clearvars Cdf

problem = 0 ;

if nansum(nansum(oCdf - nCdf)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cdf test 2_1: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCdf - nCdf)))
else
    disp('Cdf test 2_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;
clearvars oCdf nCdf problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_2.mat','Cdf'); oCdf = Cdf;
load('New_2_2.mat','Cdf'); nCdf = Cdf;

clearvars Cdf

problem = 0 ;

if nansum(nansum(oCdf - nCdf)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cdf test 2_2: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCdf - nCdf)))
else
    disp('Cdf test 2_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars oCdf nCdf problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_1.mat','Cdf'); oCdf = Cdf;
load('New_3_1.mat','Cdf'); nCdf = Cdf;

clearvars Cdf

problem = 0 ;

if nansum(nansum(oCdf - nCdf)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cdf test 3_1: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCdf - nCdf)))
else
    disp('Cdf test 3_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars oCdf nCdf problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_2.mat','Cdf'); oCdf = Cdf;
load('New_3_2.mat','Cdf'); nCdf = Cdf;

clearvars Cdf

problem = 0 ;

if nansum(nansum(oCdf - nCdf)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cdf test 3_2: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCdf - nCdf)))
else
    disp('Cdf test 3_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars -except xplot yplot m n passTest numTest

disp('Cdf tests done')

% ===========================================================================================================================================================================================================================
% -- checking that Cds is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_1.mat','Cds'); oCds = Cds;
load('New_1_1.mat','Cds'); nCds = Cds;

clearvars Cds

problem = 0 ;

if nansum(nansum(oCds - nCds)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cds test 1_1: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCds - nCds)))
else
    disp('Cds test 1_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;
clearvars oCds nCds problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_1_2.mat','Cds'); oCds = Cds;
load('New_1_2.mat','Cds'); nCds = Cds;

clearvars Cds

problem = 0 ;

if nansum(nansum(oCds - nCds)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cds test 1_2: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCds - nCds)))
else
    disp('Cds test 1_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars oCds nCds problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_1.mat','Cds'); oCds = Cds;
load('New_2_1.mat','Cds'); nCds = Cds;

clearvars Cds

problem = 0 ;

if nansum(nansum(oCds - nCds)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cds test 2_1: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCds - nCds)))
else
    disp('Cds test 2_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars oCds nCds problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_2.mat','Cds'); oCds = Cds;
load('New_2_2.mat','Cds'); nCds = Cds;

clearvars Cds

problem = 0 ;

if nansum(nansum(oCds - nCds)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cds test 2_2: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCds - nCds)))
else
    disp('Cds test 2_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars oCds nCds problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_1.mat','Cds'); oCds = Cds;
load('New_3_1.mat','Cds'); nCds = Cds;

clearvars Cds

problem = 0 ;

if nansum(nansum(oCds - nCds)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cds test 3_1: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCds - nCds)))
else
    disp('Cds test 3_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars oCds nCds problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_2.mat','Cds'); oCds = Cds;
load('New_3_2.mat','Cds'); nCds = Cds;

clearvars Cds

problem = 0 ;

if nansum(nansum(oCds - nCds)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('Cds test 3_2: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oCds - nCds)))
else
    disp('Cds test 3_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars -except xplot yplot m n passTest numTest

disp('Cds tests done')

% ===========================================================================================================================================================================================================================
% -- checking that Sin is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_1.mat','Sin'); oSin = Sin;
load('New_1_1.mat','Sin'); nSin = Sin;

problem = 0;

for o = 1:25
    for p = 1:287
        A = oSin(:,:,o,p);
        B = nSin(:,:,o,p);
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


if ~problem
    disp('Sin Test 1_1: Pass')
    passTest = passTest + 1;
else
    disp('Sin Test 1_1: FAIL < ----------------------------------------')
end

numTest = numTest + 1;


clearvars oSin nSin
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_1_2.mat','Sin'); oSin = Sin;
load('New_1_2.mat','Sin'); nSin = Sin;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSin(:,:,o,p);
        B = nSin(:,:,o,p);
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


if ~problem
    disp('Sin Test 1_2: Pass')
    passTest = passTest + 1;
else
    disp('Sin Test 1_2: FAIL < ----------------------------------------')
end

numTest = numTest + 1;


clearvars oSin nSin
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_1.mat','Sin'); oSin = Sin;
load('New_2_1.mat','Sin'); nSin = Sin;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSin(:,:,o,p);
        B = nSin(:,:,o,p);
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


if ~problem
    disp('Sin Test 2_1: Pass')
    passTest = passTest + 1;
else
    disp('Sin Test 2_1: FAIL < ----------------------------------------')
end

numTest = numTest + 1;


clearvars oSin nSin
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_2.mat','Sin'); oSin = Sin;
load('New_2_2.mat','Sin'); nSin = Sin;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSin(:,:,o,p);
        B = nSin(:,:,o,p);
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


if ~problem
    disp('Sin Test 2_2: Pass')
    passTest = passTest + 1;
else
    disp('Sin Test 2_2: FAIL < ----------------------------------------')
end


numTest = numTest + 1;


clearvars oSin nSin
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_1.mat','Sin'); oSin = Sin;
load('New_3_1.mat','Sin'); nSin = Sin;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSin(:,:,o,p);
        B = nSin(:,:,o,p);
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


if ~problem
    disp('Sin Test 3_1: Pass')
    passTest = passTest + 1;
else
    disp('Sin Test 3_1: FAIL < ----------------------------------------')
end

numTest = numTest + 1;


clearvars oSin nSin
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_2.mat','Sin'); oSin = Sin;
load('New_3_2.mat','Sin'); nSin = Sin;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSin(:,:,o,p);
        B = nSin(:,:,o,p);
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


if ~problem
    disp('Sin Test 3_2: Pass')
    passTest = passTest + 1;
else
    disp('Sin Test 3_2: FAIL < ----------------------------------------')
end

numTest = numTest + 1;

clearvars -except xplot yplot m n passTest numTest

disp('Sin tests done')

% ===========================================================================================================================================================================================================================
% -- checking that Sds_wc is the same 
% ===========================================================================================================================================================================================================================
load('Old_1_1.mat','Sds_wc'); oSds_wc = Sds_wc;
load('New_1_1.mat','Sds_wc'); nSds_wc = Sds_wc;

problem = 0;

for o = 1:25
    for p = 1:287
        A = oSds_wc(:,:,o,p);
        B = nSds_wc(:,:,o,p);
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


if ~problem
    disp('Sds_wc Test 1_1: Pass')
    passTest = passTest + 1;
else
    disp('Sds_wc Test 1_1: FAIL < ----------------------------------------')
end

numTest = numTest + 1;

clearvars oSds_wc nSds_wc
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_1_2.mat','Sds_wc'); oSds_wc = Sds_wc;
load('New_1_2.mat','Sds_wc'); nSds_wc = Sds_wc;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSds_wc(:,:,o,p);
        B = nSds_wc(:,:,o,p);
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


if ~problem
    disp('Sds_wc Test 1_2: Pass')
    passTest = passTest + 1;
else
    disp('Sds_wc Test 1_2: FAIL < ----------------------------------------')
end

numTest = numTest + 1;


clearvars oSds_wc nSds_wc
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_1.mat','Sds_wc'); oSds_wc = Sds_wc;
load('New_2_1.mat','Sds_wc'); nSds_wc = Sds_wc;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSds_wc(:,:,o,p);
        B = nSds_wc(:,:,o,p);
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


if ~problem
    disp('Sds_wc Test 2_1: Pass')
    passTest = passTest + 1;
else
    disp('Sds_wc Test 2_1: FAIL < ----------------------------------------')
end

numTest = numTest + 1;


clearvars oSds_wc nSds_wc
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_2.mat','Sds_wc'); oSds_wc = Sds_wc;
load('New_2_2.mat','Sds_wc'); nSds_wc = Sds_wc;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSds_wc(:,:,o,p);
        B = nSds_wc(:,:,o,p);
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


if ~problem
    disp('Sds_wc Test 2_2: Pass')
    passTest = passTest + 1;
else
    disp('Sds_wc Test 2_2: FAIL < ----------------------------------------')
end

numTest = numTest + 1;


clearvars oSds_wc nSds_wc
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_1.mat','Sds_wc'); oSds_wc = Sds_wc;
load('New_3_1.mat','Sds_wc'); nSds_wc = Sds_wc;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSds_wc(:,:,o,p);
        B = nSds_wc(:,:,o,p);
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


if ~problem
    disp('Sds_wc Test 3_1: Pass')
    passTest = passTest + 1;
else
    disp('Sds_wc Test 3_1: FAIL < ----------------------------------------')
end
numTest = numTest + 1;


clearvars oSds_wc nSds_wc
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_2.mat','Sds_wc'); oSds_wc = Sds_wc;
load('New_3_2.mat','Sds_wc'); nSds_wc = Sds_wc;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSds_wc(:,:,o,p);
        B = nSds_wc(:,:,o,p);
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


if ~problem
    disp('Sds_wc Test 3_2: Pass')
    passTest = passTest + 1;
else
    disp('Sds_wc Test 3_2: FAIL < ----------------------------------------')
end
numTest = numTest + 1;

clearvars -except xplot yplot m n passTest numTest

disp('Sin tests done')

% ===========================================================================================================================================================================================================================
% -- checking that Sds is the same 
% ===========================================================================================================================================================================================================================
% ---------------------------------------------- %
load('Old_1_1.mat','Sds'); oSds = Sds;
load('New_1_1.mat','Sds'); nSds = Sds;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSds(:,:,o,p);
        B = nSds(:,:,o,p);
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


if ~problem
    disp('Sds Test 1_1: Pass')
    passTest = passTest + 1;
else
    disp('Sds Test 1_1: FAIL < ----------------------------------------')
end

numTest = numTest + 1;

clearvars oSds nSds
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_1_2.mat','Sds'); oSds = Sds;
load('New_1_2.mat','Sds'); nSds = Sds;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSds(:,:,o,p);
        B = nSds(:,:,o,p);
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


if ~problem
    disp('Sds Test 1_2: Pass')
    passTest = passTest + 1;
else
    disp('Sds Test 1_2: FAIL < ----------------------------------------')
end


numTest = numTest + 1;

clearvars oSds nSds
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_1.mat','Sds'); oSds = Sds;
load('New_2_1.mat','Sds'); nSds = Sds;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSds(:,:,o,p);
        B = nSds(:,:,o,p);
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


if ~problem
    disp('Sds Test 2_1: Pass')
    passTest = passTest + 1;
else
    disp('Sds Test 2_1: FAIL < ----------------------------------------')
end


numTest = numTest + 1;

clearvars oSds nSds
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_2.mat','Sds'); oSds = Sds;
load('New_2_2.mat','Sds'); nSds = Sds;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSds(:,:,o,p);
        B = nSds(:,:,o,p);
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


if ~problem
    disp('Sds Test 2_2: Pass')
    passTest = passTest + 1;
else
    disp('Sds Test 2_2: FAIL < ----------------------------------------')
end


numTest = numTest + 1;

clearvars oSds nSds
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_1.mat','Sds'); oSds = Sds;
load('New_3_1.mat','Sds'); nSds = Sds;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSds(:,:,o,p);
        B = nSds(:,:,o,p);
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


if ~problem
    disp('Sds Test 3_1: Pass')
    passTest = passTest + 1;
else
    disp('Sds Test 3_1: FAIL < ----------------------------------------')
end


numTest = numTest + 1;

clearvars oSds nSds
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_2.mat','Sds'); oSds = Sds;
load('New_3_2.mat','Sds'); nSds = Sds;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSds(:,:,o,p);
        B = nSds(:,:,o,p);
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


if ~problem
    disp('Sds Test 3_2: Pass')
    passTest = passTest + 1;
else
    disp('Sds Test 3_2: FAIL < ----------------------------------------')
end


numTest = numTest + 1;

clearvars -except xplot yplot m n passTest numTest

disp('Sds tests done')

% ===========================================================================================================================================================================================================================
% -- checking that Snl is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_1.mat','Snl'); oSnl = Snl;
load('New_1_1.mat','Snl'); nSnl = Snl;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSnl(:,:,o,p);
        B = nSnl(:,:,o,p);
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


if ~problem
    disp('Snl Test 1_1: Pass')
    passTest = passTest + 1;
else
    disp('Snl Test 1_1: FAIL < ----------------------------------------')
end


numTest = numTest + 1;

clearvars oSnl nSnl
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_1_2.mat','Snl'); oSnl = Snl;
load('New_1_2.mat','Snl'); nSnl = Snl;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSnl(:,:,o,p);
        B = nSnl(:,:,o,p);
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


if ~problem
    disp('Snl Test 1_2: Pass')
    passTest = passTest + 1;
else
    disp('Snl Test 1_2: FAIL < ----------------------------------------')
end

numTest = numTest + 1;

clearvars oSnl nSnl
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_1.mat','Snl'); oSnl = Snl;
load('New_2_1.mat','Snl'); nSnl = Snl;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSnl(:,:,o,p);
        B = nSnl(:,:,o,p);
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


if ~problem
    disp('Snl Test 2_1: Pass')
    passTest = passTest + 1;
else
    disp('Snl Test 2_1: FAIL < ----------------------------------------')
end

numTest = numTest + 1;

clearvars oSnl nSnl
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_2.mat','Snl'); oSnl = Snl;
load('New_2_2.mat','Snl'); nSnl = Snl;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSnl(:,:,o,p);
        B = nSnl(:,:,o,p);
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


if ~problem
    disp('Snl Test 2_2: Pass')
    passTest = passTest + 1;
else
    disp('Snl Test 2_2: FAIL < ----------------------------------------')
end

numTest = numTest + 1;

clearvars oSnl nSnl
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_1.mat','Snl'); oSnl = Snl;
load('New_3_1.mat','Snl'); nSnl = Snl;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSnl(:,:,o,p);
        B = nSnl(:,:,o,p);
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


if ~problem
    disp('Snl Test 3_1: Pass')
    passTest = passTest + 1;
else
    disp('Snl Test 3_1: FAIL < ----------------------------------------')
end

numTest = numTest + 1;

clearvars oSnl nSnl
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_2.mat','Snl'); oSnl = Snl;
load('New_3_2.mat','Snl'); nSnl = Snl;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSnl(:,:,o,p);
        B = nSnl(:,:,o,p);
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


if ~problem
    disp('Snl Test 3_2: Pass')
    passTest = passTest + 1;
else
    disp('Snl Test 3_2: FAIL < ----------------------------------------')
end
numTest = numTest + 1;

clearvars -except xplot yplot m n passTest numTest

disp('Snl tests done')

% ===========================================================================================================================================================================================================================
% -- checking that Sdt is the same 
% ===========================================================================================================================================================================================================================
load('Old_1_1.mat','Sdt'); oSdt = Sdt;
load('New_1_1.mat','Sdt'); nSdt = Sdt;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSdt(:,:,o,p);
        B = nSdt(:,:,o,p);
        if ~isequaln(A,B)
            problem = problem + 1;
        end
        %surf(xplot',yplot',abs(nD(:,:,o,p)-oD(:,:,o,p)))
        %title("o "+o+"p "+p)
        %drawnow;
        %pause(0.1)
    end
end


if ~problem
    disp('Sdt Test 1_1: Pass')
    passTest = passTest + 1;
else
    disp('Sdt Test 1_1: FAIL < ----------------------------------------')
end

numTest = numTest + 1;

clearvars oSdt nSdt
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_1_2.mat','Sdt'); oSdt = Sdt;
load('New_1_2.mat','Sdt'); nSdt = Sdt;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSdt(:,:,o,p);
        B = nSdt(:,:,o,p);
        if ~isequaln(A,B)
            problem = problem + 1;
        end
        %surf(xplot',yplot',abs(nD(:,:,o,p)-oD(:,:,o,p)))
        %title("o "+o+"p "+p)
        %drawnow;
        %pause(0.1)
    end
end


if ~problem
    disp('Sdt Test 1_2: Pass')
    passTest = passTest + 1;
else
    disp('Sdt Test 1_2: FAIL < ----------------------------------------')
end


numTest = numTest + 1;

clearvars oSdt nSdt
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_1.mat','Sdt'); oSdt = Sdt;
load('New_2_1.mat','Sdt'); nSdt = Sdt;


problem = 0;
for o = 1:25
    for p = 1:287
        A = oSdt(:,:,o,p);
        B = nSdt(:,:,o,p);
        if ~isequaln(A,B)
            problem = problem + 1;
        end
        %surf(xplot',yplot',abs(nD(:,:,o,p)-oD(:,:,o,p)))
        %title("o "+o+"p "+p)
        %drawnow;
        %pause(0.1)
    end
end


if ~problem
    disp('Sdt Test 2_1: Pass')
    passTest = passTest + 1;
else
    disp('Sdt Test 2_1: FAIL < ----------------------------------------')
end

numTest = numTest + 1;

clearvars oSdt nSdt
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_2.mat','Sdt'); oSdt = Sdt;
load('New_2_2.mat','Sdt'); nSdt = Sdt;


problem = 0;
for o = 1:25
    for p = 1:287
        A = oSdt(:,:,o,p);
        B = nSdt(:,:,o,p);
        if ~isequaln(A,B)
            problem = problem + 1;
        end
        %surf(xplot',yplot',abs(nD(:,:,o,p)-oD(:,:,o,p)))
        %title("o "+o+"p "+p)
        %drawnow;
        %pause(0.1)
    end
end


if ~problem
    disp('Sdt Test 2_2: Pass')
    passTest = passTest + 1;
else
    disp('Sdt Test 2_2: FAIL < ----------------------------------------')
end

numTest = numTest + 1;

clearvars oSdt nSdt
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_1.mat','Sdt'); oSdt = Sdt;
load('New_3_1.mat','Sdt'); nSdt = Sdt;


problem = 0;
for o = 1:25
    for p = 1:287
        A = oSdt(:,:,o,p);
        B = nSdt(:,:,o,p);
        if ~isequaln(A,B)
            problem = problem + 1;
        end
        %surf(xplot',yplot',abs(nD(:,:,o,p)-oD(:,:,o,p)))
        %title("o "+o+"p "+p)
        %drawnow;
        %pause(0.1)
    end
end

if ~problem
    disp('Sdt Test 3_1: Pass')
    passTest = passTest + 1;
else
    disp('Sdt Test 3_1: FAIL < ----------------------------------------')
end


numTest = numTest + 1;

clearvars oSdt nSdt
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_2.mat','Sdt'); oSdt = Sdt;
load('New_3_2.mat','Sdt'); nSdt = Sdt;


problem = 0;
for o = 1:25
    for p = 1:287
        A = oSdt(:,:,o,p);
        B = nSdt(:,:,o,p);
        if ~isequaln(A,B)
            problem = problem + 1;
        end
        %surf(xplot',yplot',abs(nD(:,:,o,p)-oD(:,:,o,p)))
        %title("o "+o+"p "+p)
        %drawnow;
        %pause(0.1)
    end
end


if ~problem
    disp('Sdt Test 3_2: Pass')
    passTest = passTest + 1;
else
    disp('Sdt Test 3_2: FAIL < ----------------------------------------')
end
numTest = numTest + 1;


clearvars -except xplot yplot m n passTest numTest

disp('Sdt tests done')


% ===========================================================================================================================================================================================================================
% -- checking that Sbf is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_1.mat','Sbf'); oSbf = Sbf;
load('New_1_1.mat','Sbf'); nSbf = Sbf;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSbf(:,:,o,p);
        B = nSbf(:,:,o,p);
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


if ~problem
    disp('Sbf Test 1_1: Pass')
    passTest = passTest + 1;
else
    disp('Sbf Test 1_1: FAIL < ----------------------------------------')
end
numTest = numTest + 1;
clearvars oSbf nSbf
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_1_2.mat','Sbf'); oSbf = Sbf;
load('New_1_2.mat','Sbf'); nSbf = Sbf;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSbf(:,:,o,p);
        B = nSbf(:,:,o,p);
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


if ~problem
    disp('Sbf Test 1_2: Pass')
    passTest = passTest + 1;
else
    disp('Sbf Test 1_2: FAIL < ----------------------------------------')
end
numTest = numTest + 1;
clearvars oSbf nSbf

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_2_1.mat','Sbf'); oSbf = Sbf;
load('New_2_1.mat','Sbf'); nSbf = Sbf;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSbf(:,:,o,p);
        B = nSbf(:,:,o,p);
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


if ~problem
    disp('Sbf Test 2_1: Pass')
    passTest = passTest + 1;
else
    disp('Sbf Test 2_1: FAIL < ----------------------------------------')
end
numTest = numTest + 1;
clearvars oSbf nSbf
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_2_2.mat','Sbf'); oSbf = Sbf;
load('New_2_2.mat','Sbf'); nSbf = Sbf;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSbf(:,:,o,p);
        B = nSbf(:,:,o,p);
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


if ~problem
    disp('Sbf Test 2_2: Pass')
    passTest = passTest + 1;
else
    disp('Sbf Test 2_2: FAIL < ----------------------------------------')
end
numTest = numTest + 1;
clearvars oSbf nSbf
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_1.mat','Sbf'); oSbf = Sbf;
load('New_3_1.mat','Sbf'); nSbf = Sbf;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSbf(:,:,o,p);
        B = nSbf(:,:,o,p);
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


if ~problem
    disp('Sbf Test 3_1: Pass')
    passTest = passTest + 1;
else
    disp('Sbf Test 3_1: FAIL < ----------------------------------------')
end
numTest = numTest + 1;
clearvars oSbf nSbf
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Old_3_2.mat','Sbf'); oSbf = Sbf;
load('New_3_2.mat','Sbf'); nSbf = Sbf;

problem = 0;
for o = 1:25
    for p = 1:287
        A = oSbf(:,:,o,p);
        B = nSbf(:,:,o,p);
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


if ~problem
    disp('Sbf Test 3_2: Pass')
    passTest = passTest + 1;
else
    disp('Sbf Test 3_2: FAIL < ----------------------------------------')
end
numTest = numTest + 1;

clearvars -except xplot yplot m n passTest numTest

disp('Sbf tests done')

% ===========================================================================================================================================================================================================================
% -- checking that ms is the same 
% ===========================================================================================================================================================================================================================

load('Old_1_1.mat','ms'); oms = ms;
load('New_1_1.mat','ms'); nms = ms;

clearvars ms

problem = 0 ;

if nansum(nansum(oms - nms)) ~= 0 
    problem = problem + 1;SigH
end

if problem ~= 0
    disp('ms test 1_1: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oms - nms)))
else
    disp('ms test 1_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars oms nms problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_1_2.mat','ms'); oms = ms;
load('New_1_2.mat','ms'); nms = ms;

clearvars ms

problem = 0 ;

if nansum(nansum(oms - nms)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('ms test 1_2: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oms - nms)))
else
    disp('ms test 1_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars oms nms problem


% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_2_1.mat','ms'); oms = ms;
load('New_2_1.mat','ms'); nms = ms;

clearvars ms

problem = 0 ;

if nansum(nansum(oms - nms)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('ms test 2_1: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oms - nms)))
else
    disp('Cd test 2_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars oms nms problem


% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_2_2.mat','ms'); oms = ms;
load('New_2_2.mat','ms'); nms = ms;

clearvars ms

problem = 0 ;

if nansum(nansum(oms - nms)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('ms test 2_2: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oms - nms)))
else
    disp('ms test 2_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars oms nms problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_3_1.mat','ms'); oms = ms;
load('New_3_1.mat','ms'); nms = ms;

clearvars ms

problem = 0 ;

if nansum(nansum(oms - nms)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('ms test 3_1: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oms - nms)))
else
    disp('ms test 3_1: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;

clearvars oms nms problem


% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

load('Old_3_2.mat','ms'); oms = ms;
load('New_3_2.mat','ms'); nms = ms;

clearvars ms

problem = 0 ;

if nansum(nansum(oms - nms)) ~= 0 
    problem = problem + 1;
end

if problem ~= 0
    disp('ms test 3_2: FAIL < ----------------------------------------')
    fprintf('difference total: %.4f\n',nansum(nansum(oms - nms)))
else
    disp('ms test 3_2: PASS')
    passTest = passTest + 1;
end

numTest = numTest + 1;


clearvars -except xplot yplot m n passTest numTest

disp('ms tests done')

% ===========================================================================================================================================================================================================================
% ===========================================================================================================================================================================================================================

% report pass rate 
passRate = (passTest/numTest)*100;
fprintf('%.2f percent of tests passed\n',passRate)
disp('testing_changes_results.m completed')

% ===========================================================================================================================================================================================================================
% ===========================================================================================================================================================================================================================