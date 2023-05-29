% clc
% clear
% close all
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
%% Test
% m: number of gridpoints in x-direction 
load('Old_1_ref.mat','m'); om = m;
load('New_Reference.mat','m'); nm = m;
clearvars m
assert(isequal(om,nm),'m-dimension for original model has changed. Need to re-run original model for unit testing')

%% Test
% n: number of gridpoints in y-direction
load('Old_1_ref.mat','n'); on = n;
load('New_Reference.mat','n'); nn = n;
clearvars n
assert(isequal(on,nn),'n-dimension for original model has changed. Need to re-run original model for unit testing')

% ===========================================================================================================================================================================================================================
% -- checking that E is the same
% ===========================================================================================================================================================================================================================
%% Test
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

assert(problem==0,'E Test 1_1: FAIL')

clearvars oE nE

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'E Test 1_2: FAIL')

clearvars oE nE

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'E Test 2_1: FAIL')

clearvars oE nE

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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

assert(problem==0,'E Test 2_2: FAIL')

clearvars oE nE

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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



assert(problem==0,'E Test 3_1: FAIL')


clearvars oE nE
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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



assert(problem==0,'E Test 3_2: FAIL')


% ===========================================================================================================================================================================================================================
% -- checking that ht is the same 
% ===========================================================================================================================================================================================================================
%% Test
load('Old_1_1.mat','ht'); oht = ht;
load('New_1_1.mat','ht'); nht = ht;

clearvars ht

problem = 0 ;

if nansum(nansum(oht - nht)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'ht test 1_1: FAIL')

clearvars oht nht problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_1_2.mat','ht'); oht = ht;
load('New_1_2.mat','ht'); nht = ht;

clearvars ht

problem = 0 ;

if nansum(nansum(oht - nht)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'ht test 1_2: FAIL')

clearvars oht nht problem


% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_2_1.mat','ht'); oht = ht;
load('New_2_1.mat','ht'); nht = ht;

clearvars ht

problem = 0 ;

if nansum(nansum(oht - nht)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'ht test 2_1: FAIL')


clearvars oht nht problem


% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_2_2.mat','ht'); oht = ht;
load('New_2_2.mat','ht'); nht = ht;

clearvars ht

problem = 0 ;

if nansum(nansum(oht - nht)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'ht test 2_2: FAIL')


clearvars oht nht problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_3_1.mat','ht'); oht = ht;
load('New_3_1.mat','ht'); nht = ht;

clearvars ht

problem = 0 ;

if nansum(nansum(oht - nht)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'ht test 3_1: FAIL')

clearvars oht nht problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_3_2.mat','ht'); oht = ht;
load('New_3_2.mat','ht'); nht = ht;

clearvars ht

problem = 0 ;

if nansum(nansum(oht - nht)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'ht test 3_2: FAIL')


% ===========================================================================================================================================================================================================================
% -- checking that freqs is the same 
% ===========================================================================================================================================================================================================================
%% Test
load('Old_1_1.mat','freqs'); ofreqs = freqs;
load('New_1_1.mat','freqs'); nfreqs = freqs;

clearvars freqs

assert(isequal(nfreqs,ofreqs),'freqs test 1_1: FAIL')

clearvars ofreqs nfreqs
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_1_2.mat','freqs'); ofreqs = freqs;
load('New_1_2.mat','freqs'); nfreqs = freqs;

clearvars freqs

assert(isequal(nfreqs,ofreqs),'freqs test 1_2: FAIL')


clearvars ofreqs nfreqs
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_2_1.mat','freqs'); ofreqs = freqs;
load('New_2_1.mat','freqs'); nfreqs = freqs;

clearvars freqs

assert(isequal(nfreqs,ofreqs),'freqs test 2_1: FAIL')

clearvars ofreqs nfreqs
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_2_2.mat','freqs'); ofreqs = freqs;
load('New_2_2.mat','freqs'); nfreqs = freqs;

clearvars freqs

assert(isequal(nfreqs,ofreqs),'freqs test 2_2: FAIL')

clearvars ofreqs nfreqs
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_3_1.mat','freqs'); ofreqs = freqs;
load('New_3_1.mat','freqs'); nfreqs = freqs;

clearvars freqs

assert(isequal(nfreqs,ofreqs),'freqs test 3_1: FAIL')


clearvars ofreqs nfreqs
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_3_2.mat','freqs'); ofreqs = freqs;
load('New_3_2.mat','freqs'); nfreqs = freqs;

clearvars freqs

assert(isequal(nfreqs,ofreqs),'freqs test 3_2: FAIL')


% ===========================================================================================================================================================================================================================
% -- checking that Cd is the same -- %
% ===========================================================================================================================================================================================================================
%% Test
load('Old_1_1.mat','Cd'); oCd = Cd;
load('New_1_1.mat','Cd'); nCd = Cd;

clearvars Cd

problem = 0 ;

if nansum(nansum(oCd - nCd)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cd test 1_1: FAIL')

clearvars oCd nCd problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_1_2.mat','Cd'); oCd = Cd;
load('New_1_2.mat','Cd'); nCd = Cd;

clearvars Cd

problem = 0 ;

if nansum(nansum(oCd - nCd)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cd test 1_2: FAIL')


clearvars oCd nCd problem


% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_2_1.mat','Cd'); oCd = Cd;
load('New_2_1.mat','Cd'); nCd = Cd;

clearvars Cd

problem = 0 ;

if nansum(nansum(oCd - nCd)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cd test 2_1: FAIL')

clearvars oCd nCd problem


% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_2_2.mat','Cd'); oCd = Cd;
load('New_2_2.mat','Cd'); nCd = Cd;

clearvars Cd

problem = 0 ;

if nansum(nansum(oCd - nCd)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cd test 2_2: FAIL')

clearvars oCd nCd problem
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_3_1.mat','Cd'); oCd = Cd;
load('New_3_1.mat','Cd'); nCd = Cd;

clearvars Cd

problem = 0 ;

if nansum(nansum(oCd - nCd)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cd test 3_1: FAIL')

clearvars oCd nCd problem


% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_3_2.mat','Cd'); oCd = Cd;
load('New_3_2.mat','Cd'); nCd = Cd;

clearvars Cd

problem = 0 ;

if nansum(nansum(oCd - nCd)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cd test 3_2: FAIL')


% ===========================================================================================================================================================================================================================
% -- checking that Cdf is the same
% ===========================================================================================================================================================================================================================
%% Test
load('Old_1_1.mat','Cdf'); oCdf = Cdf;
load('New_1_1.mat','Cdf'); nCdf = Cdf;

clearvars Cdf

problem = 0 ;

if nansum(nansum(oCdf - nCdf)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cdf test 1_1: FAIL')

clearvars oCdf nCdf problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_1_2.mat','Cdf'); oCdf = Cdf;
load('New_1_2.mat','Cdf'); nCdf = Cdf;

clearvars Cdf

problem = 0 ;

if nansum(nansum(oCdf - nCdf)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cdf test 1_2: FAIL')

clearvars oCdf nCdf problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_2_1.mat','Cdf'); oCdf = Cdf;
load('New_2_1.mat','Cdf'); nCdf = Cdf;

clearvars Cdf

problem = 0 ;

if nansum(nansum(oCdf - nCdf)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cdf test 2_1: FAIL')

clearvars oCdf nCdf problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_2_2.mat','Cdf'); oCdf = Cdf;
load('New_2_2.mat','Cdf'); nCdf = Cdf;

clearvars Cdf

problem = 0 ;

if nansum(nansum(oCdf - nCdf)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cdf test 2_2: FAIL')


clearvars oCdf nCdf problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_3_1.mat','Cdf'); oCdf = Cdf;
load('New_3_1.mat','Cdf'); nCdf = Cdf;

clearvars Cdf

problem = 0 ;

if nansum(nansum(oCdf - nCdf)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cdf test 3_1: FAIL')

clearvars oCdf nCdf problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_3_2.mat','Cdf'); oCdf = Cdf;
load('New_3_2.mat','Cdf'); nCdf = Cdf;

clearvars Cdf

problem = 0 ;

if nansum(nansum(oCdf - nCdf)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cdf test 1_2: FAIL')

% ===========================================================================================================================================================================================================================
% -- checking that Cds is the same 
% ===========================================================================================================================================================================================================================
%% Test
load('Old_1_1.mat','Cds'); oCds = Cds;
load('New_1_1.mat','Cds'); nCds = Cds;

clearvars Cds

problem = 0 ;

if nansum(nansum(oCds - nCds)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cds test 1_1: FAIL')

clearvars oCds nCds problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_1_2.mat','Cds'); oCds = Cds;
load('New_1_2.mat','Cds'); nCds = Cds;

clearvars Cds

problem = 0 ;

if nansum(nansum(oCds - nCds)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cds test 1_2: FAIL')


clearvars oCds nCds problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_2_1.mat','Cds'); oCds = Cds;
load('New_2_1.mat','Cds'); nCds = Cds;

clearvars Cds

problem = 0 ;

if nansum(nansum(oCds - nCds)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cds test 2_1: FAIL')


clearvars oCds nCds problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_2_2.mat','Cds'); oCds = Cds;
load('New_2_2.mat','Cds'); nCds = Cds;

clearvars Cds

problem = 0 ;

if nansum(nansum(oCds - nCds)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cds test 2_2: FAIL')


clearvars oCds nCds problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_3_1.mat','Cds'); oCds = Cds;
load('New_3_1.mat','Cds'); nCds = Cds;

clearvars Cds

problem = 0 ;

if nansum(nansum(oCds - nCds)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cds test 3_1: FAIL')


clearvars oCds nCds problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_3_2.mat','Cds'); oCds = Cds;
load('New_3_2.mat','Cds'); nCds = Cds;

clearvars Cds

problem = 0 ;

if nansum(nansum(oCds - nCds)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'Cds test 3_2: FAIL')

% ===========================================================================================================================================================================================================================
% -- checking that Sin is the same 
% ===========================================================================================================================================================================================================================
%% Test
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


assert(problem==0,'Sin Test 1_1: FAI')

clearvars oSin nSin
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sin Test 1_2: FAI')

clearvars oSin nSin
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sin Test 2_1: FAI')


clearvars oSin nSin
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sin Test 2_2: FAI')

clearvars oSin nSin
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sin Test 3_1: FAI')



clearvars oSin nSin
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sin Test 3_2: FAI')





% ===========================================================================================================================================================================================================================
% -- checking that Sds_wc is the same 
% ===========================================================================================================================================================================================================================
%% Test
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


assert(problem==0,'Sds_wc Test 1_1: FAIL')


clearvars oSds_wc nSds_wc
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sds_wc Test 1_2: FAIL')


clearvars oSds_wc nSds_wc
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sds_wc Test 2_1: FAIL')



clearvars oSds_wc nSds_wc
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sds_wc Test 2_2: FAIL')



clearvars oSds_wc nSds_wc
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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

assert(problem==0,'Sds_wc Test 3_1: FAIL')



clearvars oSds_wc nSds_wc
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sds_wc Test 3_2: FAIL')




% ===========================================================================================================================================================================================================================
% -- checking that Sds is the same 
% ===========================================================================================================================================================================================================================
% ---------------------------------------------- %
%% Test
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

assert(problem==0,'Sds Test 1_1: FAIL')


clearvars oSds nSds
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sds Test 1_2: FAIL')


clearvars oSds nSds
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sds Test 2_1: FAIL')


clearvars oSds nSds
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sds Test 2_2: FAIL')

clearvars oSds nSds
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sds Test 3_1: FAIL')


clearvars oSds nSds
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sds Test 3_2: FAIL')

% ===========================================================================================================================================================================================================================
% -- checking that Snl is the same 
% ===========================================================================================================================================================================================================================
%% Test
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


assert(problem==0,'Snl Test 1_1: FAIL')


clearvars oSnl nSnl
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Snl Test 1_2: FAIL')


clearvars oSnl nSnl
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Snl Test 2_1: FAIL')


clearvars oSnl nSnl
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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

assert(problem==0,'Snl Test 2_2: FAIL')


clearvars oSnl nSnl
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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

assert(problem==0,'Snl Test 3_1: FAIL')


clearvars oSnl nSnl
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Snl Test 3_2: FAIL')


% ===========================================================================================================================================================================================================================
% -- checking that Sdt is the same 
% ===========================================================================================================================================================================================================================
%% Test
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


assert(problem==0,'Sdt Test 1_1: FAIL')



clearvars oSdt nSdt
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sdt Test 1_2: FAIL')


clearvars oSdt nSdt
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sdt Test 2_1: FAIL')


clearvars oSdt nSdt
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sdt Test 2_2: FAIL')


clearvars oSdt nSdt
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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

assert(problem==0,'Sdt Test 3_1: FAIL')


clearvars oSdt nSdt
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sdt Test 3_2: FAIL')





% ===========================================================================================================================================================================================================================
% -- checking that Sbf is the same 
% ===========================================================================================================================================================================================================================
%% Test
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


assert(problem==0,'Sbf Test 1_1: FAIL')


clearvars oSbf nSbf
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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

assert(problem==0,'Sbf Test 1_2: FAIL')

clearvars oSbf nSbf

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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

assert(problem==0,'Sbf Test 2_1: FAIL')

clearvars oSbf nSbf
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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

assert(problem==0,'Sbf Test 2_2: FAIL')

clearvars oSbf nSbf
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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

assert(problem==0,'Sbf Test 3_1: FAIL')

clearvars oSbf nSbf
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
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


assert(problem==0,'Sbf Test 3_2: FAIL')

% ===========================================================================================================================================================================================================================
% -- checking that ms is the same 
% ===========================================================================================================================================================================================================================
%% Test
load('Old_1_1.mat','ms'); oms = ms;
load('New_1_1.mat','ms'); nms = ms;

clearvars ms

problem = 0 ;

if nansum(nansum(oms - nms)) ~= 0 
    problem = problem + 1;SigH
end

assert(problem==0,'ms test 1_1: FAIL')



clearvars oms nms problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_1_2.mat','ms'); oms = ms;
load('New_1_2.mat','ms'); nms = ms;

clearvars ms

problem = 0 ;

if nansum(nansum(oms - nms)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'ms test 1_2: FAIL')


clearvars oms nms problem


% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_2_1.mat','ms'); oms = ms;
load('New_2_1.mat','ms'); nms = ms;

clearvars ms

problem = 0 ;

if nansum(nansum(oms - nms)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'ms test 2_1: FAIL')


clearvars oms nms problem


% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_2_2.mat','ms'); oms = ms;
load('New_2_2.mat','ms'); nms = ms;

clearvars ms

problem = 0 ;

if nansum(nansum(oms - nms)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'ms test 2_2: FAIL')


clearvars oms nms problem

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_3_1.mat','ms'); oms = ms;
load('New_3_1.mat','ms'); nms = ms;

clearvars ms

problem = 0 ;

if nansum(nansum(oms - nms)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'ms test 3_1: FAIL')

clearvars oms nms problem


% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%% Test
load('Old_3_2.mat','ms'); oms = ms;
load('New_3_2.mat','ms'); nms = ms;

clearvars ms

problem = 0 ;

if nansum(nansum(oms - nms)) ~= 0 
    problem = problem + 1;
end

assert(problem==0,'ms test 3_2: FAIL')

disp('ms tests done')

% ===========================================================================================================================================================================================================================
% ===========================================================================================================================================================================================================================

