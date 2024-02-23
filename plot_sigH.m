function [sigH] = plot_sigH(folderName,num_winds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Makes plots of sig wave height from saved prev runs saved in '\Titan'
% EX.) H = plot_sigH('cutoff_5') where cutoff_5 is folder in current directory
% INPUT:
%       folderName    = e.g. 'cutoff_1'
%       
%       optional:
%           num_winds = number of winds to do (defaults to all)
% OUTPUT:
%       sigH  (mxn)   = signifigant wave height at center of grid for each speed (m) and timestep (n) 
%
% ALSO MAKES A PLOT OF SIGNIFIGANT WAVE HEIGHT VERSUS TIME STEP
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

subfolder = '\Results';
ignoreFileName = 'New_Reference.mat';

% move to file with logfile
originalDirectory = pwd;
currentDirectory = strcat(pwd,'\',folderName);
cd(currentDirectory)

% get the legend from the logfile
legend_titles = get_log_winds();

% check if number of winds to plot is specified, defaults to all otherwise
if ~exist('num_winds','var')
    num_winds = length(legend_titles);
end

for qq = 1:num_winds
    fileList{qq} = strcat('New_',num2str(qq));
end

% move to file with data
currentDirectory = strcat(pwd,subfolder);
cd(currentDirectory);

% strip filename for figure title
folderName = strrep(folderName, '_', ' ');

sortedFilenames = sort_files();

% initialize array
sigH = NaN(numel(fileList),numel(sortedFilenames));
k = ones(1, length(fileList));

% Loop through each file
for ii = 1:numel(sortedFilenames)
    if strcmp(sortedFilenames{ii}, ignoreFileName)
        fprintf('Ignoring file: %s\n', ignoreFileName);
        continue; 
    end

    for jj = 1:numel(fileList)
    
        % Load the .mat file
        matData = load(fullfile(currentDirectory, sortedFilenames{ii}));
    
        % Check if 'ht' variable exists in the loaded data
        if isfield(matData, 'ht')
            % Access the 'ht' variable
            ht_loaded = matData.ht;
            [xLen,yLen] = size(ht_loaded);
    
            if startsWith(sortedFilenames{ii}, fileList{jj})
                sigH(jj, k(jj)) = ht_loaded(round(xLen/2), round(yLen/2));
                k(jj) = k(jj) + 1;
                fprintf('Loaded ''ht'' from %s\n', sortedFilenames{ii});
                break;  % Exit the inner loop 
            end
        else
            error('''ht'' variable not found in %s\n', sortedFilenames{ii});
        end
    end
    
end

make_figure()

cd(originalDirectory)

% sub-routines
    function sortedFilenames = sort_files()
        % files to load
        files = dir(fullfile(currentDirectory, '*.mat'));
        filenames = {files.name};
        % Sort the filenames using natural numbers (i.e. so A_3_4.mat comes before A_3_30.mat)
        sortedFilenames = natsort(filenames);
    end
    function make_figure()
        figure('units','normalized','outerposition',[0 0 1 1]);
        for pp = 1:length(fileList)
            plot(1:numel(sortedFilenames),sigH(pp,:),'-*','LineWidth',2)
            hold on;
        end
        legend(legend_titles,'location','best')
        grid on;
        xlabel('\Deltat')
        ylabel('H_{sig}')
        title(folderName)
        moveupdir = fullfile(pwd,'..','..');
        outputfolder = fullfile(moveupdir,'output_figs');
        if ~exist(outputfolder,'dir')
            mkdir(outputfolder)
        end
        savepath = fullfile(outputfolder,strcat(strrep(folderName, ' ', ''),'.png'));
        saveas(gcf,savepath);
    end

    function legend_titles = get_log_winds()
    % Read the entire log file 
    files = dir('*.txt');
    if isempty(files)
        disp('No .txt files found in the current directory.');
    else
        % Read the first .txt file
        fileName = files(1).name;
        filePath = fullfile(pwd, fileName);  % Full path to the file
        logContent = fileread(filePath);
    end
    
    % extract all numeric values after 'speed:' until \n
    matches = regexp(logContent, 'speed:\s*([\d\.\s,]+)\n', 'tokens', 'once');
    
    % Check if there is a match
    if ~isempty(matches)
        % Extracted values as a string
        extractedValuesString = strtrim(matches{1});  % Trim leading/trailing whitespace
        % Convert the string to a numeric array
        windSpeeds_legend = str2num(extractedValuesString);
    else
        error('Could not extract wind speeds from logfile\n');
    end
    
    % put wind speeds from logfile into legend
    for ii = 1:length(windSpeeds_legend)
        legend_titles{ii} = sprintf('u = %.1f m/s',windSpeeds_legend(ii));
    end
    
    end

end