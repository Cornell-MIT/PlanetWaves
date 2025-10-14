function test_LakeSuperior()
% Unit test to check if the results of makeWaves.m has changed from baseline established with UMWM-fortran version
    logDir = 'planetwaves/unit_testing/logs';
    if ~isfolder(logDir)
        mkdir(logDir);
    end
    timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
    logFile = fullfile(logDir, ['log_' timestamp '.txt']);
    fid = fopen(logFile, 'w');

    try
        fprintf(fid, '[%s] makeWaves.m has changed. Running Earth_LakeSuperior...\n', timestamp);

        % Git metadata (optional: suppress warnings if git isn't available)
        [~, hash] = system('git rev-parse HEAD');
        [~, author] = system('git log -1 --pretty=format:"%an <%ae>"');
        [~, message] = system('git log -1 --pretty=format:%s');

        fprintf(fid, 'Commit Hash: %s\n', strtrim(hash));
        fprintf(fid, 'Author     : %s\n', strtrim(author));
        fprintf(fid, 'Message    : %s\n\n', strtrim(message));

        run('applications/WaveModelling/Earth_LakeSuperior.m');

        if ~exist('output', 'var')
            error('The script Earth_LakeSuperior.m must assign variable "output".');
        end

        baselineFile = 'planetwaves/unit_testing/baselines/Earth_LakeSuperior_output.mat';
        saveNeeded = false;

        if exist(baselineFile, 'file')
            baseline = load(baselineFile);
            outputOld = baseline.output;

            tolerance = 1e-6;
            [isEqual, maxDiff] = compare_outputs(outputOld, output, tolerance);

            if isEqual
                msg = sprintf('Outputs match within tolerance (max diff: %.2e <= %.2e)', maxDiff, tolerance);
                saveNeeded = true;
            else
                msg = sprintf('Output has changed. Max absolute difference: %.2e exceeds tolerance %.2e', maxDiff, tolerance);
            end
        else
            msg = 'No baseline found. Saving current results.';
            saveNeeded = true;
        end

        fprintf('%s\n', msg);
        fprintf(fid, '%s\n', msg);

        if saveNeeded
            if ~isfolder('planetwaves/unit_testing/baselines')
                mkdir('planetwaves/unit_testing/baselines');
            end
            save(baselineFile, 'output');
            fprintf(fid, 'Output saved as new baseline.\n');
        end
    catch ME
        fprintf(fid, 'ERROR: %s\n', ME.message);
        rethrow(ME);
    end

    fclose(fid);
end

function [equal, maxDiff] = compare_outputs(a, b, tol)
    maxDiff = 0;
    equal = true;

    if ~isequal(class(a), class(b))
        equal = false;
        maxDiff = inf;
        return;
    end

    if iscell(a)
        if ~isequal(size(a), size(b))
            equal = false;
            maxDiff = inf;
            return;
        end
        for i = 1:numel(a)
            [eq, diff] = compare_outputs(a{i}, b{i}, tol);
            if ~eq
                equal = false;
            end
            maxDiff = max(maxDiff, diff);
        end
    elseif isnumeric(a)
        if ~isequal(size(a), size(b))
            equal = false;
            maxDiff = inf;
            return;
        end
        diff = abs(a - b);
        maxDiff = max(diff(:));
        if maxDiff > tol
            equal = false;
        end
    else
        warning('Unsupported type in output: %s', class(a));
        equal = false;
        maxDiff = inf;
    end
end
