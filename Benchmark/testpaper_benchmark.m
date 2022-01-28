%% This file tests the rhmc package on benchmark instances
% It generates files "result_(modelname)_(dim)" with a struct with the following 
% properties in the folder "rhmc_test".
% 1. dim
% 2. ess
% 3. preparation time
% 4. (total) sample time
% 5. (total) # steps

curFolder = fileparts(mfilename('fullpath'));
instances = {'box', 'simplex', 'birkhoff'};
dimBoxSimplex = [10:10:90, 100:100:900, 1000:1000:10000, 50000:50000:100000, 500000:500000:1000000];
dimBirkhoff = [10:10:90, 100:100:900, 1000:1000:10000, 50000:50000:100000, 500000:500000:1000000];
numSamples = 500;

%numSamples = 10;

for c = [1, 2, 3]
    inst = instances{c};
    if strcmp('box', inst) || strcmp('simplex', inst)
        dimInput = dimBoxSimplex;
    elseif strcmp('birkhoff', inst)
        dimInput = dimBirkhoff;
    else
        error('Not supported')
    end

    for dim = dimInput
        fprintf('%s: %d\n', inst, dim);
        disp(datetime('now'));
        P = loadProblem(['basic\' inst '@' int2str(dim)]);

        opts = default_options();
        opts.module = {'MixingTimeEstimator', 'MemoryStorage', 'DynamicRegularizer', 'DynamicStepSize', 'DynamicWeight', 'ProgressBar', 'DebugLogger'};
        opts.MemoryStorage.memoryLimit = 8*1024*1024*1024;
        opts.simdLen = 1;
        opts.ProgressBar.refreshInterval = 600;
        opts.maxTime = 3600*24;
        opts.seed = 1;
        opts.logging = 'testpaper_benchmark.log'; % Output the debug log to testpaper_benchmark.log
        
        o = sample(P, numSamples, opts);

        exps = struct; 
        exps.dim = size(P.Aeq, 2); 
        exps.ess = size(o.samples, 2);
        exps.numelA = numel(P.Aeq);
        exps.nnzA = nnz(P.Aeq);
        exps.roundTime = o.prepareTime;
        exps.sampleTime = o.sampleTime;
        exps.step = o.totalStep * o.opts.simdLen;
        fprintf('ESS: %f\n', exps.ess);

        save(fullfile(curFolder, strcat(inst, '/result_rhmc_', inst, int2str(dim))), 'exps');
    end
end