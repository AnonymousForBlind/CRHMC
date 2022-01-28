maxNumCompThreads(1);
curFolder = fileparts(mfilename('fullpath'));
matfiles = dir(fullfile(fullfile(curFolder, '/SavedPoints/'), '*.mat'));

%%
idx = {[1:21]};
for c = idx{1}
     try
        fprintf('%d: %s\n',c, matfiles(c).name);
        model = load(fullfile(curFolder, '/SavedPoints', matfiles(c).name));
        model = cellfun(@(x)(model.(x)),fieldnames(model));
        name = matfiles(c).name(8:end); name(2:3) = 'hr';
        load(fullfile(fileparts(curFolder), 'Instances/1chrr', name));

        tic;
        model.points = P.N*model.points + P.p_shift;
        samples = P.T*model.points + P.y;
        transTime = toc;
        
        ess = effectiveSampleSize(samples);
        nSamples = min(ess);
        gap = round(size(samples,2) / nSamples);
        samples = samples(:,gap:gap:end);
        %assert(norm(P.originalb - P.originalA * samples, Inf) <= 1e-8)

        exps = struct; exps.dim = model.dim; 
        exps.ess = min(effective_sample_size(samples));
        exps.sampleTime = model.time + transTime;
        exps.step = model.step;

        save(fullfile(curFolder, strcat('/volesti_test/result_', matfiles(c).name(8:end))), 'exps')
    catch e
        fprintf('Error occurs. Skip model %d\n', c);
        fprintf(1,'The identifier was:\n%s\n', e.identifier);
        fprintf(1,'There was an error! The message was:\n%s\n', e.message);
        continue;
    end
end



