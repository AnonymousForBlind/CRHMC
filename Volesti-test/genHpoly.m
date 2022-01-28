initSampler

%matfiles = dir(fullfile('~/Dropbox/fork-cobratoolbox/test/models/mat/', '*.mat'));
curFolder = fileparts(mfilename('fullpath'));
matfiles = dir(fullfile(fullfile(curFolder, '/mat'), '*.mat'));
%%
idx = {[2]};
for c = idx{1}
    try
        fprintf('%d: %s\n',c, matfiles(c).name);
        %model = load(fullfile('~/Dropbox/fork-cobratoolbox/test/models/mat/', matfiles(c).name));
        model = load(fullfile(curFolder, '/mat', matfiles(c).name));
        model = cellfun(@(x)(model.(x)),fieldnames(model));
        %[modelSampling,samples] = sampleCbModel(model, [], 'RHMC', opt);
        P = struct;        
        if (~isfield(model,'S') || ~isfield(model,'b'))
            error('You need to define both model.S and model.b');
        else
            P.Aeq = model.S;
            P.beq = model.b;
        end
        if isfield(model,'lb')
            P.lb = model.lb;
        end
        if isfield(model,'ub')
            P.ub = model.ub;
        end
        if isfield(model,'dsense')
            I = (model.dsense == 'E');
            P.Aeq = [P.Aeq; model.C(I,:)];
            P.beq = [P.beq; model.d(I)];
            P.Aineq = model.C(~I,:);
            P.bineq = model.d(~I,:);
            flip = 1-2*(model.dsense(~I) == 'G');
            P.Aineq = flip.*P.Aineq;
            P.bineq = flip.*P.bineq;
        end

        o = sample(P, 1);

        %% Convert to H polytope
        % Ax=b, lb <= x <= ub
        A = o.problem.A;
        b = o.problem.b;
        lb = o.problem.barrier.lb';
        ub = o.problem.barrier.ub';
        x = o.problem.center;

        % 1) shift the center to 0.
        % the problem becomes Ax = 0, lb <= x <= ub
        lb = lb - x;
        ub = ub - x;

        % 2) rescale the polytope
        s = min(-lb, ub);
        A = A .* s';
        lb = lb ./ s;
        ub = ub ./ s;

        % 3) convert to {Hx<=b}
        % reason: {Ax = 0, lb <= x <= ub} = {lb <= Cx <= ub}
        C = null(full(A));
        H = [C; -C];
        b1 = [ub; -lb];

        %% Save Polytope
        poly=cell(2,1);
        poly{1} = H;
        poly{2} = b1;
        %save(strcat('HPolytope_', matfiles(c).name),'poly')
        save(strcat(fullfile(curFolder, '/Polytopes/'), 'HPolytope_', matfiles(c).name), 'poly');
    catch
        fprintf('Error occurs. Skip model %d\n', c);
        continue;
    end
end



