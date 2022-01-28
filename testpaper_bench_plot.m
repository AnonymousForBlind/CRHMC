curFolder = fullfile(fileparts(mfilename('fullpath')), 'Benchmark');

%plotter = str2func('semilogx');
plotter = str2func('loglog');
%plotter = str2func('plot');

%%%%%%%%%%%%%%%%%%%%%
% %% 1. Mixing Time vs Dim
%subplot(1,2,1)
quantVSdim(curFolder, 'box', plotter, 'go', "mixing")
quantVSdim(curFolder, 'simplex', plotter, 'b+', "mixing")
quantVSdim(curFolder, 'birkhoff', plotter, 'r*', "mixing")
% quantVSdim(curFolder, 'Volesti-test/volesti_test', plotter, polytopes, 'go', 'originalSize', "mixing")
legend('Cube', 'Simplex', 'Birkhoff', 'Location', 'northwest', 'FontSize', 14)
 xlim([10 1e6]); %ylim([10 1e9]);
% 
% %% 2. Sample Time vs Dim
%subplot(1,2,2)
figure;
quantVSdim(curFolder, 'box', plotter, 'go', "sample")
quantVSdim(curFolder, 'simplex', plotter, 'b+', "sample")
quantVSdim(curFolder, 'birkhoff', plotter, 'r*', "sample")
% quantVSdim(curFolder, 'Volesti-test/volesti_test', plotter, polytopes, 'go', 'originalSize', "sample")
legend('Cube', 'Simplex', 'Birkhoff', 'Location', 'northwest', 'FontSize', 14)
 xlim([10 1e6]); %ylim([10 1e9]);
%
% %% 3. Mixing Time vs NNZ
% figure;
% quantVSnnz(curFolder, 'RHMC-test/rhmc_test', plotter, 'ro', "mixing")
% legend('CRHMC', 'Location', 'northwest')
% 
% %% 4. Sample Time vs NNZ
% figure;
% quantVSnnz(curFolder, 'RHMC-test/rhmc_test', plotter, 'ro', "sample")
% legend('CRHMC', 'Location', 'northwest')

%% 5. Preparation Time (RHMC: presolve & CHRR: Fullify and MVE)
% figure;
% roundResult1(curFolder, 'RHMC-test/rhmc_test', plotter, polytopes, 'ro', 'originalSize')
% roundResult2(curFolder, 'CHRR-test/chrr_test', plotter, polytopes, 'bo', 'originalSize')
% legend('CRHMC', 'CHRR', 'Location', 'northwest')
% hold off;
% xlim([50 2*1e5]);














%%%%%%%%%%%%%%%%%%%%%

%% Quantities vs Dimension 
% figure;
% plotResult(curFolder, 'RHMC-test/rhmc_test', plotter, polytopes, 'r--o', 'processedSize')
% plotResult(curFolder, 'CHRR-test/chrr_test', plotter, polytopes, 'b--o', 'processedSize')
% plotResult(curFolder, 'Volesti-test/volesti_test', plotter, polytopes, 'g--o', 'processedSize')
% subplot(1,2,1)
% legend('RHMC', 'CHRR (Ben)', 'CDHR (Volesti)', 'Location', 'northwest')
% subplot(1,2,2)
% legend('RHMC', 'CHRR (Ben)', 'CDHR (Volesti)', 'Location', 'northwest')
% hold off;
% 
% figure;
% plotResult(curFolder, 'RHMC-test/rhmc_test', plotter, polytopes, 'r--o', 'fullSize')
% plotResult(curFolder, 'CHRR-test/chrr_test', plotter, polytopes, 'b--o', 'fullSize')
% plotResult(curFolder, 'Volesti-test/volesti_test', plotter, polytopes, 'g--o', 'fullSize')
% subplot(1,2,1)
% legend('RHMC', 'CHRR (Ben)', 'CDHR (Volesti)', 'Location', 'northwest')
% subplot(1,2,2)
% legend('RHMC', 'CHRR (Ben)', 'CDHR (Volesti)', 'Location', 'northwest')
% hold off;
%%%%%%%

%% Quantities vs nnz(A)
% figure;
% plotResult2(curFolder, 'RHMC-test/rhmc_test', plotter, 'ro')
% %plotResult2(curFolder, 'CHRR-test/chrr_test', plotter, 'bo')
% %plotResult2(curFolder, 'Volesti-test/volesti_test', plotter, 'go')
% subplot(1,2,1)
% %legend('RHMC', 'CHRR (Ben)', 'CDHR (Volesti)', 'Location', 'northwest')
% subplot(1,2,2)
% %legend('RHMC', 'CHRR (Ben)', 'CDHR (Volesti)', 'Location', 'northwest')

%%%%%%%
%% Prepare Time (Presolve for RHMC & MVE for CHRR and CDHR)
% figure;
% roundResult1(curFolder, 'RHMC-test/rhmc_test', plotter, polytopes, 'ro', 'fullSize')
% roundResult2(curFolder, 'CHRR-test/chrr_test', plotter, polytopes, 'bo', 'fullSize')
% legend('RHMC', 'CHRR', 'Location', 'northwest')
% hold off;

function quantVSnnz(curFolder, pathdir, plotter, color, quant)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));
    
    datapath = fullfile(curFolder, '/Instances/1chrr/', '*.mat');
    matfiles2 = dir(datapath);

    numPoly = length(matfiles);
    vnnz = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    step = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        if idx < 22
            load(fullfile(curFolder, 'Instances/1chrr/', matfiles2(idx).name));

            if result.exps.ess >= 10
                vnnz(idx) = nnz(P.A_eq);
                time(idx) = result.exps.sampleTime/result.exps.ess;
                step(idx) = result.exps.step/result.exps.ess;
            else
                fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
            end
        else
            if result.exps.ess >= 10
                vnnz(idx) = 295946;
                time(idx) = result.exps.sampleTime/result.exps.ess;
                step(idx) = result.exps.step/result.exps.ess;
            else
                fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
            end
        end
    end

    [vnnz, seq] = sort(vnnz); time = time(seq); step = step(seq);
    filter = vnnz>0;
    vnnz= vnnz(filter); time = time(filter); step = step(filter);
    
    if quant == "sample"
        h = plotter(vnnz, time, color);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title('Sampling Time', 'FontSize', 15);
        xlabel('NNZ', 'FontSize', 15); ylabel('Time/Sample (s)', 'FontSize', 15);
        hold on;
        fit = polyfit(log(vnnz), log(time), 1);
        vnnz = [100, 1000, 10000, 100000, 1e6];
        z = polyval(fit, log(vnnz));
        plotter(vnnz, exp(z), color(1));
        
        algo = pathdir(1:4);
        fprintf(strcat(algo, ": Time/NNZ = %f\n"), fit(1));
    elseif quant == "mixing"
        h = plotter(vnnz, step, color);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title('Mixing Time', 'FontSize', 15);
        xlabel('NNZ', 'FontSize', 15); ylabel('Step/Sample', 'FontSize', 15);
        hold on;
        fit = polyfit(log(vnnz), log(step), 1);
        vnnz = [100, 1000, 10000, 100000, 1e6];
        z = polyval(fit, log(vnnz));
        plotter(vnnz, exp(z), color(1));
        
        algo = pathdir(1:4);
        fprintf(strcat(algo, ": Step/NNZ= %f\n"), fit(1));
    end
end

function quantVSdim(curFolder, pathdir, plotter, color, quant)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));

    numPoly = length(matfiles);
    dim = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    step = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));

        if result.exps.ess >= 1
            dim(idx) = result.exps.dim;
            time(idx) = result.exps.sampleTime/result.exps.ess;
            step(idx) = result.exps.step/result.exps.ess;
        else
            fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
        end
    end

    [dim, seq] = sort(dim); time = time(seq); step = step(seq);
    filter = dim>0;
    dim = dim(filter); time = time(filter); step = step(filter);
    
    if quant == "sample"
        h = plotter(dim, time, color);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title('Sampling Time', 'FontSize', 15);
        xlabel('Dimension', 'FontSize', 15); ylabel('Time/Sample (s)', 'FontSize', 15);
        hold on;
        fit = polyfit(log(dim), log(time), 1);
        newdim = dim; newdim(end+1)=1000000; newdim(1)=10;
        %newdim = [50, 10, 100, 1000, 10000, 100000, 200000];
        z = polyval(fit, log(newdim));
        plotter(newdim, exp(z), color(1));
        fprintf(strcat(pathdir, ": Time/Dim = %f & Coeff = %f\n"), fit(1), exp(fit(2)));
    elseif quant == "mixing"
        h = plotter(dim, step, color);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        title('Mixing Rate', 'FontSize', 15);
        xlabel('Dimension', 'FontSize', 15); ylabel('Step/Sample', 'FontSize', 15);
        hold on;
        fit = polyfit(log(dim), log(step), 1);
        newdim = dim; newdim(end+1)=1000000; newdim(1)=10;
        %newdim = [50, 10, 100, 1000, 10000, 100000, 200000];
        z = polyval(fit, log(newdim));
        plotter(newdim, exp(z), color(1));
        fprintf(strcat(pathdir, ": Step/Dim = %f & Coeff = %f\n"), fit(1), exp(fit(2)));
    end
end

function roundResult2(curFolder, pathdir, plotter, polytopes, color, dimOpt)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));
    datapath = fullfile(curFolder, '/Instances/1chrr/', '*.mat');
    matfiles2 = dir(datapath);
    
    numPoly = length(matfiles2);
    dim = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    
    for idx = 1:numPoly
        load(fullfile(curFolder, 'Instances/1chrr/', matfiles2(idx).name));
        name = matfiles(idx).name(17:end-4); 
        dim(idx) = polytopes.(name).(dimOpt)(2);
        time(idx) = P.roundTime;
    end

    [dim, seq] = sort(dim); time = time(seq);
    filter = dim>0;
    dim = dim(filter); time = time(filter);
    h = plotter(dim, time, color);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    title('Preparation Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Time', 'FontSize', 15);
    hold on;
    fit = polyfit(log(dim), log(time), 1);
    newdim = [60, 100, 1000, 10000, 100000, 1e6];
    z = polyval(fit, log(newdim));
    plotter(newdim, exp(z), color(1));
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": RoundTime/Dim = %f & Coeff = %f\n"), fit(1), exp(fit(2)));
end

function roundResult1(curFolder, pathdir, plotter, polytopes, color, dimOpt)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));

    numPoly = length(matfiles);
    dim = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        name = matfiles(idx).name(17:end-4); 

        if result.exps.ess >= 10
            dim(idx) = polytopes.(name).(dimOpt)(2);
            time(idx) = result.exps.roundTime;
        else
            fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
        end
    end

    [dim, seq] = sort(dim); time = time(seq);
    filter = dim>0;
    dim = dim(filter); time = time(filter);
    h = plotter(dim, time, color);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    title('Preparation Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Time', 'FontSize', 15);
    hold on;
    fit = polyfit(log(dim), log(time), 1);
    newdim = [60, 100, 1000, 10000, 100000, 1e6];
    z = polyval(fit, log(newdim));
    plotter(newdim, exp(z), color(1));
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": RoundTime/Dim = %f & Coeff = %f\n"), fit(1), exp(fit(2)));
end

function plotResult2(curFolder, pathdir, plotter, color)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));
    
    datapath = fullfile(curFolder, '/Instances/1chrr/', '*.mat');
    matfiles2 = dir(datapath);

    numPoly = length(matfiles);
    vnnz = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    step = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        if idx < 22
            load(fullfile(curFolder, 'Instances/1chrr/', matfiles2(idx).name));

            if result.exps.ess >= 10
                vnnz(idx) = nnz(P.A_eq);
                time(idx) = result.exps.sampleTime/result.exps.ess;
                step(idx) = result.exps.step/result.exps.ess;
            else
                fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
            end
        else
            if result.exps.ess >= 10
                vnnz(idx) = 295946;
                time(idx) = result.exps.sampleTime/result.exps.ess;
                step(idx) = result.exps.step/result.exps.ess;
            else
                fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
            end
        end
    end

    [vnnz, seq] = sort(vnnz); time = time(seq); step = step(seq);
    filter = vnnz>0;
    vnnz= vnnz(filter); time = time(filter); step = step(filter);
    subplot(1,2,1);
    plotter(vnnz, time, color)
    title('Sampling Time', 'FontSize', 15);
    xlabel('NNZ', 'FontSize', 15); ylabel('Time/Sample (s)', 'FontSize', 15);
    hold on;

    subplot(1,2,2)
    plotter(vnnz, step, color)
    title('Mixing Time', 'FontSize', 15);
    xlabel('NNZ', 'FontSize', 15); ylabel('Step/Sample', 'FontSize', 15);
    hold on;

    fit1 = polyfit(log(vnnz), log(time), 1);
    fit2 = polyfit(log(vnnz), log(step), 1);
    
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": Time/NNZ = %f, Step/NNZ= %f\n"), fit1(1), fit2(1));
end

function plotResult_RHMC(curFolder, pathdir, plotter, polytopes, color)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));

    numPoly = length(matfiles);
    dim = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    step = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        name = matfiles(idx).name(17:end-4); 

        if result.exps.ess >= 10
            dim(idx) = polytopes.(name).processedSize(2);
            time(idx) = result.exps.sampleTime/result.exps.ess;
            step(idx) = result.exps.step/result.exps.ess;
        else
            fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
        end
    end

    [dim, seq] = sort(dim); time = time(seq); step = step(seq);
    subplot(1,2,1);
    plotter(dim, time, color)
    title('Sampling Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Time/Sample (s)', 'FontSize', 15);
    hold on;

    subplot(1,2,2)
    plotter(dim, step, color)
    title('Mixing Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Step/Sample', 'FontSize', 15);
    hold on;

    fit1 = polyfit(log(dim), log(time), 1);
    fit2 = polyfit(log(dim), log(step), 1);
    
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": Time/Dim = %f, Step/Dim = %f\n"), fit1(1), fit2(1));
end

function plotResult(curFolder, pathdir, plotter, polytopes, color, dimOpt)
    matfiles = dir(fullfile(fullfile(curFolder, pathdir), '*.mat'));

    numPoly = length(matfiles);
    dim = zeros(numPoly, 1);
    time = zeros(numPoly, 1);
    step = zeros(numPoly, 1);
    for idx = 1:numPoly
        result = load(strcat(matfiles(1).folder,'/',matfiles(idx).name));
        name = matfiles(idx).name(17:end-4); 

        if result.exps.ess >= 10
            dim(idx) = polytopes.(name).(dimOpt)(2);
            time(idx) = result.exps.sampleTime/result.exps.ess;
            step(idx) = result.exps.step/result.exps.ess;
        else
            fprintf("%s: Ess %d is too small\n", matfiles(idx).name, result.exps.ess) 
        end
    end

    [dim, seq] = sort(dim); time = time(seq); step = step(seq);
    filter = dim>0;
    dim = dim(filter); time = time(filter); step = step(filter);
    subplot(1,2,1);
    plotter(dim, time, color)
    title('Sampling Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Time/Sample (s)', 'FontSize', 15);
    hold on;

    subplot(1,2,2)
    plotter(dim, step, color)
    title('Mixing Time', 'FontSize', 15);
    xlabel('Dimension', 'FontSize', 15); ylabel('Step/Sample', 'FontSize', 15);
    hold on;

    fit1 = polyfit(log(dim), log(time), 1);
    fit2 = polyfit(log(dim), log(step), 1);
    
    algo = pathdir(1:4);
    fprintf(strcat(algo, ": Time/Dim = %f, Step/Dim = %f\n"), fit1(1), fit2(1));
end

