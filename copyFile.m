clear;
clc;
for i = 21:35
    sourceFile1 = fullfile('Biomodels', sprintf('Model%d',i), 'f.m');
    sourceFile2 = fullfile('Biomodels', sprintf('Model%d',i), 'Initialize_ODE.m');
    fprintf('i = %d   ',i);
    for j = 1:500
        dest = fullfile('Config', sprintf('Model%d',i), sprintf('Copy%d',j));
        copyfile(sourceFile1, dest);
        copyfile(sourceFile2, dest);
        fprintf('j = %d\n',j);
    end
end