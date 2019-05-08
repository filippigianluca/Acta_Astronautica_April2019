function callMex(funName, libs, outputDir, slash, obj)

%%
mexName = cell(1,length(funName));
for i = 1:length(funName)
    mexName{i} = ['mex' upper(funName{i}(1)) funName{i}(2:end)];
end

%%
fprintf('\n');
for i = 1:length(funName)
    fprintf([funName{i} ': ']);
    str = ['mex src' slash mexName{i} '.c '];
    for j = 1:length(libs)
        str = [str libs{j} obj ' '];
    end
    str = [str '-outdir ' outputDir ' -output ' funName{i} ' -DMEXCOMPILE'];
    eval(str)
    fprintf('done\n');
end