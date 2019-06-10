function add_path()

dirc = [matlabroot filesep 'toolbox' filesep 'maxwellfdfd-master'];

% dirc = [pwd filesep 'maxwellfdfd-master'];

if ~contains(path, dirc)
    addpath(dirc);
    addpath([dirc, filesep, 'base']);
    addpath([dirc, filesep, 'diff']);
    addpath([dirc, filesep, 'dielconst']);
    addpath([dirc, filesep, 'grid']);
    addpath([dirc, filesep, 'integ']);
    addpath([dirc, filesep, 'io']);
    addpath([dirc, filesep, 'material']);
    addpath([dirc, filesep, 'modesolver']);
    addpath([dirc, filesep, 'petsc']);
    addpath([dirc, filesep, 'shape']);
    addpath([dirc, filesep, 'source']);
    addpath([dirc, filesep, 'vis']);
end

dirc = [pwd, filesep, 'classes'];
if ~contains(path, dirc); addpath(dirc); end

dir = [pwd, filesep, 'functions'];
if ~contains(path, dir); addpath(dir); end

dir = [pwd, filesep, 'FresnelData'];
if ~contains(path, dir); addpath(dir); end

dir = [pwd, filesep, 'FresnelData', filesep, 'The_first_opus'];
if ~contains(path, dir); addpath(dir); end

dir = [pwd, filesep, 'FresnelData', filesep, 'The_second_opus'];
if ~contains(path, dir); addpath(dir); end

dir = [pwd, filesep, 'FresnelData', filesep, 'The_third_opus'];
if ~contains(path, dir); addpath(dir); end

dir = [pwd, filesep, 'ndSparse_G3_2013_03_13'];
if ~contains(path, dir); addpath(dir); end

dir = [pwd, filesep, 'DGradient'];
if ~contains(path, dir); addpath(dir); end


