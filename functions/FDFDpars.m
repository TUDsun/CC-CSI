function [osc, grid3d, M_S, A, b] = FDFDpars(varargin)
DEFAULT_METHOD = 'direct';  % 'direct', 'gpu', 'aws', 'inputfile'
% Set solver options.
iarg = nargin; arg = varargin{iarg};
% inspect_only = false;
% if istypesizeof(arg, 'logical')
%     inspect_only = arg;
%     iarg = iarg - 1; arg = varargin{iarg};
% end

is_solveropts = false;
if isenum(arg) || ischar(arg)
    polar = arg;
    iarg = iarg - 1;
    arg = varargin{iarg};
end
if istypesizeof(arg, 'struct')
    solveropts = arg;
    is_solveropts = true;
    iarg = iarg - 1;  % arg = varargin{iarg};
end

if ~is_solveropts || ~isfield(solveropts, 'showstruct')
    solveropts.showstruct = true;
end

if ~is_solveropts || ~isfield(solveropts, 'method')
    solveropts.method = DEFAULT_METHOD;
end

if is_solveropts && isequal(solveropts.method, 'inputfile')
    chkarg(isfield(solveropts, 'filenamebase'), '"solveropts" should have "filenamebase" field.');
end

if ~is_solveropts || ~isfield(solveropts, 'maxit')
    % 		solveropts.maxit = intmax;
    solveropts.maxit = 1e8;
else
    chkarg(istypesizeof(solveropts.maxit, 'real') && solveropts.maxit > 0, ...
        'solveropts.maxit should be positive.');
end

if ~is_solveropts || ~isfield(solveropts, 'tol')
    solveropts.tol = 1e-6;
else
    chkarg(istypesizeof(solveropts.tol, 'real') && solveropts.tol > 0, ...
        'solveropts.tol should be positive.');
end

if ~is_solveropts || ~isfield(solveropts, 'eqtype')
    solveropts.eqtype = EquationType(FT.e, GT.prim);
else
    chkarg(istypesizeof(solveropts.eqtype, 'EquationType'), ...
        'solveropts.eqtype should be instance of EquationType.');
end

if ~is_solveropts || ~isfield(solveropts, 'pml')
    solveropts.pml = PML.sc;
else
    chkarg(istypesizeof(solveropts.pml, 'PML'), ...
        'solveropts.pml should be instance of PML.');
end

if ~is_solveropts || ~isfield(solveropts, 'returnAandb')
    solveropts.returnAandb = false;
else
    chkarg(istypesizeof(solveropts.returnAandb, 'logical'), ...
        'solveropts.returnAandb should be logical.');
end

if ~is_solveropts || ~isfield(solveropts, 'returnDiv')
    solveropts.returnDiv = false;
else
    chkarg(istypesizeof(solveropts.returnDiv, 'logical'), ...
        'solveropts.returnDiv should be logical.');
end

chkarg(iarg > 0, 'first argument is not correct.');

fprintf('E-field grid type: %s\n', char(solveropts.eqtype.ge));
pm = ProgMark();

% Build the system.
% Make sure to pass the first consecutive elements of varargin to
% build_system() for correct error reports.
[osc, grid3d, s_factor, Eps, mu, J, M, M_S, ~, ~, ~, ~, ~] = ...
    csi_build_system(solveropts.eqtype.ge, solveropts.pml, varargin{1:iarg}, pm, solveropts.returnAandb);
% [osc, grid3d, s_factor, eps, mu, J, M, obj_array, src_array, mat_array, eps_node, mu_node] = ...
%     mybuild_system(solveropts.eqtype.ge, solveropts.pml, varargin{3:iarg}, pm);
% [E, H] = solve_eq_direct(solveropts.eqtype, solveropts.pml, omega, eps, mu, s_factor, J, M, grid3d);

omega = osc.in_omega0();
if solveropts.returnAandb
    switch polar
        
        case PT.TM
            [A, b, ~, ~, ~] = mycreate_eqTM(solveropts.eqtype, solveropts.pml, omega, Eps, mu, s_factor, J, M, grid3d);
            ry = 1 : size(M_S, 2);
            ry = reshape(ry, 3, length(ry) / 3);
            ry(1 : 2, :) = [];
            ry = ry(:);
            rx = 1 : size(M_S, 1);
            rx = reshape(rx, 3, length(rx) / 3);
            rx(1 : 2, :) = [];
            rx = rx(:);
            M_S = M_S(rx, ry);
        case PT.TE
            [A, b, ~, ~, ~] = mycreate_eqTE(solveropts.eqtype, solveropts.pml, omega, Eps, mu, s_factor, J, M, grid3d);
            ry = 1 : size(M_S, 2);
            ry = reshape(ry, 3, length(ry) / 3);
            ry(3, :) = [];
            ry = ry(:);
            rx = 1 : size(M_S, 1);
            rx = reshape(rx, 3, length(rx) / 3);
            rx(3, :) = [];
            rx = rx(:);
            M_S = M_S(rx, ry);
        otherwise
            %             case PT.FULL
            [A, b, ~, ~, ~] = mycreate_eq(solveropts.eqtype, solveropts.pml, omega, Eps, mu, s_factor, J, M, grid3d);
            %                 exception = MException('Wrong polarization');
            %                 throwAsCaller(exception);
    end
    
else
    A = [];
    b = [];
end
