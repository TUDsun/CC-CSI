function [A, b, hfcn_GfromF, hfcn_Op, D] = mycreate_eqTE(eqtype, pml, omega, eps_cell, mu_cell, s_factor_cell, J_cell, M_cell, grid3d)

chkarg(istypesizeof(eqtype, 'EquationType'), '"eqtype" should be instance of EquationType.');
chkarg(istypesizeof(pml, 'PML'), '"pml" should be instance of PML.');
chkarg(istypesizeof(omega, 'complex'), '"omega" should be complex.');
chkarg(istypesizeof(grid3d, 'Grid3d'), '"grid3d" should be instance of Grid3d.');

N = grid3d.N;
src_n = size(J_cell,1);
if src_n ==0
    src_n = 1;
end
chkarg(istypesizeof(eps_cell, 'complexcell', [1 Axis.count], N), ...
    '"eps_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements', ...
    Axis.count, N(Axis.x), N(Axis.y), N(Axis.z));
chkarg(istypesizeof(mu_cell, 'complexcell', [1 Axis.count], N), ...
    '"mu_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements', ...
    Axis.count, N(Axis.x), N(Axis.y), N(Axis.z));
chkarg(istypesizeof(s_factor_cell, 'complexcell', [Axis.count GT.count], [1 0]), ...
    '"mu_cell" should be %d-by-%d cell array whose each element is row vector with complex elements', ...
    Axis.count, GT.count);
chkarg(istypesizeof(J_cell, 'ndSparsecell', [src_n Axis.count], N), ...
    '"J_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements', ...
    Axis.count, N(Axis.x), N(Axis.y), N(Axis.z));
chkarg(istypesizeof(M_cell, 'ndSparsecell', [src_n Axis.count], N), ...
    '"M_cell" should be length-%d row cell array whose each element is %d-by-%d-by-%d array with complex elements', ...
    Axis.count, N(Axis.x), N(Axis.y), N(Axis.z));

% Reorder the indices of the elements of matrices and vectors to reduce the bandwidth of A.
r = reordering_indices(Axis.count, N);
r = reshape(r,3,length(r)/3);
r(3,:) = [];
r = r(:);
% Construct curls
dl_factor_cell = [];
if pml == PML.sc
    dl_factor_cell = s_factor_cell;
end

ge = eqtype.ge;
[Ce, Cm] = create_curls(ge, dl_factor_cell, grid3d);

% Construct material parameters
if pml == PML.u
    [smx, smy, smz] = ndgrid(s_factor_cell{Axis.x, alter(ge)}, s_factor_cell{Axis.y, alter(ge)}, s_factor_cell{Axis.z, alter(ge)});
    [sex, sey, sez] = ndgrid(s_factor_cell{Axis.x, ge}, s_factor_cell{Axis.y, ge}, s_factor_cell{Axis.z, ge});
    sm = {smx, smy, smz};
    se = {sex, sey, sez};
    mu_cell = mult_vec(mu_cell, sm([Axis.y Axis.z Axis.x]));
    mu_cell = mult_vec(mu_cell, sm([Axis.z Axis.x Axis.y]));
    mu_cell = div_vec(mu_cell, se([Axis.x Axis.y Axis.z]));
    
    eps_cell = mult_vec(eps_cell, se([Axis.y Axis.z Axis.x]));
    eps_cell = mult_vec(eps_cell, se([Axis.z Axis.x Axis.y]));
    eps_cell = div_vec(eps_cell, sm([Axis.x Axis.y Axis.z]));
end

mu = [mu_cell{Axis.x}(:) ; mu_cell{Axis.y}(:) ; mu_cell{Axis.z}(:)];
eps = [eps_cell{Axis.x}(:) ; eps_cell{Axis.y}(:) ; eps_cell{Axis.z}(:)];
% j = cell(1,src_n); m = cell(1,src_n);
j = sparse([]);
m = sparse([]);
for ii = 1:src_n
    tmp1 = [J_cell{ii,Axis.x}(:) ; J_cell{ii,Axis.y}(:) ; J_cell{ii,Axis.z}(:)];
    tmp2 = [M_cell{ii,Axis.x}(:) ; M_cell{ii,Axis.y}(:) ; M_cell{ii,Axis.z}(:)];
    j = [j, sparse(tmp1)];
    m = [m, sparse(tmp2)];
end
% j(:,1) = [];
% m(:,1) = [];

% Mask elements corresponding to PEC.
ind_pec = isinf(abs(eps));
eps(ind_pec) = 1;
pm = ones(size(ind_pec));  % PEC mask
pm(ind_pec) = 0;

% Ce = Ce(r,r);
% Cm = Cm(r,r);

% mu = mu(r);
% eps = eps(r);

% pm = pm(r);

PM = create_spdiag(pm);

test_vector_identity = false;
if test_vector_identity
    [Dive, Divm] = create_divs(ge, dl_factor_cell, grid3d);
    Dive = Dive(:,r);
    Divm = Divm(:,r);
    fprintf('norm(Dive * Cm, 1) = %e\n', norm(Dive * Cm, 1));
    fprintf('norm(Divm * Ce, 1) = %e\n', norm(Divm * Ce, 1));
end
if eqtype.f == FT.e
    INV_MU = create_spdiag(1./mu);  % when mu has Inf, "MU \ Mat" complains about singularity
    D = create_spdiag(eps);
    
    A = PM * (Cm * INV_MU * Ce) * PM - omega^2 * D;
    hfcn_A = @(e) pm .* (Cm * ((Ce * (pm .* e)) ./ mu)) - omega^2 * (eps .* e);
    hfcn_Atr = @(e) pm .* (Ce.' * ((Cm.' * (pm .* e)) ./ mu)) - omega^2 * (eps .* e);
    
    b = -1i*omega*j - Cm*(m./repmat(mu,1,src_n));
    
    hfcn_GfromF = @(e) (Ce * e + m) ./ (-1i*omega*repmat(mu,1,src_n));
else  % eqtype.f == FT.h
    INV_EPS = create_spdiag(1./eps);  % when mu has Inf, "MU \ Mat" complains about singularity
    D = create_spdiag(mu);
    
    A = (Ce * INV_EPS * Cm) - omega^2 * D;
    hfcn_A = @(h) Ce * ((Cm * h) ./ eps) - omega^2 * (mu .* h);
    hfcn_Atr = @(h) Cm.' * ((Ce.' * h) ./ eps) - omega^2 * (mu .* h);
    
    b = -1i*omega*m + Ce*(j./repmat(mu,1,src_n));
    
    hfcn_GfromF = @(h) (Cm * h - j) ./ (1i*omega*repmat(mu,1,src_n));
end

    function y = Op(x, transp_flag)
        if strcmp(transp_flag,'transp')
            y = hfcn_Atr(x);
        elseif strcmp(transp_flag,'notransp')
            y = hfcn_A(x);
        end
    end

hfcn_Op = @Op;
ind = size(A,1)/3;
A(end - ind + 1:end,:) = [];
A(:,end - ind + 1:end) = [];
b(end - ind + 1:end,:) = [];
A = A(r, r);
b = b(r, :);
end
