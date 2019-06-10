%% CSI
% function of the constrast source inversion method




function [Mchi, time, chi, vJ, eTot] = CSI(Emea, Phi, A, chi, vJ, eIncInv, eTot, pars, epsinv)

interval    = 64;
Mchi        = cell(ceil(pars.itenum / interval) + 1, 1);
kk          = 1;

if isempty(epsinv); epsinv  = chi{1}; end

NP          = 1 : 3;
Nfre        = length(eIncInv);
N           = pars.Ninv;
ri          = pars.rd;
rchi        = pars.rchi;
ptype       = pars.ptype;
% epsinv      = forward_clb(epsinv, N, 1, ptype);
% epsinv      = sum(epsinv, 1) / NP(ptype);
% epsinv      = epsinv.';
bgchi       = pars.bgchi;
omega       = pars.omega;
r           = pars.rmea;
disflag     = pars.disflag;
terflag     = pars.terflag;

mchi        = cellfun(@(x) reshape(x, N), chi, ...
    'UniformOutput', false);

Dchi        = cellfun(@(x) back_clb(x, N, ptype), mchi, ...
    'UniformOutput', false);

Nsrc        = size(Emea{1}, 2);

vnu         = cell(1, Nfre);
vg          = cell(1, Nfre);



for ii = 1 : Nfre
    %     eTot{ii}(ri, :) 	= 0;
    vnu{ii}     = zeros(size(eIncInv{1}));
    vg{ii}      = zeros(size(eIncInv{1}));
    
end

etaS    = cellfun(@(x) norm(x, 'fro'), Emea, ...
    'UniformOutput',  false);

etaS    = cellfun(@(x) 1 / (x ^ 2), etaS, ...
    'UniformOutput', false);

chieTot = cellfun(@(x, y) x * y, Dchi, eTot, ...
    'UniformOutput', false);

chieInc = cellfun(@(x, y) x * y, Dchi, eIncInv, ...
    'UniformOutput', false);

vrhotmp = cellfun(@(x, y, z) x - y * z, Emea, Phi, vJ, ...
    'UniformOutput', false);

vrho    = cellfun(@(x) x - x, Emea, ...
    'UniformOutput', false);

% vxi         = vrho;
for ii = 1 : Nfre; vrho{ii}(r{ii}) = vrhotmp{ii}(r{ii}); end

vr      	= cellfun(@(x, y) x - y, chieTot, vJ, ...
    'UniformOutput', false);

for ii= 1 : Nfre; vr{ii}(ri, :) = 0; end

ctmp     	= cellfun(@(x) norm(x, 'fro') ^ 2, chieInc, ...
    'UniformOutput', false);

etaD     	= cellfun(@(x) x / (x * x + eps), ctmp, ...
    'UniformOutput', false);

err         = 1e3; tot = 5;
sb0         = 0; sb = 0;  Rn = 0; jj = 1; FC = 1;
logH        = '%5s%10s%10s%10s%13s%13s%13s%13s\n';
logB        = '%5i%10.2e%10.2e%10.2e%13.4e%13.4e%13.4e%13.4e\n';
disp('CSI-MF iteration starts \n')
fprintf(logH, 'iters', 'FS', 'FD', 'FC', 'err', 'ChiErr', 'sb0', 'sb')

%% Iteration starts
if pars.recflag
    diary([pars.str '.txt'])
    diary on
end
startTime   = tic;
while jj < pars.itenum + 1
    % +++++++++++++++++ update contrast source +++++++++++++++++
    % vgs
    vgs     = cellfun(@(y, z, a) -a * y' * z, Phi, vrho, etaS, ...
        'UniformOutput', false);
    
    chir    = cellfun(@(x, y, z) x' \ (y' * z), A, Dchi, vr, ...
        'UniformOutput', false); % A' \ (Dchi' * vr)
    
    vgOld   = vg;
    
    vg      = cellfun(@(x, y, z, a) x - y * (z - a), vgs, etaD, vr, chir, ...
        'UniformOutput', false); %vgs - etaD * (vr - chir)
    
    for ii = 1 : Nfre; vg{ii}(ri, :) = 0; end
    % vNu
    ctmp    = cellfun(@(x) norm(x, 'fro') ^ 2, vgOld, ...
        'UniformOutput', false);
    
    yeta    = cellfun(@(x) x / (x + eps), ctmp, ...
        'UniformOutput', false);
    
    t_v     = cellfun(@(x, y, z, a) x * sum(sum((y - z) .* conj(y))) / (a + eps), yeta, vg, vgOld, ctmp, ...
        'UniformOutput', false);
    
    vnuOld  = vnu;
    
    vnu     = cellfun(@(x, y, a) x + y * real(a), vg, vnuOld, t_v, ...
        'UniformOutput', false);
    
    eNu     = cellfun(@(x, y) x \ y, A, vnu, ...
        'UniformOutput', false);
    %     for ii = 1 : Nfre; eNu{ii}(ri,:) = 0; end
    % va
    
    nuchie  = cellfun(@(x, y, z) x - y * z, vnu, Dchi, eNu, ...
        'UniformOutput', false);
    
    va      = cellfun(@(x, y, z, a, b, c) -real(sum(conj(x) .* y))./(z * sum(abs(a * x) .^ 2) + b * sum(abs(c) .^ 2)), ...
        vnu, vg, etaS, Phi, etaD, nuchie, ...
        'UniformOutput', false);
    
    vJ      = cellfun(@(x, y, z) x + y * diag(z), vJ, vnu, va, ...
        'UniformOutput', false);
    
    eTot    = cellfun(@(x, y, z) x + y * diag(z), eTot, eNu, va, ...
        'UniformOutput', false);
    
    % +++++++++++++ Update Chi with Positivity +++++++++++++++++
    chi{1}          = get_chiMF(vJ, eTot, N, Nsrc, omega, ptype);
    chi{1}          = OrthgPro(chi{1}, bgchi);
    chi{1}(rchi)    = 0;
    
    chi     = cellfun(@(x) real(chi{1}) + 1j * omega{1} / x * imag(chi{1}), omega, ...
        'UniformOutput', false);
    
    mchi    = cellfun(@(x) reshape(x, N), chi, ...
        'UniformOutput', false);
    
    Dchi    = cellfun(@(x) back_clb(x, N, ptype), mchi, ...
        'UniformOutput', false);
    
    chieTot = cellfun(@(x, y) x * y, Dchi, eTot, ...
        'UniformOutput', false);
    
    chieInc = cellfun(@(x, y) x * y, Dchi, eIncInv, ...
        'UniformOutput', false);
    
    vrhotmp = cellfun(@(x, y, z) x - y * z, Emea, Phi, vJ, ...
        'UniformOutput', false);
    
    for ii = 1 : Nfre; vrho{ii}(r{ii}) = vrhotmp{ii}(r{ii}); end
    
    vr      = cellfun(@(x, y) x - y, chieTot, vJ, ...
        'UniformOutput', false);
    
    for ii= 1 : Nfre; vr{ii}(ri, :) = 0; end
    
    ctmp    = cellfun(@(x) norm(x, 'fro') ^ 2, chieInc, ...
        'UniformOutput', false);
    
    etaD    = cellfun(@(x) x / (x * x + eps), ctmp, ...
        'UniformOutput', false);
    
    FSC         = cellfun(@(x, a) a * norm(x, 'fro') ^ 2, vrho, etaS, ...
        'UniformOutput', false);
    
    FdC     = cellfun(@(x) norm(x, 'fro') ^ 2, vr, ...
        'UniformOutput', false);
    
    FDC     = cellfun(@(x, y) x * y, etaD, FdC, ...
        'UniformOutput', false);
    
    FS      = sum(cell2mat(FSC)) / Nfre;
    
    FD      = sum(cell2mat(FDC)) / Nfre;
    %     FC              = sum(cell2mat(FCC)) / Nfre;
    errOld  = err;
    err     = FS + FD;
    
    ChiErr  = norm(epsinv - chi{1}, 'fro') / norm(epsinv, 'fro');
    fprintf(logB, jj, FS, FD, FC, err, ChiErr, sb0, sb)
    
    if mod(jj, interval) == 0
        Mchi{kk} = mchi{end};
        kk = kk + 1;
    end
    if (abs(errOld - err) <= 1e-5 && terflag)
        tot = tot - 1;
        if tot == 0
            itenum = jj;
        end
    else
        tot = 5;
    end
    if (mod(jj, 32) == 0 && disflag)
        subplot(1, 2, 1); imagesc(real(mchi{end})); axis equal tight; colormap jet; colorbar; drawnow;
        subplot(1, 2, 2); imagesc(-imag(mchi{end})); axis equal tight; colormap jet; colorbar; drawnow;
    end
    jj = jj + 1;
    
end
if pars.recflag; diary off; end
time        = toc(startTime);
Mchi{end}   = mchi{end};
























