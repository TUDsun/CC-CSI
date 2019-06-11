function [Emea, Phi, A, chiinv, vJ, Einc, eTotInv, grid3d, pars] = Pre_InvFresnel_Conf(rawdat, fre, regSize, centre, xs, ys, para)
%%
NP              = 1 : 3;
regSize         = round(regSize / para.m_unit) * para.m_unit;
xs              = round(xs * 1e-3 / para.m_unit) * para.m_unit;
ys              = round(ys * 1e-3 / para.m_unit) * para.m_unit;
centre          = round(centre * 1e-3 / para.m_unit) * para.m_unit;
invdom(1)       = centre(1) - xs;
invdom(2)       = centre(1) + xs;
invdom(3)       = centre(2) - ys;
invdom(4)       = centre(2) + ys;
invdom(5)       = centre(3);
invdom(6)       = centre(3);

Nfre            = length(fre);
A               = cell(1, Nfre);
dat             = cell(1, Nfre);
Emea            = cell(1, Nfre);
Einc            = cell(1, Nfre);
err             = cell(1, Nfre);
grid3d          = cell(1, Nfre);
omega           = cell(1, Nfre);
Phi             = cell(1, Nfre);
rmea            = cell(1, Nfre);
XBP             = cell(1, Nfre);

para.disflag    = false;
para.showconfg  = false;

% Tr              = Tr(fre);
% Rr              = Rr(fre);
for ii = 1 : Nfre
%     para.Tr     = Tr(ii);
%     para.Rr     = Rr(ii);
    
    [dat{ii}, Phi{ii}, A{ii}, Einc{ii}, grid3d{ii}, Rindex, Runiq, pars] ...
        = ReadFreData(fre(ii), rawdat, invdom, regSize, para);
    
    omega{ii}   = pars.omega;
    
    r           = Rindex + repmat(...
        (0 : (size(dat{ii}, 2) - 1)) * length(Runiq) * NP(para.pt), ...
        size(dat{ii}, 1), 1);
    
    rmea{ii}    = r(:);
    tmp         = zeros(size(Phi{ii}, 1), size(dat{ii}, 2));
    tmp(r(:))   = dat{ii};
    Emea{ii}    = tmp;
    
    % Backpropagation, multiplied by a weight to ensure that 
    % the data error is minimized ...
    phimea      = Phi{ii}' * tmp;
    pphimea     = Phi{ii} * phimea;
    kappa       = sum(abs(phimea) .^ 2) ./ sum(abs(pphimea) .^ 2);
    XBP{ii}     = phimea * diag(kappa);
    err{ii}     = tmp - Phi{ii} * XBP{ii};
    %     norm(err{ii}, 'fro')./norm(tmp, 'fro')
end

XBP	= cell2mat(XBP);

%% Linear sampling method
% ILSM            = zeros(1, size(Phi{1}, 2));
% for ii = 1 : Nfre
%     [U, S, V]   = svd(Emea{ii});
%     Spse        = S';
%     Spse        = Spse ./ (Spse .^ 2 + (0.01 * S(1)) ^ 2);
%     ILSM        = ILSM + sum(abs(Spse * U' * Phi{ii}) .^ 2, 1);
% end
% ILSM            = ILSM(:, 1 : NP(ptype) : end) + ILSM(:, NP(ptype) : NP(ptype) : end);
% ILSM            = reshape(1 ./ ILSM, pars.Ninv(1), pars.Ninv(2));
% ILSM            = ILSM / max(ILSM(:));
% ILSMdB          = db(ILSM, 'power');

%% vJ, chiinv, eTotInv

[NRX, NTX]    	= size(dat{1});
bgchi        	= 1;
vJ           	= cell(1, Nfre);
eTotInv        	= cell(1, Nfre);
chiinv        	= cell(1, Nfre);
for ii = 1 : Nfre
    sind                = NTX * (ii - 1) + 1 : NTX * ii;
    vJ{ii}              = XBP(:, sind) / omega{ii};
    vJ{ii}(pars.rd, :)  = 0;
    eSctInv             = A{ii} \ vJ{ii};
    eTotInv{ii}         = Einc{ii} + eSctInv;
end
chiinv{1}    	= get_chiMF(vJ, eTotInv, pars.Ninv, NTX, omega, para.pt);
chiinv{1}    	= OrthgPro(chiinv{1}, bgchi);
chiinv{1}(...
    pars.rchi)  = 0;

chiinv          = cellfun(@(x) ...
    real(chiinv{1}) + 1j * imag(chiinv{1}) * omega{1} / x, omega, ...
    'UniformOutput', false);

pars.rmea      	= rmea;
pars.omega      = omega;
pars.ptype      = para.pt;