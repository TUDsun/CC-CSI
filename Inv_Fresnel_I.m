clearvars; restoredefaultpath; add_path

%% Load data, incident fields, measurement matrices, stiffness matrices

% str = 'rectTM_cent'; ptype = PT.TM; centre = [0, 0, 0]; xs = 40; ys = 40; m_unit = 1.0e-3; fre = 4 : 4 : 16;

% str = 'rectTM_dece'; ptype = PT.TM; centre = [0, 40, 0]; xs = 40; ys = 40; m_unit = 1.0e-3; fre = 4 : 4 : 16;

str = 'uTM_shaped'; para.pt = PT.TM; centre = [0, 0, 0]; xs = 70; ys = 70; para.m_unit = 1.5e-3; fre = 2 : 2 : 16;

% str = 'rectTE_8f'; ptype = PT.TE; centre = [0, 0, 0]; xs = 40; ys = 40; m_unit = 1.0e-3; fre = 4 : 4 : 16;

% str = 'dielTM_dec8f'; ptype = PT.TM;

% str = 'twodielTM_8f'; ptype = PT.TM; centre = [0, 0, 0]; xs = 70; ys = 100; m_unit = 1.5e-3; fre = 2 : 2 : 8;

regSize         = [-1.2, 1.2; -1.2, 1.2]; % 12GHz
para.NTX        = 36;
TxInterval      = 360 / para.NTX;
% rho             = 0.20066;
% DTr             = rho/0.6964*ones(1,18);
DTr             = 0; % zeros(1, 18);
para.Tr      	= 0.72 + DTr / 2; % para.Tr      	= 0.72 * ones(1, 18) + DTr / 2;
para.Rr       	= 0.76 + DTr / 2;
rawdat          = load([str '.txt']);
rawdat(:, 1)    = TxInterval * (rawdat(:, 1) - 1) - 2.5;
rawdat(:, 2)    = 5 * (rawdat(:, 2) - 1);

%% Configuration

[Emea, Phi, A, chiinv, vJ, Einc, eTotInv, grid3d, pars] ...
    = Pre_InvFresnel_Conf(rawdat, fre, regSize, centre, xs, ys, para);

%% Display

run show_InitGuess.m

%% Iteration starts

Phi                 = cellfun(@(x, y) x * y, Phi, pars.omega, 'UniformOutput', false);
rf                  = 1 : length(fre);      % select the frequency components for inversion
pars.bgchi          = 1;
pars.disflag        = true;                 % show the output every N iterations?
pars.terflag        = false;                %
pars.recflag        = false;                % save the records?
pars.omega          = pars.omega(rf);
pars.itenum        	= 16;                   % set iteration number

pars.recstr	= ['RecMF' str '_CCCSI'];

[mchi1, time1, chi1, vJ1, eTot1] = CC_CSI_Brent(...
    Emea(rf), Phi(rf), A(rf), chiinv(rf), vJ(rf), Einc(rf), eTotInv(rf), pars, []);

pars.recstr	= ['RecMF' str '_MRCSI'];

[mchi2, time2, chi2, vJ2, eTot2] = MR_CSI_Brent(...
    Emea(rf), Phi(rf), A(rf), chiinv(rf), vJ(rf), Einc(rf), eTotInv(rf), pars, []);

pars.recstr	= ['RecMF' str '_CSI'];

[mchi3, time3, chi3, vJ3, eTot3] = CSI(...
    Emea(rf), Phi(rf), A(rf), chiinv(rf), vJ(rf), Einc(rf), eTotInv(rf), pars, []);

%% save the inversion output

% save(['OMF' str  '.mat'], 'xinv', 'yinv', 'invdom', 'omega',...
%     'mchi1', 'time1', 'chi1', 'vJ1', 'eTot1', ...
%     'mchi2', 'time2', 'chi2', 'vJ2', 'eTot2', ...
%     'mchi3', 'time3', 'chi3', 'vJ3', 'eTot3', ...
%     'grid3d', 'ii', 'pars')

%% Show the inversion output

run show_output.m

