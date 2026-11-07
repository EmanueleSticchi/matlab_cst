%%  Fit the CST parameters of an airfoil
% Reads an airfoil file (Selig or Lednicer), fits the CST coefficients with
% Airfoil2CST for increasing Bernstein order n and shows the fit, the
% residuals and the convergence with n (as Kulfan 2008, Fig. 9-11).
clc; clear; close all
%% DATA -------------------------------------------------------------------
addpath(fullfile('..','src'))                         % matlab_cst functions

airfoilFile  = fullfile('data','naca0012.dat');       % any Selig or Lednicer file
% airfoilFile = Airfoil_class.naca4("2412");          % ... or an Airfoil_class
nList        = 2:12;                                  % Bernstein orders of the convergence study
nShow        = 6;                                     % order shown in detail
N1           = 0.5;                                   % class function: round nose
N2           = 1.0;                                   %                 sharp trailing edge
LEcontinuity = true;                                  % Al(1) = -Au(1)
saveflag     = false;                                 % write the fitted airfoil and its CST

% plot settings
axFS        = 15;
axFN        = "Times New Roman";
lw          = 2;
axlw        = 1.25;

%% PROCESS ----------------------------------------------------------------
AF = Airfoil_class(airfoilFile);

% convergence with the Bernstein order
rmsErr = zeros(size(nList)); maxErr = rmsErr;
for i = 1:numel(nList)
    CST       = Airfoil2CST(AF,n=nList(i),N1=N1,N2=N2,LEcontinuity=LEcontinuity);
    rmsErr(i) = CST.err.rms;
    maxErr(i) = CST.err.max;
end

% detail at n = nShow
CST   = Airfoil2CST(AF,n=nShow,N1=N1,N2=N2,LEcontinuity=LEcontinuity);
AFfit = CST2Airfoil(CST,201,name=AF.name + " CST n = " + nShow);
fprintf('%s, n = %d\n',AF.name,nShow);
fprintf('  Au  = %s\n  Al  = %s\n',mat2str(CST.Au',5),mat2str(CST.Al',5));
fprintf('  dzu = %.6f, dzl = %.6f, r_LE/c = %.5f\n',CST.dzu,CST.dzl,CST.rLE);
fprintf('  z/c error: rms %.2e, max %.2e\n',CST.err.rms,CST.err.max);

% residuals at the data points (the fitted airfoil is the normalized one)
[xu,zu,xl,zl] = CST.AF.split();
resU = zu - cstCurve(min(max(xu,0),1),CST.Au,CST.dzu,N1,N2);
resL = zl - cstCurve(min(max(xl,0),1),CST.Al,CST.dzl,N1,N2);

if saveflag
    AFfit.write('airfoil_cst.dat');
    save('airfoil_cst.mat','CST');
end

%% PLOT -------------------------------------------------------------------
set(0, 'defaultAxesTickLabelInterpreter','latex', ...
       'defaultTextInterpreter','latex', ...
       'DefaultLegendInterpreter','latex', ...
       'DefaultColorbarTickLabelInterpreter','latex');

% airfoil and fit
figure('Name','fit_airfoil_shape');
plot(CST.AF.x,CST.AF.z,'ko','MarkerSize',4,'DisplayName','data'); hold on
AFfit.plot(gca,'r-','LineWidth',lw,'DisplayName',sprintf('CST, $n = %d$',nShow));
xlabel('$x/c$'); ylabel('$z/c$'); legend('Location','northeast');
styleAx(gca,axFS,axFN,axlw);

% residuals: lower side on -x/c, upper side on x/c (as Kulfan Fig. 9)
figure('Name','fit_airfoil_residual');
plot(-xl,resL*1e4,'r.-','DisplayName','lower'); hold on
plot( xu,resU*1e4,'b.-','DisplayName','upper');
yline(0,'k-','HandleVisibility','off');
xlabel('$-x/c$ (lower) $\quad x/c$ (upper)'); ylabel('$\Delta z/c \times 10^4$');
legend('Location','best'); styleAx(gca,axFS,axFN,axlw);

% convergence
figure('Name','fit_airfoil_convergence');
semilogy(nList,rmsErr,'ko-','LineWidth',lw,'DisplayName','rms'); hold on
semilogy(nList,maxErr,'rs-','LineWidth',lw,'DisplayName','max');
xlabel('Bernstein order $n$'); ylabel('$|\Delta z/c|$');
legend('Location','northeast'); styleAx(gca,axFS,axFN,axlw);

%% FUNCTIONS --------------------------------------------------------------
function styleAx(ax,axFS,axFN,axlw)
% axes style of the script skeleton
ax.FontName       = axFN;
ax.FontSize       = axFS;
ax.Color          = [1,1,1];
ax.LineWidth      = axlw;
ax.XMinorTick     = "on";
ax.YMinorTick     = "on";
ax.Layer          = "top";
ax.GridAlpha      = 0.1;
ax.MinorGridAlpha = 0.1;
box(ax,'on'); grid(ax,'on');
end
