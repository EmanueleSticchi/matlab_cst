%%  Fit the CST parameters of a cranked wing
% A wing is defined in the physical domain with Wing3D_class (the NUIST
% Wing_class arguments plus the leading edge position), its spanwise CST is
% fitted with Wing2CST and the wing is rebuilt with CST2Wing.
% The same wing is fitted with and without the planform break, to show
% what the break buys: without it a straight panel cannot follow the kink
% of the leading and trailing edge.
clc; clear; close all
%% DATA -------------------------------------------------------------------
addpath(fullfile('..','src'))                         % matlab_cst functions

% wing geometry (half wing), NUIST-like inputs
eta     = [0 0.15 0.35 0.55 0.75 0.9 1]';             % y/b of the sections
b       = 18;                                         % m, half span
c       = interp1([0 0.35 1],[7 4.5 1.6],eta);        % m, chord: kink at the break
twist   = deg2rad([3.0 2.4 1.0 0.2 -0.6 -1.1 -1.5]'); % rad, nose up about c/4
Xle     = interp1([0 0.35 1],[0 3 9.5],eta);          % m, untwisted LE: kink at the break
Zle     = interp1([0 0.35 1],[0 0.3 1.2],eta);        % m, dihedral
codes   = ["2415" "2414" "2413" "2412" "2411" "2410" "2410"];   % NACA airfoils
eta_break = 0.35;                                     % planform break (kink)

% CST orders
nx      = 8;                                          % chordwise Bernstein order
ny      = 2;                                          % spanwise order of the airfoil coefficients (per segment)
ny_plan = 1;                                          % spanwise order of c, twist, Xle, Zle (per segment)
saveflag = false;                                     % write the CST and a Tecplot file of the wing

% plot settings
axFS        = 15;
axFN        = "Times New Roman";
lw          = 2;
axlw        = 1.25;

%% PROCESS ----------------------------------------------------------------
Airfoil = arrayfun(@(s) Airfoil_class.naca4(s),codes,'UniformOutput',false);
Wing    = Wing3D_class(eta,b,c,twist,Airfoil,Xle=Xle,Zle=Zle,eta_breaks=eta_break,name="input wing");
% the same wing from 3D sections (e.g. cut from a CAD model) would be
% Wing = Wing3D_class.fromSections(Y,XZ,eta_breaks=eta_break);

% fit with and without the break
CST   = Wing2CST(Wing,nx=nx,ny=ny,ny_plan=ny_plan);
CST0  = Wing2CST(Wing,nx=nx,ny=ny,ny_plan=ny_plan,eta_breaks=[]);
WingC = CST2Wing(CST,name="CST wing");
WingC0 = CST2Wing(CST0,name="CST wing, no break");

fprintf('                         with break   without break\n');
fprintf('max airfoil z/c error    %10.2e   %10.2e\n',max(CST.err.airfoil_rms),max(CST0.err.airfoil_rms));
fprintf('max chord error [m]      %10.2e   %10.2e\n',max(abs(CST.err.c)),max(abs(CST0.err.c)));
fprintf('max Xle error [m]        %10.2e   %10.2e\n',max(abs(CST.err.Xle)),max(abs(CST0.err.Xle)));
fprintf('max twist error [deg]    %10.2e   %10.2e\n',rad2deg(max(abs(CST.err.twist))),rad2deg(max(abs(CST0.err.twist))));
fprintf('CST parameters: Bu, Bl %d x %d, planform %d values each\n',size(CST.Bu),numel(CST.c));
fprintf('planform control values (root, break, tip): c = %s m, Xle = %s m\n',mat2str(CST.c',4),mat2str(CST.Xle',4));

if saveflag
    save('wing_cst.mat','CST');
    writeTecplot('wing_cst.dat',WingC.surface(),title="CST wing");
end

%% PLOT -------------------------------------------------------------------
set(0, 'defaultAxesTickLabelInterpreter','latex', ...
       'defaultTextInterpreter','latex', ...
       'DefaultLegendInterpreter','latex', ...
       'DefaultColorbarTickLabelInterpreter','latex');

% 3D: CST wing and input sections
figure('Name','fit_wing_3D');
WingC.plot(gca,sections=false);
for k = 1:numel(eta)                                  % input sections on top
    [Xu,Zu,Xl,Zl,Y] = Wing.section3D(k,cosineSpacing(61));
    plot3([flip(Xu) Xl],[flip(Y) Y],[flip(Zu) Zl],'r-','LineWidth',lw,'HandleVisibility','off');
end
view([-140 25]); axis off

% planform with and without the break
figure('Name','fit_wing_planform');
Wing.plotPlanform(gca,'ko-','LineWidth',lw,'DisplayName','input');
WingC.plotPlanform(gca,'b-','LineWidth',1,'DisplayName','CST, break at $\eta = 0.35$');
WingC0.plotPlanform(gca,'r--','LineWidth',1,'DisplayName','CST, no break');
legend('Location','southoutside'); styleAx(gca,axFS,axFN,axlw);

% airfoil error along the span
figure('Name','fit_wing_error');
semilogy(eta,CST.err.airfoil_rms,'bo-','LineWidth',lw,'DisplayName','with break'); hold on
semilogy(eta,CST0.err.airfoil_rms,'rs--','LineWidth',lw,'DisplayName','without break');
xlabel('$\eta = y/b$'); ylabel('rms $|\Delta z/c|$');
legend('Location','best'); styleAx(gca,axFS,axFN,axlw);

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
