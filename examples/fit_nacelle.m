%%  Fit the CST parameters of a nacelle
% A nacelle is defined in the physical domain with Nacelle_class: five
% cowl profiles in the meridional planes from the crown (phi = 0) to the
% keel (phi = 180 deg), each one an airfoil (upper side = outer cowl, lower
% side = inner duct) with its highlight radius, length, x position and
% rotation (Kulfan 2008, Sec. X). The CST is fitted with Nacelle2CST and
% the nacelle rebuilt with CST2Nacelle. An axisymmetric nacelle is just a
% single profile (end of the PROCESS section).
clc; clear; close all
%% DATA -------------------------------------------------------------------
addpath(fullfile('..','src'))                         % matlab_cst functions

phi   = deg2rad([0 45 90 135 180])';                  % rad, crown -> keel
Rle   = [1.00 1.01 1.04 1.08 1.10]';                  % m, highlight radius
c     = [3.00 2.97 2.90 2.84 2.80]';                  % m, profile length
Xle   = [0.00 -0.02 -0.07 -0.12 -0.15]';              % m, highlight x (scarfed inlet)
twist = deg2rad([0 -0.5 -1 -1.6 -2])';                % rad, positive = lip outward
codes = ["4412" "4412" "4413" "4414" "4415"];         % NACA cowl profiles

% CST orders
nx      = 6;                                          % chordwise Bernstein order
ny      = 2;                                          % order in t = (1-cos(phi))/2 of the airfoil coefficients
ny_plan = 2;                                          % order in t of Rle, Xle, c, twist

% plot settings
axFS        = 15;
axFN        = "Times New Roman";
lw          = 2;
axlw        = 1.25;

%% PROCESS ----------------------------------------------------------------
Airfoil = arrayfun(@(s) Airfoil_class.naca4(s),codes,'UniformOutput',false);
Nac     = Nacelle_class(phi,Rle,c,twist,Airfoil,Xle=Xle,name="input nacelle");
CST     = Nacelle2CST(Nac,nx=nx,ny=ny,ny_plan=ny_plan);
NacC    = CST2Nacelle(CST,name="CST nacelle");

fprintf('max airfoil z/c error %.2e | max Rle error %.2e m | max Xle error %.2e m\n', ...
        max(CST.err.airfoil_rms),max(abs(CST.err.Rle)),max(abs(CST.err.Xle)));
fprintf('CST parameters: Bu, Bl %d x %d, profile quantities %d values each\n',size(CST.Bu),numel(CST.Rle));

% axisymmetric nacelle: one profile, ny = ny_plan = 0
NacAx  = Nacelle_class(0,1.0,3.0,0,Airfoil_class.naca4("4412"),name="axisymmetric");
CSTax  = Nacelle2CST(NacAx,nx=nx);
NacAxC = CST2Nacelle(CSTax);
fprintf('axisymmetric: ny = %d, ny_plan = %d, airfoil z/c error %.2e\n',CSTax.ny,CSTax.ny_plan,CSTax.err.airfoil_rms);

%% PLOT -------------------------------------------------------------------
set(0, 'defaultAxesTickLabelInterpreter','latex', ...
       'defaultTextInterpreter','latex', ...
       'DefaultLegendInterpreter','latex', ...
       'DefaultColorbarTickLabelInterpreter','latex');

% 3D: CST nacelle and input profiles
figure('Name','fit_nacelle_3D');
NacC.plot(gca,profiles=false);
psi = cosineSpacing(81);
for k = 1:numel(phi)
    [xo,ro,xi,ri] = Nac.profile(k,psi);
    r = [flip(ro); ri(2:end)];
    plot3([flip(xo); xi(2:end)],r*sin(phi(k)),r*cos(phi(k)),'r-','LineWidth',lw,'HandleVisibility','off');
end
view([-125 20]); axis off

% profiles: input vs CST
figure('Name','fit_nacelle_profiles');
Nac.plotProfiles(gca,'k-','LineWidth',lw);
NacP = CST2Nacelle(CST,phi=phi);                      % CST profiles at the input angles
NacP.plotProfiles(gca,'r--','LineWidth',1);
ylabel('$r$ [m]'); xlabel('$x$ [m]'); title('input (black) and CST (red dashed)');
styleAx(gca,axFS,axFN,axlw);

% axisymmetric nacelle
figure('Name','fit_nacelle_axisymmetric');
NacAxC.plot(gca); view([-125 20]); axis off

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
