%%  Fit the CST parameters of a fuselage
% An airliner-like fuselage is given as cross sections (point clouds in the
% planes x = const, as they would be cut from a CAD model). The
% longitudinal lines and the lobe exponents are extracted with
% Fuselage_class.fromSections, the CST is fitted with Fuselage2CST for
% increasing Bernstein order and the fuselage rebuilt with CST2Fuselage.
clc; clear; close all
%% DATA -------------------------------------------------------------------
addpath(fullfile('..','src'))                         % matlab_cst functions

% synthetic cross sections
L        = 40;                                        % m, length
R        = 2;                                         % m, cabin radius
xNose    = 6;                                         % m, end of the nose
xTail    = 26;                                        % m, start of the tail cone
upsweep  = 1.5;                                       % m, tail centre height
Nst      = 60;                                        % number of sections
Ntheta   = 61;                                        % points per section

% CST orders
nList    = [4 6 8 10 12];                             % Bernstein orders of hU, hL, w
nz       = 6;                                         % order of the max half-breadth line
nN       = 3;                                         % order of the lobe exponents

% plot settings
axFS        = 15;
axFN        = "Times New Roman";
lw          = 2;
axlw        = 1.25;

%% PROCESS ----------------------------------------------------------------
% sections: elliptic nose, cylinder, upswept tail cone; the lower lobe of the
% nose is flatter (NL < 0.5, cockpit-like) and the lower lobe is 10% shallower
x  = L*cosineSpacing(Nst);
th = linspace(0,2*pi,Ntheta);
Ysec = cell(Nst,1); Zsec = cell(Nst,1);
for k = 1:Nst
    [r,zm,NL] = sectionLaw(x(k),L,R,xNose,xTail,upsweep);
    up      = sin(th) >= 0;                           % upper lobe: elliptic (NU = 0.5)
    Ysec{k} = r*cos(th);
    Zsec{k} = zeros(size(th));
    Zsec{k}(up)  = zm + r*(sin(th(up)).^2).^0.5;
    Zsec{k}(~up) = zm - 0.9*r*(sin(th(~up)).^2).^NL;
end
F = Fuselage_class.fromSections(x,Ysec,Zsec,name="input fuselage");

% CST for increasing order
errMax = zeros(numel(nList),3);
for i = 1:numel(nList)
    CST = Fuselage2CST(F,n=nList(i),nz=nz,nN=nN);
    errMax(i,:) = [max(abs(CST.err.hU)) max(abs(CST.err.hL)) max(abs(CST.err.w))];
    fprintf('n = %2d: max error hU %.3f m, hL %.3f m, w %.3f m, zm %.3f m\n',nList(i),errMax(i,:),max(abs(CST.err.zm)));
end
FC = CST2Fuselage(CST,name="CST fuselage");          % last (highest) order

%% PLOT -------------------------------------------------------------------
set(0, 'defaultAxesTickLabelInterpreter','latex', ...
       'defaultTextInterpreter','latex', ...
       'DefaultLegendInterpreter','latex', ...
       'DefaultColorbarTickLabelInterpreter','latex');

% 3D
figure('Name','fit_fuselage_3D');
FC.plot(gca,sections=false);
for k = 1:4:Nst
    plot3(x(k)*ones(size(th)),Ysec{k},Zsec{k},'r-','LineWidth',1,'HandleVisibility','off');
end
view([-35 20]); axis off

% longitudinal lines
figure('Name','fit_fuselage_lines');
F.plotLines(gca,'Color','k','LineWidth',lw);
FC.plotLines(gca,'Color','r','LineWidth',1);
xlabel('$x$ [m]'); ylabel('$z$, $w$ [m]'); title('input (black) and CST (red)');
styleAx(gca,axFS,axFN,axlw);

% convergence
figure('Name','fit_fuselage_convergence');
semilogy(nList,errMax,'o-','LineWidth',lw);
legend('$h_U$','$h_L$','$w$'); xlabel('Bernstein order $n$'); ylabel('max error [m]');
styleAx(gca,axFS,axFN,axlw);

%% FUNCTIONS --------------------------------------------------------------
function [r,zm,NL] = sectionLaw(x,L,R,xNose,xTail,upsweep)
% radius, centre height and lower lobe exponent of the synthetic fuselage
if x < xNose                                          % elliptic nose
    r  = R*sqrt(1 - ((xNose - x)/xNose)^2);
    zm = -0.3*((xNose - x)/xNose)^2;
    NL = 0.5 - 0.2*(xNose - x)/xNose;
elseif x <= xTail                                     % cylinder
    r = R; zm = 0; NL = 0.5;
else                                                  % tail cone
    s  = (x - xTail)/(L - xTail);
    r  = R*(1 - 0.9*s^1.5);
    zm = upsweep*s^2;
    NL = 0.5;
end
end

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
