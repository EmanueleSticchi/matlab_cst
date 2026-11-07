function CST = Wing2CST(Wing,options)
% WING2CST  Fit the spanwise CST parameters of a wing (least squares).
%
%   CST = Wing2CST(Wing)
%   CST = Wing2CST(Wing,nx=6,ny=2,ny_plan=1,eta_breaks=Wing.eta_breaks)
%
% Kulfan (2008, Sec. XII-XIII) "scalar loft": every chordwise CST
% coefficient of the sections is a Bernstein polynomial along the span,
%   Au_i(eta) = sum_j Bu_ij Sy_j(eta)          (Eq. 43, the same for Al)
% and so is every planform quantity (chord, twist, LE position, TE
% thickness). With planform breaks the spanwise polynomials are piecewise,
% one per segment, continuous at the breaks (see spanwiseBasis).
%
% The fit is done in two linear steps, after normalizing the section
% airfoils (normalizeSections: LE at (0,0), TE at (1,0), any residual
% chord, rotation or LE shift moved to c, twist, Xle, Zle; the physical
% wing does not change):
%   1) each section airfoil -> chordwise coefficients (fitSections)
%   2) coefficients and planform along the span -> control values
%      (fitSpanwise)
%
% INPUTS
%   Wing - Wing3D_class
% OPTIONS (name-value)
%   nx           - chordwise Bernstein order (default 6)
%   ny           - spanwise order of the airfoil coefficients in each
%                  segment ([] = automatic, up to 3, limited by the number
%                  of sections in each segment)
%   ny_plan      - spanwise order of c, twist, Xle, Zle, dzu, dzl in each
%                  segment (default 1: straight leading and trailing edge
%                  and linear twist in each segment, the usual wing panel)
%   eta_breaks   - planform breaks in (0,1) (default Wing.eta_breaks)
%   N1, N2       - class function exponents (default 0.5, 1)
%   LEcontinuity - impose Al(1) = -Au(1) in every section (default true)
% OUTPUT
%   CST - struct
%     type          - "wing"
%     b, b0         - span and root y, m (from Wing)
%     eta_breaks    - breaks
%     nx, ny, ny_plan, N1, N2, twist_axis
%     Bu, Bl        - (Nseg*ny+1) x (nx+1): row j = spanwise control value j,
%                     column i = chordwise coefficient i
%     c, twist, Xle, Zle, dzu, dzl - (Nseg*ny_plan+1) x 1 control values
%                     (m, rad). With ny_plan = 1 they ARE the values at the
%                     root, at the breaks and at the tip.
%     eta_ctrl, eta_ctrl_plan - eta of the control values
%     err           - struct, error of the CST wing at the input sections:
%                     airfoil_rms, airfoil_max (N x 1, z/c), c, twist, Xle,
%                     Zle (N x 1, CST - input after normalizeSections)
%
% See also: CST2Wing, Wing3D_class, sectionCST, spanwiseBasis
arguments
    Wing {mustBeA(Wing,'Wing3D_class')}
    options.nx           (1,1) double {mustBeInteger,mustBeNonnegative} = 6
    options.ny           double = []
    options.ny_plan      (1,1) double {mustBeInteger,mustBeNonnegative} = 1
    options.eta_breaks   (1,:) double = Wing.eta_breaks
    options.N1           (1,1) double {mustBeNonnegative} = 0.5
    options.N2           (1,1) double {mustBeNonnegative} = 1
    options.LEcontinuity (1,1) logical = true
end
Wing.check();
Wing   = Wing.normalizeSections();
eta    = Wing.eta;
breaks = options.eta_breaks;

% 1) chordwise CST of every section
[Au,Al,dzu,dzl] = fitSections(Wing.Airfoil,options.nx,options.N1,options.N2,options.LEcontinuity);

% 2) spanwise control values
[Bu,ny] = fitSpanwise(eta,Au,breaks,options.ny);
Bl      = fitSpanwise(eta,Al,breaks,ny);
plan    = fitSpanwise(eta,[Wing.c Wing.twist Wing.Xle Wing.Zle dzu dzl],breaks,options.ny_plan);

% output
CST.type          = "wing";
CST.b             = Wing.b;
CST.b0            = Wing.b0;
CST.eta_breaks    = breaks;
CST.nx            = options.nx;
CST.ny            = ny;
CST.ny_plan       = options.ny_plan;
CST.N1            = options.N1;
CST.N2            = options.N2;
CST.twist_axis    = Wing.twist_axis;
CST.Bu            = Bu;
CST.Bl            = Bl;
CST.c             = plan(:,1);
CST.twist         = plan(:,2);
CST.Xle           = plan(:,3);
CST.Zle           = plan(:,4);
CST.dzu           = plan(:,5);
CST.dzl           = plan(:,6);
[~,CST.eta_ctrl]      = spanwiseBasis(0,breaks,ny);
[~,CST.eta_ctrl_plan] = spanwiseBasis(0,breaks,options.ny_plan);

% error of the CST wing at the input sections
sec   = sectionCST(CST,eta);
Bp    = spanwiseBasis(eta,breaks,CST.ny_plan);
[CST.err.airfoil_rms,CST.err.airfoil_max] = sectionError(Wing.Airfoil,sec);
CST.err.c     = Bp*CST.c     - Wing.c;
CST.err.twist = Bp*CST.twist - Wing.twist;
CST.err.Xle   = Bp*CST.Xle   - Wing.Xle;
CST.err.Zle   = Bp*CST.Zle   - Wing.Zle;
end
