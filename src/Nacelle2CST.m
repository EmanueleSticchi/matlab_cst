function CST = Nacelle2CST(Nac,options)
% NACELLE2CST  Fit the CST parameters of a nacelle (least squares).
%
%   CST = Nacelle2CST(Nac)
%   CST = Nacelle2CST(Nac,nx=6,ny=2,ny_plan=2)
%
% Same scalar loft as a wing (Wing2CST), with the meridional angle phi in
% place of the span. Every chordwise coefficient and every profile
% quantity (Rle, Xle, c, twist, TE thickness) is a Bernstein polynomial in
%   t = (1 - cos(phi))/2,     t = 0 crown, t = 1/2 side, t = 1 keel
% A smooth function of cos(phi) is even about phi = 0 and phi = pi, so the
% nacelle surface is smooth across the symmetry plane (no ridge on the
% crown and keel lines), whatever the coefficients.
%
% INPUTS
%   Nac - Nacelle_class
% OPTIONS (name-value)
%   nx           - chordwise Bernstein order (default 6)
%   ny           - order in t of the airfoil coefficients ([] = automatic,
%                  up to 2: crown, side, keel)
%   ny_plan      - order in t of Rle, Xle, c, twist, dzu, dzl ([] =
%                  automatic, up to 2)
%   N1, N2       - class function exponents (default 0.5, 1)
%   LEcontinuity - impose Al(1) = -Au(1) (default true)
% OUTPUT
%   CST - struct
%     type             - "nacelle"
%     nx, ny, ny_plan, N1, N2, twist_axis, eta_breaks (= [], no breaks)
%     Bu, Bl           - (ny+1) x (nx+1) control values (rows: crown -> keel)
%     Rle, Xle, c, twist, dzu, dzl - (ny_plan+1) x 1 control values
%     err              - struct, error at the input profiles: airfoil_rms,
%                        airfoil_max (z/c), Rle, Xle, c, twist (CST - input)
%   An axisymmetric nacelle (one profile) gives ny = ny_plan = 0.
%
% See also: CST2Nacelle, Nacelle_class, Wing2CST, sectionCST
arguments
    Nac {mustBeA(Nac,'Nacelle_class')}
    options.nx           (1,1) double {mustBeInteger,mustBeNonnegative} = 6
    options.ny           double = []
    options.ny_plan      double = []
    options.N1           (1,1) double {mustBeNonnegative} = 0.5
    options.N2           (1,1) double {mustBeNonnegative} = 1
    options.LEcontinuity (1,1) logical = true
end
Nac.check();
Nac = Nac.normalizeSections();
t   = (1 - cos(Nac.phi))/2;
nob = zeros(1,0);                                    % no breaks around the nacelle

% 1) chordwise CST of every profile
[Au,Al,dzu,dzl] = fitSections(Nac.Airfoil,options.nx,options.N1,options.N2,options.LEcontinuity);

% 2) control values in t
[Bu,ny]        = fitSpanwise(t,Au,nob,options.ny,maxOrder=2);
Bl             = fitSpanwise(t,Al,nob,ny);
[plan,ny_plan] = fitSpanwise(t,[Nac.Rle Nac.Xle Nac.c Nac.twist dzu dzl],nob,options.ny_plan,maxOrder=2);

% output
CST.type       = "nacelle";
CST.eta_breaks = nob;
CST.nx         = options.nx;
CST.ny         = ny;
CST.ny_plan    = ny_plan;
CST.N1         = options.N1;
CST.N2         = options.N2;
CST.twist_axis = Nac.twist_axis;
CST.Bu         = Bu;
CST.Bl         = Bl;
CST.Rle        = plan(:,1);
CST.Xle        = plan(:,2);
CST.c          = plan(:,3);
CST.twist      = plan(:,4);
CST.dzu        = plan(:,5);
CST.dzl        = plan(:,6);

% error at the input profiles
sec = sectionCST(CST,t);
Bp  = spanwiseBasis(t,nob,ny_plan);
[CST.err.airfoil_rms,CST.err.airfoil_max] = sectionError(Nac.Airfoil,sec);
CST.err.Rle   = Bp*CST.Rle   - Nac.Rle;
CST.err.Xle   = Bp*CST.Xle   - Nac.Xle;
CST.err.c     = Bp*CST.c     - Nac.c;
CST.err.twist = Bp*CST.twist - Nac.twist;
end
