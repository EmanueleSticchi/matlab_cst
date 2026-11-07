function CST = Airfoil2CST(AF,options)
% AIRFOIL2CST  Fit the CST coefficients of an airfoil (least squares).
%
%   CST = Airfoil2CST(AF)
%   CST = Airfoil2CST(AF,n=8,N1=0.5,N2=1,LEcontinuity=true,normalize=true)
%
% Each side of the airfoil is written as (Kulfan 2008, Eq. 17, 19)
%   z_u/c = C(psi) * sum_r Au_r S_r(psi) + psi*dzu
%   z_l/c = C(psi) * sum_r Al_r S_r(psi) + psi*dzl
% with psi = x/c, C the class function and S_r the Bernstein polynomials of
% order n. For fixed N1, N2 the unknowns [Au; Al; dzu; dzl] enter linearly,
% so the fit is a linear least squares problem solved with backslash: no
% optimizer, no initial guess, unique solution.
%
% INPUTS
%   AF  - Airfoil_class (or anything Airfoil_class accepts: [x z], file)
% OPTIONS (name-value)
%   n            - Bernstein polynomial order, n+1 coefficients per side
%                  (default 6: from about 6 on the airfoil is "statistically
%                  identical" to the original, Kulfan 2008 Sec. VII)
%   N1, N2       - class function exponents (default 0.5, 1: round nose,
%                  sharp trailing edge)
%   LEcontinuity - if true (default) impose Al(1) = -Au(1): the same
%                  leading edge radius on the two sides, so curvature is
%                  continuous around the nose
%   normalize    - if true (default) move the LE to (0,0) and the TE
%                  midpoint to (1,0) before the fit (see Airfoil_class.normalize)
% OUTPUT
%   CST - struct with the semantics of RASO/WiMOT/NUIST:
%     Au, Al   - (n+1) x 1 upper and lower coefficients
%     dzu, dzl - upper and lower z/c at the trailing edge (half TE thickness)
%     n, N1, N2
%     rLE      - leading edge radius r_LE/c = Au(1)^2/2 (Kulfan Eq. 4)
%     err      - struct: rms and max of the z/c error on the input points
%     AF       - the (normalized) airfoil that was fitted
%
% See also: CST2Airfoil, Airfoil_class, fitCstCurve
arguments
    AF
    options.n            (1,1) double {mustBeInteger,mustBeNonnegative} = 6
    options.N1           (1,1) double {mustBeNonnegative} = 0.5
    options.N2           (1,1) double {mustBeNonnegative} = 1
    options.LEcontinuity (1,1) logical = true
    options.normalize    (1,1) logical = true
end
if ~isa(AF,'Airfoil_class'); AF = Airfoil_class(AF); end
if options.normalize; AF = AF.normalize(); end
n  = options.n;
N1 = options.N1;
N2 = options.N2;

% data: the two sides from LE to TE, psi clamped to [0,1] (round-off of normalize)
[xu,zu,xl,zl] = AF.split();
psiu = min(max(xu,0),1);
psil = min(max(xl,0),1);

% linear system  M * [Au; Al; dzu; dzl] = [zu; zl]
Mu = classFunction(psiu,N1,N2).*bernsteinBasis(psiu,n);
Ml = classFunction(psil,N1,N2).*bernsteinBasis(psil,n);
nu = numel(psiu); nl = numel(psil); nA = n + 1;
M  = [Mu,           zeros(nu,nA), psiu,         zeros(nu,1); ...
      zeros(nl,nA), Ml,           zeros(nl,1),  psil];
T  = eye(2*nA + 2);                    % X = T*Y, Y = free unknowns
if options.LEcontinuity                % Al(1) = -Au(1): drop Al(1) from the unknowns
    T(nA+1,:) = 0;
    T(nA+1,1) = -1;
    T(:,nA+1) = [];
end
Y = (M*T) \ [zu; zl];
X = T*Y;

% output struct (same field names as RASO / WiMOT / NUIST CST)
CST.Au  = X(1:nA);
CST.Al  = X(nA+1:2*nA);
CST.dzu = X(end-1);
CST.dzl = X(end);
CST.n   = n;
CST.N1  = N1;
CST.N2  = N2;
CST.rLE = CST.Au(1)^2/2;
res     = [zu - cstCurve(psiu,CST.Au,CST.dzu,N1,N2); ...
           zl - cstCurve(psil,CST.Al,CST.dzl,N1,N2)];
CST.err = struct('rms',sqrt(mean(res.^2)),'max',max(abs(res)));
CST.AF  = AF;
end
