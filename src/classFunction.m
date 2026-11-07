function C = classFunction(psi,N1,N2)
% CLASSFUNCTION  Kulfan class function C(psi) = psi^N1 * (1-psi)^N2.
%
%   C = classFunction(psi,N1,N2)
%
% The class function fixes the family of the shape, the shape function
% (Bernstein polynomial, see bernsteinBasis) the particular member of the
% family (Kulfan 2008, Eq. 6 and Sec. III):
%   N1 = 0.5,   N2 = 1      round nose, sharp aft end (NACA type airfoil)
%   N1 = 0.5,   N2 = 0.5    elliptic airfoil, ellipsoid of revolution
%   N1 = 1,     N2 = 1      biconvex airfoil, ogive
%   N1 = 0.75,  N2 = 0.75   Sears-Haack body radius
%   N1 = 0.75,  N2 = 0.25   low drag projectile
%   N1 = 1,     N2 = 0.001  cone, wedge
%   N1 = 0.001, N2 = 0.001  rectangle, circular duct
%
% INPUTS
%   psi - non-dimensional coordinate in [0,1] (x/c for an airfoil), any size
%   N1  - exponent at psi = 0 (leading edge, nose)
%   N2  - exponent at psi = 1 (trailing edge, tail)
% OUTPUT
%   C   - class function, same size as psi
%
% See also: bernsteinBasis, cstCurve
arguments
    psi double {mustBeInRange(psi,0,1)}
    N1  (1,1) double {mustBeNonnegative}
    N2  (1,1) double {mustBeNonnegative}
end
C = psi.^N1.*(1-psi).^N2;   % 0^0 = 1 in MATLAB, so N = 0 gives a constant
end
