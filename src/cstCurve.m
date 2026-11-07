function zeta = cstCurve(psi,A,dz,N1,N2)
% CSTCURVE  One CST curve (one side of an airfoil, a radius distribution, ...).
%
%   zeta = cstCurve(psi,A,dz,N1,N2)
%
%   zeta(psi) = C(psi) * S(psi) + psi*dz              (Kulfan 2008, Eq. 7, 17)
%   C(psi)    = psi^N1 (1-psi)^N2                     class function
%   S(psi)    = sum_r A_r K_r psi^r (1-psi)^(n-r)     shape function
%
% INPUTS
%   psi - non-dimensional coordinate in [0,1] (x/c)
%   A   - n+1 Bernstein coefficients (n = order of the polynomial)
%   dz  - value at psi = 1 (half trailing edge thickness for an airfoil side)
%   N1  - class function exponent at psi = 0 (0.5: round nose)
%   N2  - class function exponent at psi = 1 (1: sharp aft end)
% OUTPUT
%   zeta - column, same number of points as psi (z/c for an airfoil)
%
% Physical meaning of the end coefficients (round nose, N1 = 0.5, N2 = 1):
%   A(1)   = sqrt(2*R_LE/c)         leading edge radius
%   A(end) = tan(beta) + dz         beta = boat-tail angle
%
% See also: classFunction, bernsteinBasis, fitCstCurve, CST2Airfoil
arguments
    psi (:,1) double {mustBeInRange(psi,0,1)}
    A   (:,1) double {mustBeFinite}
    dz  (1,1) double {mustBeFinite} = 0
    N1  (1,1) double {mustBeNonnegative} = 0.5
    N2  (1,1) double {mustBeNonnegative} = 1
end
n    = numel(A) - 1;
zeta = classFunction(psi,N1,N2).*(bernsteinBasis(psi,n)*A) + psi*dz;
end
