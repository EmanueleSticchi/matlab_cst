function [A,dz,err] = fitCstCurve(psi,zeta,n,N1,N2,options)
% FITCSTCURVE  Least squares CST coefficients of one curve.
%
%   [A,dz,err] = fitCstCurve(psi,zeta,n,N1,N2)
%   [A,dz,err] = fitCstCurve(psi,zeta,n,N1,N2,dz=0)    % end value imposed
%
% For fixed N1, N2 the CST curve is LINEAR in the coefficients A and in dz:
%   zeta = [C.*B, psi] * [A; dz]
% so the best fit is a plain linear least squares problem (backslash), no
% optimizer and no initial guess are needed and the solution is unique.
%
% INPUTS
%   psi  - points in [0,1] (x/c)
%   zeta - curve values at psi (z/c)
%   n    - order of the Bernstein polynomial (n+1 coefficients)
%   N1   - class function exponent at psi = 0
%   N2   - class function exponent at psi = 1
%   dz   - (name-value) value at psi = 1; empty (default) = fitted
% OUTPUTS
%   A    - (n+1) x 1 coefficients
%   dz   - value at psi = 1
%   err  - residual zeta - cstCurve(psi,A,dz,N1,N2)
%
% See also: cstCurve, Airfoil2CST
arguments
    psi  (:,1) double {mustBeInRange(psi,0,1)}
    zeta (:,1) double {mustBeFinite}
    n    (1,1) double {mustBeInteger,mustBeNonnegative}
    N1   (1,1) double {mustBeNonnegative} = 0.5
    N2   (1,1) double {mustBeNonnegative} = 1
    options.dz double = []
end
if numel(psi) ~= numel(zeta)
    error('\nfitCstCurve: psi and zeta must have the same number of points (%d vs %d)\n',numel(psi),numel(zeta));
end
M = classFunction(psi,N1,N2).*bernsteinBasis(psi,n);   % columns: C(psi)*S_r(psi)
if isempty(options.dz)
    X  = [M, psi] \ zeta;
    A  = X(1:end-1);
    dz = X(end);
else
    dz = options.dz;
    A  = M \ (zeta - psi*dz);
end
err = zeta - cstCurve(psi,A,dz,N1,N2);
end
