function B = bernsteinBasis(t,n)
% BERNSTEINBASIS  Bernstein polynomials of order n evaluated at t.
%
%   B = bernsteinBasis(t,n)
%
%   B(:,r+1) = K_r * t^r * (1-t)^(n-r),   K_r = n!/(r!(n-r)!),   r = 0..n
%
% The n+1 polynomials sum to 1 at every t ("partition of unity", Kulfan
% 2008, Sec. V), so each one is a scalable piece of the unit shape function.
% Only the first is non zero at t = 0 and only the last at t = 1: for an
% airfoil the first coefficient sets the leading edge radius and the last
% one the boat-tail angle, the others are pure "shaping" terms.
%
% INPUTS
%   t - points in [0,1] (any size, used as a column)
%   n - order of the polynomial (n = 0 gives a column of ones)
% OUTPUT
%   B - numel(t) x (n+1) matrix, column r+1 is the r-th polynomial
%
% See also: classFunction, cstCurve, spanwiseBasis
arguments
    t double {mustBeInRange(t,0,1)}
    n (1,1) double {mustBeInteger,mustBeNonnegative}
end
t = t(:);
r = 0:n;                                               % row: one column per polynomial
K = factorial(n)./(factorial(r).*factorial(n-r));      % binomial coefficients
B = K.*t.^r.*(1-t).^(n-r);
end
