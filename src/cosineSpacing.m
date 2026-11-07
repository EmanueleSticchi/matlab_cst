function psi = cosineSpacing(N)
% COSINESPACING  N points in [0,1] clustered at both ends.
%
%   psi = cosineSpacing(N)
%
%   psi = (1 - cos(pi*(0:N-1)/(N-1)))/2      (column)
%
% Standard chordwise distribution of an airfoil: dense at the leading edge,
% where the curvature is high, and at the trailing edge.
arguments
    N (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(N,2)}
end
psi = (1 - cos(pi*(0:N-1)'/(N-1)))/2;
psi([1 end]) = [0 1];   % exact end points
end
