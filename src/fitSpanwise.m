function [Q,ny,err] = fitSpanwise(eta,q,eta_breaks,ny,options)
% FITSPANWISE  Least squares control values of a spanwise distribution.
%
%   [Q,ny,err] = fitSpanwise(eta,q,eta_breaks,ny)
%   [Q,ny,err] = fitSpanwise(eta,q,eta_breaks,[],maxOrder=3)   % automatic ny
%
% Finds Q such that q(eta) ~ spanwiseBasis(eta,eta_breaks,ny)*Q in the least
% squares sense. q can have several columns (e.g. all the CST coefficients
% of the sections, one column each): every column is fitted on its own.
%
% INPUTS
%   eta        - N x 1 stations in [0,1] where q is known
%   q          - N x m values at the stations
%   eta_breaks - interior breaks in (0,1) ([] for none)
%   ny         - order in each segment; [] = automatic: the largest order
%                the stations can determine in every segment, capped at
%                options.maxOrder (default 3)
% OUTPUTS
%   Q          - (Nseg*ny+1) x m control values
%   ny         - order actually used
%   err        - N x m residual q - B*Q
%
% Each segment needs at least ny+1 stations (both its ends included) for a
% unique fit, otherwise an error explains what to change.
%
% See also: spanwiseBasis, Wing2CST, Nacelle2CST
arguments
    eta        (:,1) double {mustBeInRange(eta,0,1)}
    q          double {mustBeFinite}
    eta_breaks double = zeros(1,0)
    ny         double = []
    options.maxOrder (1,1) double {mustBeInteger,mustBeNonnegative} = 3
end
if size(q,1) ~= numel(eta)
    error('\nfitSpanwise: q must have one row per station (%d rows expected, got %d)\n',numel(eta),size(q,1));
end
% stations in each segment (ends included)
edges   = [0; sort(eta_breaks(:)); 1];
Nseg    = numel(edges) - 1;
perSeg  = zeros(Nseg,1);
for k = 1:Nseg
    perSeg(k) = nnz(eta >= edges(k)-1e-12 & eta <= edges(k+1)+1e-12);
end
if isempty(ny)
    ny = max(0,min(options.maxOrder,min(perSeg)-1));
end
B = spanwiseBasis(eta,eta_breaks,ny);
if rank(B) < size(B,2)
    error(['\nfitSpanwise: a spanwise order ny = %d needs at least %d stations in each segment ' ...
           '(ends included), the segments have [%s] stations.\nLower ny or add sections.\n'], ...
           ny,ny+1,num2str(perSeg'));
end
Q   = B \ q;
err = q - B*Q;
end
