function [B,eta_ctrl] = spanwiseBasis(eta,eta_breaks,ny)
% SPANWISEBASIS  Piecewise Bernstein basis along the span, with breaks.
%
%   [B,eta_ctrl] = spanwiseBasis(eta,eta_breaks,ny)
%
%   q(eta) = B(eta) * Q
%
% Kulfan (2008, Sec. XII) writes every CST coefficient of the wing sections
% as a Bernstein polynomial in the spanwise coordinate eta:
%   Au_i(eta) = sum_j Bu_ij Sy_j(eta)
% Planform breaks (kinks of the leading or trailing edge, Sec. XIII) make
% the spanwise distributions only piecewise smooth. Here the span [0,1] is
% cut by eta_breaks into Nseg segments and in each segment q is a Bernstein
% polynomial of order ny in the local coordinate
%   t = (eta - eta_a)/(eta_b - eta_a),   eta_a, eta_b = segment ends.
% A Bernstein polynomial is equal to its first control value at t = 0 and
% to its last one at t = 1, so letting two consecutive segments SHARE the
% control value at the break makes q continuous there, while its slope can
% jump (the kink). The control values are stacked as
%   Q = [segment 1: ny+1 values; segment 2: ny more values; ...]
% for a total of Nseg*ny + 1 rows.
%
% Examples (ny = 1, linear in each segment):
%   eta_breaks = []      Q = [q_root; q_tip]              straight panel
%   eta_breaks = 0.35    Q = [q_root; q_break; q_tip]     cranked wing
%
% INPUTS
%   eta        - spanwise points in [0,1]
%   eta_breaks - interior breaks in (0,1), ascending ([] for no break)
%   ny         - order of the Bernstein polynomial in each segment
%                (ny = 0: q constant along the whole span)
% OUTPUTS
%   B          - numel(eta) x (Nseg*ny+1) matrix
%   eta_ctrl   - (Nseg*ny+1) x 1, eta associated with each control value
%                (segment ends and equally spaced points in between). The
%                control value is equal to q only at the segment ends.
%
% See also: bernsteinBasis, fitSpanwise, sectionCST
arguments
    eta        double {mustBeInRange(eta,0,1)}
    eta_breaks double {mustBeInRange(eta_breaks,0,1,'exclusive')} = zeros(1,0)
    ny         (1,1) double {mustBeInteger,mustBeNonnegative} = 1
end
eta   = eta(:);
edges = [0; sort(eta_breaks(:)); 1];
Nseg  = numel(edges) - 1;
if ny == 0                          % constant: one control value, breaks irrelevant
    B        = ones(numel(eta),1);
    eta_ctrl = 0;
    return
end
Nctrl    = Nseg*ny + 1;
B        = zeros(numel(eta),Nctrl);
eta_ctrl = zeros(Nctrl,1);
for k = 1:Nseg
    in = eta >= edges(k) & eta <= edges(k+1);
    if k > 1; in = in & eta > edges(k); end            % a point on a break belongs to the inner segment (same value anyway)
    t        = (eta(in) - edges(k))/(edges(k+1) - edges(k));
    cols     = (k-1)*ny + (1:ny+1);                    % the first column is shared with the previous segment
    B(in,cols)     = bernsteinBasis(t,ny);
    eta_ctrl(cols) = edges(k) + (0:ny)'/ny*(edges(k+1) - edges(k));
end
end
