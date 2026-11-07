function [Au,Al,dzu,dzl] = fitSections(AFs,n,N1,N2,LEcontinuity)
% FITSECTIONS  Airfoil CST of every section of a wing or nacelle.
%
%   [Au,Al,dzu,dzl] = fitSections(AFs,n,N1,N2,LEcontinuity)
%
% First step of Wing2CST / Nacelle2CST: each section airfoil (already in
% x/c, z/c, untwisted) is fitted on its own with Airfoil2CST, without
% normalizing it again (its chord and twist belong to the planform).
%
% INPUTS
%   AFs          - N x 1 cell of Airfoil_class
%   n            - chordwise Bernstein order
%   N1, N2       - class function exponents
%   LEcontinuity - impose Al(1) = -Au(1)
% OUTPUTS
%   Au, Al       - N x (n+1), row k = coefficients of section k
%   dzu, dzl     - N x 1 trailing edge z/c
%
% See also: Airfoil2CST, Wing2CST, Nacelle2CST
N   = numel(AFs);
Au  = zeros(N,n+1); Al = Au;
dzu = zeros(N,1);   dzl = dzu;
for k = 1:N
    sec      = Airfoil2CST(AFs{k},n=n,N1=N1,N2=N2,LEcontinuity=LEcontinuity,normalize=false);
    Au(k,:)  = sec.Au';
    Al(k,:)  = sec.Al';
    dzu(k)   = sec.dzu;
    dzl(k)   = sec.dzl;
end
end
