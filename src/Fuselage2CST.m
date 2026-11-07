function CST = Fuselage2CST(F,options)
% FUSELAGE2CST  Fit the CST parameters of a fuselage (least squares).
%
%   CST = Fuselage2CST(F)
%   CST = Fuselage2CST(F,n=8,nz=6,nN=3,N1=0.5,N2=1)
%
% With psi = (x - x0)/L along the body (Kulfan 2008, Sec. IX, Eq. 29-34):
%   lobe heights and half-breadth, closed at the nose, are CST curves
%     hU/L = C(psi)*Bern_n(psi)*Ahu + psi*dzhu
%     hL/L = C(psi)*Bern_n(psi)*Ahl + psi*dzhl
%     w/L  = C(psi)*Bern_n(psi)*Aw  + psi*dzw
%   with C(psi) = psi^N1 (1-psi)^N2 (N1 = 0.5: round nose; dz = height at
%   the tail, 0 for a closed tail)
%   the line of maximum half-breadth (camber-like line of the body) and the
%   lobe exponents are plain Bernstein polynomials (no class function):
%     zm/L = Bern_nz(psi)*Bzm
%     NU   = Bern_nN(psi)*BNU,    NL = Bern_nN(psi)*BNL
% All the unknowns enter linearly: each one is a backslash.
%
% INPUTS
%   F - Fuselage_class
% OPTIONS (name-value)
%   n      - Bernstein order of hU, hL, w (default 8)
%   nz     - Bernstein order of zm (default 6)
%   nN     - Bernstein order of NU, NL (default 3)
%   N1, N2 - class function exponents of hU, hL, w (default 0.5, 1)
% OUTPUT
%   CST - struct
%     type                 - "fuselage"
%     x0, L                - nose x and length, m
%     n, nz, nN, N1, N2
%     Ahu, Ahl, Aw         - (n+1) x 1, dzhu, dzhl, dzw (all / L)
%     Bzm                  - (nz+1) x 1 (/ L)
%     BNU, BNL             - (nN+1) x 1
%     err                  - struct, CST - input at the stations (m for
%                            hU, hL, w, zm)
%
% See also: CST2Fuselage, Fuselage_class, fitCstCurve
arguments
    F {mustBeA(F,'Fuselage_class')}
    options.n  (1,1) double {mustBeInteger,mustBeNonnegative} = 8
    options.nz (1,1) double {mustBeInteger,mustBeNonnegative} = 6
    options.nN (1,1) double {mustBeInteger,mustBeNonnegative} = 3
    options.N1 (1,1) double {mustBeNonnegative} = 0.5
    options.N2 (1,1) double {mustBeNonnegative} = 1
end
F.check();
x0  = F.x(1);
L   = F.L;
psi = (F.x - x0)/L;
n   = options.n; N1 = options.N1; N2 = options.N2;

[Ahu,dzhu,ehu] = fitCstCurve(psi,F.hU/L,n,N1,N2);
[Ahl,dzhl,ehl] = fitCstCurve(psi,F.hL/L,n,N1,N2);
[Aw, dzw, ew ] = fitCstCurve(psi,F.w/L, n,N1,N2);
Bz  = bernsteinBasis(psi,options.nz);
BN  = bernsteinBasis(psi,options.nN);
Bzm = Bz \ (F.zm/L);
BNU = BN \ F.NU;
BNL = BN \ F.NL;

CST.type = "fuselage";
CST.x0   = x0;
CST.L    = L;
CST.n    = n;
CST.nz   = options.nz;
CST.nN   = options.nN;
CST.N1   = N1;
CST.N2   = N2;
CST.Ahu  = Ahu;  CST.dzhu = dzhu;
CST.Ahl  = Ahl;  CST.dzhl = dzhl;
CST.Aw   = Aw;   CST.dzw  = dzw;
CST.Bzm  = Bzm;
CST.BNU  = BNU;
CST.BNL  = BNL;
CST.err  = struct('hU',-ehu*L,'hL',-ehl*L,'w',-ew*L,'zm',Bz*Bzm*L - F.zm, ...
                  'NU',BN*BNU - F.NU,'NL',BN*BNL - F.NL);
end
