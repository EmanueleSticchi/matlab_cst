function F = CST2Fuselage(CST,options)
% CST2FUSELAGE  Build a fuselage from its CST parameters.
%
%   F = CST2Fuselage(CST)
%   F = CST2Fuselage(CST,npsi=81,name="my fuselage")
%
%   psi = (x - x0)/L
%   hU  = L*(C(psi)*Bern(psi)*Ahu + psi*dzhu)        upper lobe height
%   hL  = L*(C(psi)*Bern(psi)*Ahl + psi*dzhl)        lower lobe height
%   w   = L*(C(psi)*Bern(psi)*Aw  + psi*dzw)         half-breadth
%   zm  = L*Bern(psi)*Bzm                            max half-breadth line
%   NU  = Bern(psi)*BNU,  NL = Bern(psi)*BNL         lobe exponents
%   zc  = zm + hU,  zk = zm - hL
% (see Fuselage2CST and Fuselage_class)
%
% INPUTS
%   CST  - fuselage CST struct (Fuselage2CST, or written by hand: see
%          examples/build_fuselage.m)
% OPTIONS (name-value)
%   npsi - number of stations, cosine spaced (default 61)
%   name - label
% OUTPUT
%   F    - Fuselage_class
%
% Negative heights or exponents (possible with coefficients set by hand)
% are clipped to 0 and 0.001 with a warning.
%
% See also: Fuselage2CST, Fuselage_class, cstCurve
arguments
    CST struct
    options.npsi (1,1) double {mustBeInteger,mustBeGreaterThan(options.npsi,2)} = 61
    options.name (1,1) string = "CST fuselage"
end
CST = completeFuselageCST(CST);
L   = CST.L;
psi = cosineSpacing(options.npsi);
hU  = L*cstCurve(psi,CST.Ahu,CST.dzhu,CST.N1,CST.N2);
hL  = L*cstCurve(psi,CST.Ahl,CST.dzhl,CST.N1,CST.N2);
w   = L*cstCurve(psi,CST.Aw, CST.dzw, CST.N1,CST.N2);
zm  = L*bernsteinBasis(psi,numel(CST.Bzm)-1)*CST.Bzm;
NU  = bernsteinBasis(psi,numel(CST.BNU)-1)*CST.BNU;
NL  = bernsteinBasis(psi,numel(CST.BNL)-1)*CST.BNL;
if any([hU; hL; w] < 0) || any([NU; NL] <= 0)
    warning('CST2Fuselage:clipped','negative heights, half-breadths or exponents clipped');
    hU = max(hU,0); hL = max(hL,0); w = max(w,0);
    NU = max(NU,1e-3); NL = max(NL,1e-3);
end
F = Fuselage_class(CST.x0 + psi*L,zm + hU,zm - hL,w,zm=zm,NU=NU,NL=NL,name=options.name);
end

% FUNCTIONs ---------------------------------------------
function CST = completeFuselageCST(CST)
% defaults of the optional fields, so that a CST written by hand needs only
% L, Ahu, Ahl, Aw
for f = ["L","Ahu","Ahl","Aw"]
    if ~isfield(CST,f); error('\nCST2Fuselage: the fuselage CST needs the field %s\n',f); end
end
defaults = struct('x0',0,'N1',0.5,'N2',1,'dzhu',0,'dzhl',0,'dzw',0,'Bzm',0,'BNU',0.5,'BNL',0.5);
f = fieldnames(defaults);
for i = 1:numel(f)
    if ~isfield(CST,f{i}) || isempty(CST.(f{i})); CST.(f{i}) = defaults.(f{i}); end
end
for f = ["Ahu","Ahl","Aw","Bzm","BNU","BNL"]
    CST.(f) = CST.(f)(:);
end
end
