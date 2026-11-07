function AF = CST2Airfoil(CST,psi,options)
% CST2AIRFOIL  Build an airfoil from its CST coefficients.
%
%   AF = CST2Airfoil(CST)            % 101 cosine spaced points per side
%   AF = CST2Airfoil(CST,201)        % N points per side
%   AF = CST2Airfoil(CST,psi)        % given x/c distribution (both sides)
%   AF = CST2Airfoil(CST,psi,name="my foil")
%
%   z_u/c = C(psi) * sum_r Au_r S_r(psi) + psi*dzu      (Kulfan 2008, Eq. 17)
%   z_l/c = C(psi) * sum_r Al_r S_r(psi) + psi*dzl      (Eq. 19)
%
% INPUTS
%   CST - struct with fields
%           Au, Al   - upper and lower Bernstein coefficients (same number)
%           dzu, dzl - (optional, default 0) z/c of the two sides at the TE
%           N1, N2   - (optional, default 0.5 and 1) class function exponents
%         e.g. CST = struct('Au',[0.17 0.16 0.15 0.14],'Al',-[0.17 0.12 0.10 0.05])
%   psi - number of points per side (scalar) or x/c distribution in [0,1]
% OUTPUT
%   AF  - Airfoil_class, Selig order (TE -> upper -> LE -> lower -> TE),
%         2*numel(psi)-1 points
%
% See also: Airfoil2CST, cstCurve, Airfoil_class
arguments
    CST struct
    psi double = 101
    options.name (1,1) string = "CST airfoil"
end
CST = completeAirfoilCST(CST);
if isscalar(psi); psi = cosineSpacing(psi); end
psi = psi(:);
zu  = cstCurve(psi,CST.Au,CST.dzu,CST.N1,CST.N2);
zl  = cstCurve(psi,CST.Al,CST.dzl,CST.N1,CST.N2);
AF  = Airfoil_class([flip(psi); psi(2:end)],[flip(zu); zl(2:end)],name=options.name);
end

% FUNCTIONs ---------------------------------------------
function CST = completeAirfoilCST(CST)
% default values of the optional fields and size checks
if ~isfield(CST,'Au') || ~isfield(CST,'Al')
    error('\nCST2Airfoil: CST must have the fields Au and Al\n');
end
CST.Au = CST.Au(:);
CST.Al = CST.Al(:);
if numel(CST.Au) ~= numel(CST.Al)
    error('\nCST2Airfoil: Au and Al must have the same number of coefficients (%d vs %d)\n',numel(CST.Au),numel(CST.Al));
end
defaults = struct('dzu',0,'dzl',0,'N1',0.5,'N2',1);
f = fieldnames(defaults);
for i = 1:numel(f)
    if ~isfield(CST,f{i}) || isempty(CST.(f{i})); CST.(f{i}) = defaults.(f{i}); end
end
end
