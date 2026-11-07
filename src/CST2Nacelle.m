function Nac = CST2Nacelle(CST,options)
% CST2NACELLE  Build a nacelle from its CST parameters.
%
%   Nac = CST2Nacelle(CST)
%   Nac = CST2Nacelle(CST,phi=linspace(0,pi,25),npsi=101,name="my nacelle")
%
% Every meridional profile at phi is rebuilt from the Bernstein expansion
% in t = (1 - cos(phi))/2 (see Nacelle2CST):
%   airfoil   Au(t) = Bs(t)*Bu, Al(t) = Bs(t)*Bl   -> CST2Airfoil
%   profile   Rle, Xle, c, twist, dzu, dzl = Bp(t)*[control values]
%
% INPUTS
%   CST  - nacelle CST struct (Nacelle2CST, or written by hand: see
%          examples/build_nacelle.m for the minimum set of fields)
% OPTIONS (name-value)
%   phi  - meridional angles of the output profiles in [0,pi] (default 19
%          equally spaced, or 0 for an axisymmetric CST)
%   npsi - chordwise points per airfoil side (default 81)
%   name - label
% OUTPUT
%   Nac  - Nacelle_class
%
% See also: Nacelle2CST, Nacelle_class, sectionCST, CST2Airfoil
arguments
    CST struct
    options.phi  double = []
    options.npsi (1,1) double {mustBeInteger,mustBeGreaterThan(options.npsi,2)} = 81
    options.name (1,1) string = "CST nacelle"
end
CST = completeNacelleCST(CST);
phi = options.phi(:);
if isempty(phi)
    if CST.ny == 0 && CST.ny_plan == 0
        phi = 0;                          % axisymmetric
    else
        phi = linspace(0,pi,19)';
    end
end
t = (1 - cos(phi))/2;

% profiles
sec = sectionCST(CST,t);
AF  = cell(numel(phi),1);
for k = 1:numel(phi)
    AF{k} = CST2Airfoil(sec(k),options.npsi,name=sprintf('phi = %.1f deg',rad2deg(phi(k))));
end
Bp  = spanwiseBasis(t,CST.eta_breaks,CST.ny_plan);
Nac = Nacelle_class(phi,Bp*CST.Rle,Bp*CST.c,Bp*CST.twist,AF,Xle=Bp*CST.Xle, ...
                    twist_axis=CST.twist_axis,name=options.name);
end

% FUNCTIONs ---------------------------------------------
function CST = completeNacelleCST(CST)
% defaults of the optional fields, so that a CST written by hand needs only
% Bu, Bl, Rle, c
for f = ["Bu","Bl","Rle","c"]
    if ~isfield(CST,f); error('\nCST2Nacelle: the nacelle CST needs the field %s\n',f); end
end
CST.eta_breaks = zeros(1,0);
CST.Rle        = CST.Rle(:);
nPlan          = numel(CST.Rle);
if ~isfield(CST,'ny');      CST.ny      = size(CST.Bu,1) - 1; end
if ~isfield(CST,'ny_plan'); CST.ny_plan = nPlan - 1;          end
if ~isfield(CST,'nx');      CST.nx      = size(CST.Bu,2) - 1; end
defaults = struct('N1',0.5,'N2',1,'twist_axis',0,'Xle',zeros(nPlan,1),'twist',zeros(nPlan,1), ...
                  'dzu',zeros(nPlan,1),'dzl',zeros(nPlan,1));
f = fieldnames(defaults);
for i = 1:numel(f)
    if ~isfield(CST,f{i}) || isempty(CST.(f{i})); CST.(f{i}) = defaults.(f{i}); end
end
for f = ["c","Xle","twist","dzu","dzl"]
    CST.(f) = CST.(f)(:);
    if isscalar(CST.(f)); CST.(f) = repmat(CST.(f),nPlan,1); end   % constant around the nacelle
    if numel(CST.(f)) ~= nPlan
        error('\nCST2Nacelle: CST.%s must have the same number of control values of CST.Rle (%d)\n',f,nPlan);
    end
end
if ~isequal(size(CST.Bl),size(CST.Bu))
    error('\nCST2Nacelle: CST.Bu and CST.Bl must have the same size\n');
end
end
