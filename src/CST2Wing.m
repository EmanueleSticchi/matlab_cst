function Wing = CST2Wing(CST,options)
% CST2WING  Build a wing from its spanwise CST parameters.
%
%   Wing = CST2Wing(CST)
%   Wing = CST2Wing(CST,eta=linspace(0,1,41),npsi=101,name="my wing")
%
% Every section at eta is rebuilt from the spanwise Bernstein expansion
% (Kulfan 2008, Eq. 46-48):
%   airfoil   Au(eta) = Bs(eta)*Bu, Al(eta) = Bs(eta)*Bl   -> CST2Airfoil
%   planform  c, twist, Xle, Zle, dzu, dzl = Bp(eta)*[control values]
% with Bs, Bp the piecewise Bernstein bases of order ny and ny_plan
% (spanwiseBasis).
%
% INPUTS
%   CST  - wing CST struct (Wing2CST, or written by hand: see
%          examples/build_wing.m for the minimum set of fields)
% OPTIONS (name-value)
%   eta  - y/b of the output sections (default 21 equally spaced); the
%          breaks are always added, so the kinks are exact
%   npsi - chordwise points per airfoil side (default 81)
%   name - label of the wing
% OUTPUT
%   Wing - Wing3D_class
%
% See also: Wing2CST, Wing3D_class, sectionCST, CST2Airfoil
arguments
    CST struct
    options.eta  double = []
    options.npsi (1,1) double {mustBeInteger,mustBeGreaterThan(options.npsi,2)} = 81
    options.name (1,1) string = "CST wing"
end
CST = completeWingCST(CST);
eta = options.eta(:);
if isempty(eta); eta = linspace(0,1,21)'; end
eta = unique([eta; CST.eta_breaks(:)]);            % the kinks are always sections

% sections
sec = sectionCST(CST,eta);
AF  = cell(numel(eta),1);
for k = 1:numel(eta)
    AF{k} = CST2Airfoil(sec(k),options.npsi,name=sprintf('eta = %.3f',eta(k)));
end
% planform
Bp    = spanwiseBasis(eta,CST.eta_breaks,CST.ny_plan);
Wing  = Wing3D_class(eta,CST.b,Bp*CST.c,Bp*CST.twist,AF,Xle=Bp*CST.Xle,Zle=Bp*CST.Zle, ...
                     b0=CST.b0,twist_axis=CST.twist_axis,eta_breaks=CST.eta_breaks,name=options.name);
end

% FUNCTIONs ---------------------------------------------
function CST = completeWingCST(CST)
% defaults of the optional fields and size checks, so that a CST written by
% hand needs only b, Bu, Bl, c and (optionally) twist, Xle, Zle, eta_breaks
required = ["b","Bu","Bl","c"];
for f = required
    if ~isfield(CST,f); error('\nCST2Wing: the wing CST needs the field %s\n',f); end
end
if ~isfield(CST,'eta_breaks') || isempty(CST.eta_breaks); CST.eta_breaks = zeros(1,0); end
Nseg = numel(CST.eta_breaks) + 1;
% orders from the number of control values
nShape = size(CST.Bu,1);
nPlan  = numel(CST.c);
if ~isfield(CST,'ny');      CST.ny      = orderFromRows(nShape,Nseg,'Bu'); end
if ~isfield(CST,'ny_plan'); CST.ny_plan = orderFromRows(nPlan,Nseg,'c');   end
if ~isfield(CST,'nx');      CST.nx      = size(CST.Bu,2) - 1;              end
CST.c = CST.c(:);
defaults = struct('b0',0,'N1',0.5,'N2',1,'twist_axis',0.25, ...
                  'twist',zeros(nPlan,1),'Xle',zeros(nPlan,1),'Zle',zeros(nPlan,1), ...
                  'dzu',zeros(nPlan,1),'dzl',zeros(nPlan,1));
f = fieldnames(defaults);
for i = 1:numel(f)
    if ~isfield(CST,f{i}) || isempty(CST.(f{i})); CST.(f{i}) = defaults.(f{i}); end
end
for f = ["twist","Xle","Zle","dzu","dzl"]
    CST.(f) = CST.(f)(:);
    if numel(CST.(f)) ~= nPlan
        error('\nCST2Wing: CST.%s must have the same number of control values of CST.c (%d)\n',f,nPlan);
    end
end
if ~isequal(size(CST.Bl),size(CST.Bu))
    error('\nCST2Wing: CST.Bu and CST.Bl must have the same size\n');
end
end

function ny = orderFromRows(nRows,Nseg,label)
% rows = Nseg*ny + 1  ->  ny  (one row = constant along the span)
if nRows == 1; ny = 0; return; end
ny = (nRows - 1)/Nseg;
if ny ~= round(ny) || ny < 1
    error('\nCST2Wing: CST.%s has %d control values, not compatible with %d segment(s): it must be Nseg*ny+1\n',label,nRows,Nseg);
end
end
