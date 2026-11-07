function sec = sectionCST(CST,t)
% SECTIONCST  Airfoil CST of the sections of a wing (or nacelle) CST.
%
%   sec = sectionCST(CST,t)
%
% Evaluates the spanwise Bernstein expansion of a wing / nacelle CST
% (Kulfan 2008, Eq. 43) at the spanwise coordinate t and returns the 2D
% CST of each section, ready for CST2Airfoil:
%   Au(t) = spanwiseBasis(t,eta_breaks,ny)*Bu        (same for Al)
%   dzu(t), dzl(t) with the planform order ny_plan
%
% INPUTS
%   CST - wing CST (Wing2CST, t = eta = y/b) or nacelle CST (Nacelle2CST,
%         t = (1-cos(phi))/2)
%   t   - spanwise points in [0,1]
% OUTPUT
%   sec - numel(t) x 1 struct array with fields Au, Al, dzu, dzl, n, N1, N2
%
% See also: CST2Wing, CST2Nacelle, CST2Airfoil
arguments
    CST struct
    t   double {mustBeInRange(t,0,1)}
end
t   = t(:);
Bs  = spanwiseBasis(t,CST.eta_breaks,CST.ny);        % shape coefficients
Bp  = spanwiseBasis(t,CST.eta_breaks,CST.ny_plan);   % trailing edge thickness (planform order)
Au  = Bs*CST.Bu;
Al  = Bs*CST.Bl;
dzu = Bp*CST.dzu;
dzl = Bp*CST.dzl;
sec = struct('Au',[],'Al',[],'dzu',[],'dzl',[],'n',[],'N1',[],'N2',[]);
sec = repmat(sec,numel(t),1);
for k = 1:numel(t)
    sec(k).Au  = Au(k,:)';
    sec(k).Al  = Al(k,:)';
    sec(k).dzu = dzu(k);
    sec(k).dzl = dzl(k);
    sec(k).n   = CST.nx;
    sec(k).N1  = CST.N1;
    sec(k).N2  = CST.N2;
end
end
