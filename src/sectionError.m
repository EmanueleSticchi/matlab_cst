function [rmsErr,maxErr] = sectionError(AFs,sec)
% SECTIONERROR  z/c error of CST sections against the original airfoils.
%
%   [rmsErr,maxErr] = sectionError(AFs,sec)
%
% INPUTS
%   AFs - N x 1 cell of Airfoil_class (x/c, z/c, the reference)
%   sec - N x 1 struct array of airfoil CST (sectionCST)
% OUTPUTS
%   rmsErr, maxErr - N x 1, rms and max |z/c| error at the points of AFs{k}
%
% See also: Wing2CST, Nacelle2CST
N      = numel(AFs);
rmsErr = zeros(N,1); maxErr = rmsErr;
for k = 1:N
    [xu,zu,xl,zl] = AFs{k}.split();
    xu  = min(max(xu,0),1);  xl = min(max(xl,0),1);
    res = [zu - cstCurve(xu,sec(k).Au,sec(k).dzu,sec(k).N1,sec(k).N2); ...
           zl - cstCurve(xl,sec(k).Al,sec(k).dzl,sec(k).N1,sec(k).N2)];
    rmsErr(k) = sqrt(mean(res.^2));
    maxErr(k) = max(abs(res));
end
end
