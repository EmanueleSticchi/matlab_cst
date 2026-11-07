function [AFn,c,twist,Xle,Zle] = normalizeSections(AFs,c,twist,Xle,Zle,twist_axis)
% NORMALIZESECTIONS  Normalized section airfoils, same physical sections.
%
%   [AFn,c,twist,Xle,Zle] = normalizeSections(AFs,c,twist,Xle,Zle,twist_axis)
%
% The CST describes an airfoil with its LE at (0,0) and its TE midpoint at
% (1,0). A section airfoil given in x/c, z/c may not be exactly like that
% (e.g. a cambered NACA airfoil, whose geometric LE is not at x = 0, or a
% file with a slightly rotated chord line). Here every airfoil is
% normalized (Airfoil_class.normalize) and its residual chord, rotation
% and LE shift are moved into c, twist, Xle, Zle, so that the 3D section
%   [X;Z] = [Xle;Zle] + rotateSection(c*[x;z], twist, twist_axis*c)
% is exactly the same before and after.
%
% INPUTS
%   AFs        - N x 1 cell of Airfoil_class (x/c, z/c)
%   c, twist   - N x 1 chord (m) and twist (rad, nose up)
%   Xle, Zle   - N x 1 untwisted LE (m), z = radius for a nacelle
%   twist_axis - x/c of the twist axis
% OUTPUTS
%   same quantities, with normalized airfoils
%
% See also: Airfoil_class.normalize, Wing2CST, Nacelle2CST
N   = numel(AFs);
AFn = cell(N,1);
xa  = twist_axis;
for k = 1:N
    [AFn{k},cn,twn,xle,zle] = AFs{k}.normalize();
    % physical position of the normalized LE
    [px,pz] = rotateSection(xle*c(k),zle*c(k),twist(k),xa*c(k));
    px = Xle(k) + px;
    pz = Zle(k) + pz;
    % new chord and twist, untwisted LE about the twist axis
    c(k)     = c(k)*cn;
    twist(k) = twist(k) + twn;
    Xle(k)   = px - xa*c(k)*(1 - cos(twist(k)));
    Zle(k)   = pz - xa*c(k)*sin(twist(k));
end
end
