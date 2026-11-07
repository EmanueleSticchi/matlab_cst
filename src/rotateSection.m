function [xr,zr] = rotateSection(x,z,twist,x_axis)
% ROTATESECTION  Rotate a 2D section by its twist (positive nose up).
%
%   [xr,zr] = rotateSection(x,z,twist,x_axis)
%
% x is positive aft and z positive up. A positive twist raises the leading
% edge (nose up) and lowers the trailing edge, as the twist of NUIST, WiMOT
% and WING2CATIA. The rotation is about the point (x_axis, 0), e.g. the
% quarter chord (x_axis = 0.25*c).
%
% INPUTS
%   x, z    - section coordinates (same size)
%   twist   - rad, positive nose up
%   x_axis  - x of the rotation axis (same units as x), default 0
% OUTPUTS
%   xr, zr  - rotated coordinates, same size as x
arguments
    x       double
    z       double
    twist   (1,1) double {mustBeFinite}
    x_axis  (1,1) double {mustBeFinite} = 0
end
dx = x - x_axis;
xr = x_axis + dx*cos(twist) + z*sin(twist);
zr =        - dx*sin(twist) + z*cos(twist);
end
