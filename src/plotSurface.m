function h = plotSurface(zones,ax,options)
% PLOTSURFACE  Plot the surface zones of a wing, nacelle or fuselage.
%
%   h = plotSurface(zones)
%   h = plotSurface(zones,ax,FaceColor=[0.8 0.8 0.85],mirror=true)
%
% INPUTS
%   zones - struct array with fields name, X, Y, Z (structured grids, same
%           size), as returned by the surface method of Wing3D_class,
%           Nacelle_class and Fuselage_class
%   ax    - axes (default gca)
% OPTIONS (name-value)
%   FaceColor - surface color (default light grey-blue)
%   EdgeColor - grid line color (default 'none')
%   mirror    - also plot the image about the y = 0 plane (default false)
%   light     - add a camera light (default true)
% OUTPUT
%   h     - surface handles
%
% See also: writeTecplot
arguments
    zones struct
    ax = []
    options.FaceColor = [0.80 0.82 0.88]
    options.EdgeColor = 'none'
    options.mirror (1,1) logical = false
    options.light  (1,1) logical = true
end
if isempty(ax); ax = gca; end
hold(ax,'on');
h = gobjects(0);
for k = 1:numel(zones)
    h(end+1) = surf(ax,zones(k).X,zones(k).Y,zones(k).Z,'FaceColor',options.FaceColor, ...
                    'EdgeColor',options.EdgeColor,'FaceLighting','gouraud', ...
                    'AmbientStrength',0.5,'DisplayName',zones(k).name); %#ok<AGROW>
    if options.mirror
        h(end+1) = surf(ax,zones(k).X,-zones(k).Y,zones(k).Z,'FaceColor',options.FaceColor, ...
                        'EdgeColor',options.EdgeColor,'FaceLighting','gouraud', ...
                        'AmbientStrength',0.5,'HandleVisibility','off'); %#ok<AGROW>
    end
end
axis(ax,'equal','tight');
view(ax,3);
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
if options.light && isempty(findobj(ax,'Type','light'))
    camlight(ax,'headlight');
end
end
