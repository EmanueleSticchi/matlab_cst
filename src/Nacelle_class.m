classdef Nacelle_class
% NACELLE_CLASS  Nacelle in the physical domain: cowl profiles around the axis.
%
%   Nac = Nacelle_class(phi,Rle,c,twist,Airfoil)
%   Nac = Nacelle_class(phi,Rle,c,twist,Airfoil,Xle=Xle,twist_axis=0,name="nacelle")
%   Nac = Nacelle_class(0,Rle,c,twist,Airfoil)            % axisymmetric
%
% Kulfan (2008, Sec. X, option 1): the nacelle is a set of airfoil-type
% profiles in meridional planes (crown line, maximum half-breadth, keel
% line, ...) distributed around the nacelle axis. Each profile is handled
% exactly as a wing section, with the radius in place of z:
%   upper side of the airfoil = outer cowl (radius increases)
%   lower side of the airfoil = inner duct / inlet (radius decreases)
%
% Axes: x along the nacelle axis (aft), the axis is the line y = z = 0.
% The meridional plane at angle phi contains the axis and the direction
% (y,z) = (sin(phi), cos(phi)): phi = 0 crown (top), pi/2 side (+y),
% pi keel (bottom). The nacelle is symmetric about the plane y = 0, so the
% profiles are given for phi in [0, pi] and mirrored. A single profile
% gives an axisymmetric nacelle.
%
% Profile k: Airfoil{k} (x/c, r/c, untwisted) scaled by c(k), rotated by
% twist(k) about x/c = twist_axis (positive = leading edge outward, i.e.
% "nose up" with the radius as the up direction) and placed with its
% untwisted leading edge (highlight) at x = Xle(k), r = Rle(k).
%
% Fit and build:
%   CST = Nacelle2CST(Nac)      CST parameters (least squares)
%   Nac = CST2Nacelle(CST)      nacelle from the CST parameters
%
% PROPERTIES
%   phi        - N x 1 meridional angle in [0,pi], ascending, rad
%   Rle        - N x 1 highlight radius, m
%   Xle        - N x 1 highlight x, m (default 0)
%   c          - N x 1 profile length (chord), m
%   twist      - N x 1 profile rotation, rad (positive = LE outward)
%   Airfoil    - N x 1 cell of Airfoil_class
%   twist_axis - x/c of the rotation point (default 0: the highlight)
%   name       - label
%
% See also: Nacelle2CST, CST2Nacelle, Wing3D_class
    properties
        phi        (:,1) double {mustBeInRange(phi,0,3.141592653589794)} = zeros(0,1) % rad
        Rle        (:,1) double {mustBePositive,mustBeFinite} = zeros(0,1) % m
        Xle        (:,1) double {mustBeFinite}                = zeros(0,1) % m
        c          (:,1) double {mustBePositive,mustBeFinite} = zeros(0,1) % m
        twist      (:,1) double {mustBeFinite}                = zeros(0,1) % rad
        Airfoil    (:,1) cell                                 = cell(0,1)
        twist_axis (1,1) double {mustBeFinite}                = 0          % x/c
        name       (1,1) string                               = "nacelle"
    end
    % ----------------------------------------------------
    methods
        % constructor
        function obj = Nacelle_class(phi,Rle,c,twist,Airfoil,options)
            arguments
                phi     (:,1) double = zeros(0,1)
                Rle     (:,1) double = zeros(0,1)
                c       (:,1) double = zeros(0,1)
                twist   (:,1) double = zeros(0,1)
                Airfoil = {}
                options.Xle        (:,1) double = 0
                options.twist_axis (1,1) double = 0
                options.name       (1,1) string = "nacelle"
            end
            if nargin == 0; return; end
            N              = numel(phi);
            obj.phi        = phi;
            obj.Rle        = expand(Rle,N,'Rle');
            obj.Xle        = expand(options.Xle,N,'Xle');
            obj.c          = expand(c,N,'c');
            obj.twist      = expand(twist,N,'twist');
            obj.Airfoil    = toAirfoilCell(Airfoil,N);
            obj.twist_axis = options.twist_axis;
            obj.name       = options.name;
            obj.check();
        end
        % function to check
        function N = check(obj)
            N = numel(obj.phi);
            if N < 1; error('\nNacelle_class: at least one profile is needed\n'); end
            for p = ["Rle","Xle","c","twist","Airfoil"]
                if numel(obj.(p)) ~= N
                    error('\nNacelle_class: property %s must have the same size of phi (%d)\n',p,N);
                end
            end
            if any(diff(obj.phi) <= 0); error('\nNacelle_class: phi must be strictly ascending\n'); end
            for i = 1:N
                if ~isa(obj.Airfoil{i},'Airfoil_class')
                    error('\nNacelle_class: Airfoil{%d} must be an Airfoil_class\n',i);
                end
            end
        end
        % true if the nacelle is axisymmetric (one profile)
        function tf = isAxisymmetric(obj)
            tf = numel(obj.phi) == 1;
        end
        % same nacelle with normalized profile airfoils (see normalizeSections)
        function obj = normalizeSections(obj)
            [obj.Airfoil,obj.c,obj.twist,obj.Xle,obj.Rle] = ...
                normalizeSections(obj.Airfoil,obj.c,obj.twist,obj.Xle,obj.Rle,obj.twist_axis);
        end
        % meridional profile k (x, r) on the chordwise distribution psi
        function [xo,ro,xi,ri] = profile(obj,k,psi)
            % outer cowl (xo, ro) and inner duct (xi, ri), from LE to TE
            [~,zu,zl] = obj.Airfoil{k}.resample(psi);
            ck        = obj.c(k);
            [xo,ro]   = rotateSection(psi*ck,zu*ck,obj.twist(k),obj.twist_axis*ck);
            [xi,ri]   = rotateSection(psi*ck,zl*ck,obj.twist(k),obj.twist_axis*ck);
            xo = obj.Xle(k) + xo;  ro = obj.Rle(k) + ro;
            xi = obj.Xle(k) + xi;  ri = obj.Rle(k) + ri;
        end
        % surface grids (outer cowl and inner duct), full 360 deg
        function zones = surface(obj,npsi,nphi)
            % rows: meridional planes from phi = 0 around to 2*pi (closed),
            % columns: chordwise points. The given profiles are mirrored
            % to phi in (pi, 2*pi); an axisymmetric nacelle is revolved
            % with nphi planes.
            arguments
                obj
                npsi (1,1) double {mustBeInteger,mustBeGreaterThan(npsi,2)} = 81
                nphi (1,1) double {mustBeInteger,mustBeGreaterThan(nphi,4)} = 73
            end
            psi = cosineSpacing(npsi);
            N   = obj.check();
            [Xo,Ro,Xi,Ri] = deal(zeros(N,npsi));
            for k = 1:N
                [Xo(k,:),Ro(k,:),Xi(k,:),Ri(k,:)] = obj.profile(k,psi);
            end
            if obj.isAxisymmetric()
                ang = linspace(0,2*pi,nphi)';
                [Xo,Ro,Xi,Ri] = deal(repmat(Xo,nphi,1),repmat(Ro,nphi,1),repmat(Xi,nphi,1),repmat(Ri,nphi,1));
            else
                % mirror: phi -> 2*pi - phi (without repeating 0 and pi)
                mir = flip(1:N);
                mir = mir(obj.phi(mir) > 0 & obj.phi(mir) < pi);
                ang = [obj.phi; 2*pi - obj.phi(mir)];
                Xo  = [Xo; Xo(mir,:)];  Ro = [Ro; Ro(mir,:)];
                Xi  = [Xi; Xi(mir,:)];  Ri = [Ri; Ri(mir,:)];
                if obj.phi(1) == 0 && obj.phi(end) == pi    % close the surface at the crown
                    ang = [ang; 2*pi];
                    Xo  = [Xo; Xo(1,:)];  Ro = [Ro; Ro(1,:)];
                    Xi  = [Xi; Xi(1,:)];  Ri = [Ri; Ri(1,:)];
                end
            end
            zones(1) = struct('name',obj.name + " outer",'X',Xo,'Y',Ro.*sin(ang),'Z',Ro.*cos(ang));
            zones(2) = struct('name',obj.name + " inner",'X',Xi,'Y',Ri.*sin(ang),'Z',Ri.*cos(ang));
        end
        % 3D plot
        function plot(obj,ax,options)
            arguments
                obj
                ax = []
                options.profiles (1,1) logical = true   % draw the given profiles
                options.FaceColor = [0.80 0.82 0.88]
            end
            if isempty(ax); ax = gca; end
            zones = obj.surface();
            plotSurface(zones,ax,FaceColor=options.FaceColor);
            if options.profiles
                psi = cosineSpacing(81);
                for k = 1:numel(obj.phi)
                    [xo,ro,xi,ri] = obj.profile(k,psi);
                    x = [flip(xo); xi(2:end)]; r = [flip(ro); ri(2:end)];
                    plot3(ax,x,r*sin(obj.phi(k)),r*cos(obj.phi(k)),'k-','LineWidth',1.5,'HandleVisibility','off');
                end
            end
        end
        % plot of the profiles in the (x, r) plane
        function plotProfiles(obj,ax,varargin)
            if nargin < 2 || isempty(ax); ax = gca; end
            psi = cosineSpacing(101);
            hold(ax,'on');
            for k = 1:numel(obj.phi)
                [xo,ro,xi,ri] = obj.profile(k,psi);
                plot(ax,[flip(xo); xi(2:end)],[flip(ro); ri(2:end)],varargin{:}, ...
                     'DisplayName',sprintf('\\phi = %.0f deg',rad2deg(obj.phi(k))));
            end
            axis(ax,'equal'); xlabel(ax,'x'); ylabel(ax,'r');
        end
    end
end

% FUNCTIONs ---------------------------------------------
function v = expand(v,N,label)
% scalar -> N x 1
if isscalar(v); v = repmat(v,N,1); end
if numel(v) ~= N
    error('\nNacelle_class: %s must be scalar or have the same size of phi (%d)\n',label,N);
end
end
