classdef Wing3D_class
% WING3D_CLASS  Wing in the physical domain: sections along the span.
%
%   W = Wing3D_class(eta,b,c,twist,Airfoil)                       % NUIST-like
%   W = Wing3D_class(eta,b,c,twist,Airfoil,Xle=Xle,Zle=Zle)       % LE position
%   W = Wing3D_class(eta,b,c,twist,Airfoil,sweep=sw,dihedral=dh)  % LE angles
%   W = Wing3D_class(...,eta_breaks=0.35,twist_axis=0.25,b0=0,name="wing")
%   W = Wing3D_class.fromSections(Y,XZ)                           % 3D sections
%
% Same arguments as Wing_class of NUIST (eta, b, c, twist, Airfoil), plus
% the leading edge position of every section, so sweep, crank and dihedral
% are free. The class is called Wing3D_class so that it does not shadow the
% NUIST Wing_class when both repositories are on the path.
%
% Axes: x aft, y along the span (root to tip), z up.
% Section k is the airfoil Airfoil{k} (x/c, z/c, untwisted) scaled by
% c(k), rotated nose up by twist(k) about x/c = twist_axis, and placed
% with its untwisted leading edge at (Xle(k), Yle(k), Zle(k)), in the plane
% y = Yle(k) (as WING2CATIA and WiMOT).
%
% Fit and build:
%   CST = Wing2CST(W)      spanwise CST of the wing (least squares)
%   W   = CST2Wing(CST)    wing from the CST parameters
%
% PROPERTIES
%   eta        - N x 1 y/b of the sections in [0,1], ascending (root 0, tip 1)
%   b          - span of the (half) wing along y, m
%   b0         - y of the root section, m (Yle = b0 + eta*b)
%   c          - N x 1 chord, m
%   twist      - N x 1 rad, positive nose up
%   Airfoil    - N x 1 cell of Airfoil_class (x/c, z/c, untwisted)
%   Xle, Zle   - N x 1 untwisted leading edge, m
%   twist_axis - x/c of the twist axis (default 0.25)
%   eta_breaks - planform breaks (kinks) in (0,1): the CST is piecewise
%                smooth along the span with these breaks (see spanwiseBasis)
%   name       - label
% DEPENDENT
%   Yle, S (planform area of this half wing), mac, AR (= 2 b^2/S, full
%   wing), sweepLE and dihedral of each panel (rad)
%
% See also: Wing2CST, CST2Wing, Airfoil_class, spanwiseBasis
    properties
        eta        (:,1) double {mustBeInRange(eta,0,1)}       = zeros(0,1)
        b          (1,1) double {mustBePositive,mustBeFinite}  = 1        % m
        b0         (1,1) double {mustBeFinite}                 = 0        % m
        c          (:,1) double {mustBePositive,mustBeFinite}  = zeros(0,1) % m
        twist      (:,1) double {mustBeFinite}                 = zeros(0,1) % rad
        Airfoil    (:,1) cell                                  = cell(0,1)
        Xle        (:,1) double {mustBeFinite}                 = zeros(0,1) % m
        Zle        (:,1) double {mustBeFinite}                 = zeros(0,1) % m
        twist_axis (1,1) double {mustBeFinite}                 = 0.25     % x/c
        eta_breaks (1,:) double {mustBeInRange(eta_breaks,0,1,'exclusive')} = zeros(1,0)
        name       (1,1) string                                = "wing"
    end
    properties (Dependent)
        Yle        % m
        S          % planform area of this half wing, m^2
        mac        % mean aerodynamic chord, m
        AR         % aspect ratio of the full wing, 2*b^2/S
        sweepLE    % (N-1) x 1 leading edge sweep of each panel, rad
        dihedral   % (N-1) x 1 dihedral of each panel, rad
    end
    % ----------------------------------------------------
    methods
        % constructor
        function obj = Wing3D_class(eta,b,c,twist,Airfoil,options)
            arguments
                eta     (:,1) double = zeros(0,1)
                b       (1,1) double = 1
                c       (:,1) double = zeros(0,1)
                twist   (:,1) double = zeros(0,1)
                Airfoil = {}
                options.Xle        (:,1) double = zeros(0,1)   % untwisted LE x, m
                options.Zle        (:,1) double = zeros(0,1)   % untwisted LE z, m
                options.sweep      (:,1) double = zeros(0,1)   % alternative to Xle: LE sweep of each panel (scalar or N-1), rad
                options.dihedral   (:,1) double = zeros(0,1)   % alternative to Zle: dihedral of each panel (scalar or N-1), rad
                options.b0         (1,1) double = 0
                options.twist_axis (1,1) double = 0.25
                options.eta_breaks (1,:) double = zeros(1,0)
                options.name       (1,1) string = "wing"
            end
            if nargin == 0; return; end
            N              = numel(eta);
            obj.eta        = eta;
            obj.b          = b;
            obj.b0         = options.b0;
            obj.c          = expand(c,N,'c');
            obj.twist      = expand(twist,N,'twist');
            obj.Airfoil    = toAirfoilCell(Airfoil,N);
            obj.twist_axis = options.twist_axis;
            obj.eta_breaks = options.eta_breaks;
            obj.name       = options.name;
            % leading edge: coordinates or panel angles
            dy      = diff(eta)*b;
            obj.Xle = leadingEdge(options.Xle,options.sweep,dy,N,'Xle','sweep');
            obj.Zle = leadingEdge(options.Zle,options.dihedral,dy,N,'Zle','dihedral');
            obj.check();
        end
        % function to check
        function N = check(obj)
            N = numel(obj.eta);
            if N < 2; error('\nWing3D_class: at least 2 sections are needed\n'); end
            for p = ["c","twist","Airfoil","Xle","Zle"]
                if numel(obj.(p)) ~= N
                    error('\nWing3D_class: property %s must have the same size of eta (%d)\n',p,N);
                end
            end
            if any(diff(obj.eta) <= 0); error('\nWing3D_class: eta must be strictly ascending\n'); end
            for i = 1:N
                if ~isa(obj.Airfoil{i},'Airfoil_class')
                    error('\nWing3D_class: Airfoil{%d} must be an Airfoil_class\n',i);
                end
            end
            missing = setdiff(obj.eta_breaks,obj.eta);
            if ~isempty(missing)
                warning('Wing3D_class:breakNotSection', ...
                        'eta_breaks = %s is not a section: the kink is not described exactly',mat2str(missing,4));
            end
        end
        % dependent properties
        function v = get.Yle(obj);      v = obj.b0 + obj.eta*obj.b;                          end
        function v = get.S(obj);        v = trapz(obj.Yle,obj.c);                            end
        function v = get.mac(obj);      v = trapz(obj.Yle,obj.c.^2)/obj.S;                   end
        function v = get.AR(obj);       v = 2*obj.b^2/obj.S;                                 end
        function v = get.sweepLE(obj);  v = atan(diff(obj.Xle)./diff(obj.Yle));             end
        function v = get.dihedral(obj); v = atan(diff(obj.Zle)./diff(obj.Yle));             end
        % same wing with normalized section airfoils (see normalizeSections)
        function obj = normalizeSections(obj)
            [obj.Airfoil,obj.c,obj.twist,obj.Xle,obj.Zle] = ...
                normalizeSections(obj.Airfoil,obj.c,obj.twist,obj.Xle,obj.Zle,obj.twist_axis);
        end
        % 3D coordinates of section k on the chordwise distribution psi
        function [Xu,Zu,Xl,Zl,Y] = section3D(obj,k,psi)
            % upper and lower side from LE to TE (row vectors), y = Yle(k)
            [~,zu,zl] = obj.Airfoil{k}.resample(psi);
            ck        = obj.c(k);
            [xu,zu]   = rotateSection(psi*ck,zu*ck,obj.twist(k),obj.twist_axis*ck);
            [xl,zl]   = rotateSection(psi*ck,zl*ck,obj.twist(k),obj.twist_axis*ck);
            Xu = obj.Xle(k) + xu';  Zu = obj.Zle(k) + zu';
            Xl = obj.Xle(k) + xl';  Zl = obj.Zle(k) + zl';
            Y  = obj.Yle(k)*ones(size(Xu));
        end
        % surface grids (upper and lower side), N sections x npsi points
        function zones = surface(obj,npsi)
            arguments
                obj
                npsi (1,1) double {mustBeInteger,mustBeGreaterThan(npsi,2)} = 81
            end
            psi = cosineSpacing(npsi);
            N   = obj.check();
            [Xu,Zu,Xl,Zl,Y] = deal(zeros(N,npsi));
            for k = 1:N
                [Xu(k,:),Zu(k,:),Xl(k,:),Zl(k,:),Y(k,:)] = obj.section3D(k,psi);
            end
            zones(1) = struct('name',obj.name + " upper",'X',Xu,'Y',Y,'Z',Zu);
            zones(2) = struct('name',obj.name + " lower",'X',Xl,'Y',Y,'Z',Zl);
        end
        % 3D plot: surface and section lines
        function plot(obj,ax,options)
            arguments
                obj
                ax = []
                options.mirror   (1,1) logical = false     % plot the other half wing too
                options.sections (1,1) logical = true      % draw the sections
                options.FaceColor = [0.80 0.82 0.88]
            end
            if isempty(ax); ax = gca; end
            zones = obj.surface();
            plotSurface(zones,ax,mirror=options.mirror,FaceColor=options.FaceColor);
            if options.sections
                for k = 1:numel(obj.eta)
                    lw = 0.5 + 1.5*any(abs(obj.eta(k) - [0 obj.eta_breaks 1]) < 1e-12);   % root, breaks, tip thicker
                    plot3(ax,[flip(zones(1).X(k,:)) zones(2).X(k,2:end)],[flip(zones(1).Y(k,:)) zones(2).Y(k,2:end)], ...
                          [flip(zones(1).Z(k,:)) zones(2).Z(k,2:end)],'k-','LineWidth',lw,'HandleVisibility','off');
                end
            end
        end
        % planform plot (x down, as usually drawn)
        function plotPlanform(obj,ax,varargin)
            if nargin < 2 || isempty(ax); ax = gca; end
            xte = obj.Xle + obj.c;     % untwisted planform
            plot(ax,[obj.Yle; flip(obj.Yle); obj.Yle(1)],[obj.Xle; flip(xte); obj.Xle(1)],varargin{:});
            hold(ax,'on');
            plot(ax,[obj.Yle obj.Yle]',[obj.Xle xte]','k:','HandleVisibility','off');
            set(ax,'YDir','reverse'); axis(ax,'equal');
            xlabel(ax,'y'); ylabel(ax,'x');
        end
    end
    % ----------------------------------------------------
    methods (Static)
        % wing from 3D sections (e.g. cut from a CAD model)
        function obj = fromSections(Y,XZ,options)
            % W = Wing3D_class.fromSections(Y,XZ)
            %   Y  - N x 1 y of the sections (planes y = const), root first
            %   XZ - N x 1 cell, XZ{k} = [x z] of section k (any orientation, m)
            % chord, twist and leading edge of each section are found with
            % Airfoil_class.normalize; Xle, Zle are the untwisted LE about
            % twist_axis, so that the wing rebuilds the given sections.
            arguments
                Y  (:,1) double
                XZ cell
                options.twist_axis (1,1) double = 0.25
                options.eta_breaks (1,:) double = zeros(1,0)
                options.name       (1,1) string = "wing"
            end
            N = numel(Y);
            [c,twist,Xle,Zle] = deal(zeros(N,1));
            AF = cell(N,1);
            xa = options.twist_axis;
            for k = 1:N
                [AF{k},c(k),twist(k),xle,zle] = Airfoil_class(XZ{k}).normalize();
                Xle(k) = xle - xa*c(k)*(1 - cos(twist(k)));     % undo the rotation about the twist axis
                Zle(k) = zle - xa*c(k)*sin(twist(k));
            end
            b   = Y(end) - Y(1);
            obj = Wing3D_class((Y - Y(1))/b,b,c,twist,AF,Xle=Xle,Zle=Zle,b0=Y(1), ...
                               twist_axis=xa,eta_breaks=options.eta_breaks,name=options.name);
        end
    end
end

% FUNCTIONs ---------------------------------------------
function v = expand(v,N,label)
% scalar -> N x 1
if isscalar(v); v = repmat(v,N,1); end
if numel(v) ~= N
    error('\nWing3D_class: %s must be scalar or have the same size of eta (%d)\n',label,N);
end
end

function q = leadingEdge(q,angle,dy,N,qname,aname)
% LE coordinate from its values or from the panel angles (q(1) = 0)
if ~isempty(q) && ~isempty(angle)
    error('\nWing3D_class: give %s or %s, not both\n',qname,aname);
end
if ~isempty(angle)
    if isscalar(angle); angle = repmat(angle,N-1,1); end
    if numel(angle) ~= N-1
        error('\nWing3D_class: %s must be scalar or have one value per panel (%d)\n',aname,N-1);
    end
    q = [0; cumsum(dy.*tan(angle))];
elseif isempty(q)
    q = zeros(N,1);
else
    q = expand(q,N,qname);
end
end
