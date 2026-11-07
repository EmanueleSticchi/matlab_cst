classdef Fuselage_class
% FUSELAGE_CLASS  Fuselage in the physical domain: longitudinal lines and section shape.
%
%   F = Fuselage_class(x,zc,zk,w)                          % elliptic sections
%   F = Fuselage_class(x,zc,zk,w,zm=zm,NU=NU,NL=NL,name="fuselage")
%   F = Fuselage_class.fromSections(x,Ysec,Zsec)           % from cross sections
%
% Kulfan (2008, Sec. VIII-IX): a body is a distribution of cross sections
% along its axis. Every cross section (x = const) is made of an upper and
% a lower lobe joined at the maximum half-breadth (Fig. 20-21):
%   y = w*eta,                         eta in [-1, 1]
%   z = zm + hU*(1 - eta^2)^NU         upper lobe, hU = zc - zm
%   z = zm - hL*(1 - eta^2)^NL         lower lobe, hL = zm - zk
% (1-eta^2)^N is the class function C(t) = t^N(1-t)^N, t = (1+eta)/2,
% scaled to 1 at the centre line: N = 0.5 gives an ellipse, N -> 0 a
% rectangle, N = 1 a parabola, larger N a pointed (ridge) lobe.
%
% Axes: x aft (nose first), y to the side, z up; symmetric about y = 0.
%
% Fit and build:
%   CST = Fuselage2CST(F)       CST parameters (least squares)
%   F   = CST2Fuselage(CST)     fuselage from the CST parameters
%
% PROPERTIES (all N x 1, at the stations x)
%   x      - stations, m, ascending (x(1) = nose, x(end) = tail)
%   zc     - crown line z (top of the section), m
%   zk     - keel line z (bottom of the section), m
%   zm     - z of the maximum half-breadth line, m (default (zc+zk)/2)
%   w      - maximum half-breadth, m
%   NU, NL - class exponents of the upper and lower lobe (default 0.5)
%   name   - label
% DEPENDENT
%   L (length), hU = zc - zm, hL = zm - zk
%
% See also: Fuselage2CST, CST2Fuselage
    properties
        x    (:,1) double {mustBeFinite}                     = zeros(0,1) % m
        zc   (:,1) double {mustBeFinite}                     = zeros(0,1) % m
        zk   (:,1) double {mustBeFinite}                     = zeros(0,1) % m
        zm   (:,1) double {mustBeFinite}                     = zeros(0,1) % m
        w    (:,1) double {mustBeNonnegative,mustBeFinite}   = zeros(0,1) % m
        NU   (:,1) double {mustBePositive,mustBeFinite}      = zeros(0,1)
        NL   (:,1) double {mustBePositive,mustBeFinite}      = zeros(0,1)
        name (1,1) string                                    = "fuselage"
    end
    properties (Dependent)
        L    % length, m
        hU   % upper lobe height zc - zm, m
        hL   % lower lobe height zm - zk, m
    end
    % ----------------------------------------------------
    methods
        % constructor
        function obj = Fuselage_class(x,zc,zk,w,options)
            arguments
                x  (:,1) double = zeros(0,1)
                zc (:,1) double = zeros(0,1)
                zk (:,1) double = zeros(0,1)
                w  (:,1) double = zeros(0,1)
                options.zm   (:,1) double = zeros(0,1)    % default (zc+zk)/2
                options.NU   (:,1) double = 0.5
                options.NL   (:,1) double = 0.5
                options.name (1,1) string = "fuselage"
            end
            if nargin == 0; return; end
            N        = numel(x);
            obj.x    = x;
            obj.zc   = zc;
            obj.zk   = zk;
            obj.w    = w;
            if isempty(options.zm); options.zm = (zc + zk)/2; end
            obj.zm   = options.zm;
            obj.NU   = expand(options.NU,N,'NU');
            obj.NL   = expand(options.NL,N,'NL');
            obj.name = options.name;
            obj.check();
        end
        % function to check
        function N = check(obj)
            N = numel(obj.x);
            if N < 3; error('\nFuselage_class: at least 3 stations are needed\n'); end
            for p = ["zc","zk","zm","w","NU","NL"]
                if numel(obj.(p)) ~= N
                    error('\nFuselage_class: property %s must have the same size of x (%d)\n',p,N);
                end
            end
            if any(diff(obj.x) <= 0); error('\nFuselage_class: x must be strictly ascending\n'); end
            tol = 1e-9*obj.L;
            if any(obj.zc < obj.zm - tol) || any(obj.zm < obj.zk - tol)
                error('\nFuselage_class: it must be zk <= zm <= zc at every station\n');
            end
        end
        % dependent properties
        function v = get.L(obj);  v = obj.x(end) - obj.x(1);  end
        function v = get.hU(obj); v = obj.zc - obj.zm;        end
        function v = get.hL(obj); v = obj.zm - obj.zk;        end
        % cross section k, closed curve (y, z), ntheta points
        function [y,z] = section(obj,k,ntheta)
            % theta from 0 (side, +y) over the crown (pi/2), the other
            % side (pi) and the keel (3/2 pi) back to 0: eta = cos(theta)
            if nargin < 3; ntheta = 73; end
            th  = linspace(0,2*pi,ntheta);
            y   = obj.w(k)*cos(th);
            s2  = sin(th).^2;                                   % 1 - eta^2
            up  = sin(th) >= 0;
            z   = zeros(size(th));
            z(up)  = obj.zm(k) + obj.hU(k)*s2(up).^obj.NU(k);
            z(~up) = obj.zm(k) - obj.hL(k)*s2(~up).^obj.NL(k);
        end
        % surface grid: rows = stations, columns = points around the section
        function zones = surface(obj,ntheta)
            arguments
                obj
                ntheta (1,1) double {mustBeInteger,mustBeGreaterThan(ntheta,4)} = 73
            end
            N = obj.check();
            [Y,Z] = deal(zeros(N,ntheta));
            for k = 1:N
                [Y(k,:),Z(k,:)] = obj.section(k,ntheta);
            end
            zones = struct('name',obj.name,'X',repmat(obj.x,1,ntheta),'Y',Y,'Z',Z);
        end
        % 3D plot
        function plot(obj,ax,options)
            arguments
                obj
                ax = []
                options.sections (1,1) logical = true   % draw the given stations
                options.FaceColor = [0.80 0.82 0.88]
            end
            if isempty(ax); ax = gca; end
            zones = obj.surface();
            plotSurface(zones,ax,FaceColor=options.FaceColor);
            if options.sections
                plot3(ax,zones.X',zones.Y',zones.Z','k-','LineWidth',0.5,'HandleVisibility','off');
            end
        end
        % side view: crown, keel and maximum half-breadth lines; top view: half-breadth
        function plotLines(obj,ax,varargin)
            if nargin < 2 || isempty(ax); ax = gca; end
            hold(ax,'on');
            plot(ax,obj.x,obj.zc,varargin{:},'DisplayName','crown');
            plot(ax,obj.x,obj.zk,varargin{:},'DisplayName','keel');
            plot(ax,obj.x,obj.zm,'--',varargin{:},'DisplayName','max half-breadth line');
            plot(ax,obj.x,obj.w,':',varargin{:},'DisplayName','half-breadth');
            axis(ax,'equal'); xlabel(ax,'x');
        end
    end
    % ----------------------------------------------------
    methods (Static)
        % fuselage from cross sections (e.g. cut from a CAD model)
        function obj = fromSections(x,Ysec,Zsec,options)
            % F = Fuselage_class.fromSections(x,Ysec,Zsec)
            %   x          - N x 1 stations (planes x = const), nose first
            %   Ysec, Zsec - N x 1 cells, points of each section (one side
            %                or both sides, any order). A section with a
            %                single point (nose, tail) is allowed.
            % zc, zk are the highest and lowest points, w the largest |y|
            % and zm its z; NU and NL are fitted (fminbnd) on the points of
            % each lobe. Where a lobe has no height NU, NL = 0.5.
            arguments
                x    (:,1) double
                Ysec cell
                Zsec cell
                options.name (1,1) string = "fuselage"
            end
            N = numel(x);
            [zc,zk,zm,w] = deal(zeros(N,1));
            [NU,NL]      = deal(0.5*ones(N,1));
            for k = 1:N
                y = abs(Ysec{k}(:)); z = Zsec{k}(:);
                zc(k) = max(z);
                zk(k) = min(z);
                [w(k),i] = max(y);
                zm(k) = z(i);
                if w(k) <= 1e-9*max(abs(x)); zm(k) = (zc(k) + zk(k))/2; continue; end   % point section
                NU(k) = fitLobe(y/w(k), z - zm(k),  zc(k) - zm(k));
                NL(k) = fitLobe(y/w(k), zm(k) - z,  zm(k) - zk(k));
            end
            obj = Fuselage_class(x,zc,zk,w,zm=zm,NU=NU,NL=NL,name=options.name);
        end
    end
end

% FUNCTIONs ---------------------------------------------
function v = expand(v,N,label)
% scalar -> N x 1
if isscalar(v); v = repmat(v,N,1); end
if numel(v) ~= N
    error('\nFuselage_class: %s must be scalar or have the same size of x (%d)\n',label,N);
end
end

function N = fitLobe(eta,dz,h)
% exponent N of the lobe dz = h*(1 - eta^2)^N on the points with dz > 0
N   = 0.5;
pts = dz > 0 & eta < 1;
if h <= 0 || nnz(pts) < 2; return; end
f   = @(n) sum((dz(pts) - h*(1 - eta(pts).^2).^n).^2);
N   = fminbnd(f,1e-3,10);
end
