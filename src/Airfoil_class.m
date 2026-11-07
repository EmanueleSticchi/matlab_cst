classdef Airfoil_class
% AIRFOIL_CLASS  Airfoil in the physical domain: a list of (x,z) points.
%
%   AF = Airfoil_class(x,z)                 % coordinates (column vectors)
%   AF = Airfoil_class([x z])               % N x 2 array
%   AF = Airfoil_class('naca0012.dat')      % Selig or Lednicer file
%   AF = Airfoil_class.naca4("2412")        % NACA 4 digit airfoil
%
% The points are stored in Selig order: from the trailing edge along the
% upper side to the leading edge and back along the lower side to the
% trailing edge. The constructor reorders them if they come the other way.
% x is positive aft, z positive up. For a wing section x, z are x/c, z/c
% of the untwisted section (normalize does that for raw coordinates).
%
% Fit and build:
%   CST = Airfoil2CST(AF)        CST coefficients of AF (least squares)
%   AF  = CST2Airfoil(CST)       airfoil from CST coefficients
%
% METHODS
%   check       - sizes and finite coordinates
%   iLE         - index of the data point farthest from the TE
%   leadingEdge - LE on a spline of the contour (need not be a data point)
%   split       - upper and lower sides, each from LE to TE
%   normalize   - LE at the origin, TE at (1,0): returns chord, twist, LE
%   resample    - same airfoil on a new x/c distribution (both sides)
%   thickness   - maximum thickness t/c and its position
%   plot, write - plot and Selig file
%
% See also: Airfoil2CST, CST2Airfoil, Wing3D_class
    properties
        x    (:,1) double {mustBeFinite} = zeros(0,1)   % x coordinates, Selig order
        z    (:,1) double {mustBeFinite} = zeros(0,1)   % z coordinates, Selig order
        name (1,1) string = ""                           % label (file header, plot legend)
    end
    % ----------------------------------------------------
    methods
        % constructor
        function obj = Airfoil_class(x,z,options)
            arguments
                x = zeros(0,1)
                z = zeros(0,1)
                options.name (1,1) string = ""
            end
            if nargin == 0; return; end
            if ischar(x) || isstring(x)                 % file name
                obj = Airfoil_class.fromFile(x);
                if options.name ~= ""; obj.name = options.name; end
                return
            end
            if isempty(z) && size(x,2) == 2             % [x z] array
                z = x(:,2);
                x = x(:,1);
            end
            obj.x    = x;
            obj.z    = z;
            obj.name = options.name;
            obj.check();
            obj = obj.toSelig();
        end
        % function to check
        function N = check(obj)
            N = numel(obj.x);
            if numel(obj.z) ~= N
                error('\nAirfoil_class: x and z must have the same number of points (%d vs %d)\n',N,numel(obj.z));
            end
            if N < 5
                error('\nAirfoil_class: at least 5 points are needed (got %d)\n',N);
            end
        end
        % leading edge index: the point farthest from the trailing edge midpoint
        function i = iLE(obj)
            xte   = (obj.x(1) + obj.x(end))/2;
            zte   = (obj.z(1) + obj.z(end))/2;
            [~,i] = max(hypot(obj.x - xte, obj.z - zte));
        end
        % leading edge between the points: spline through the airfoil
        function [xle,zle,sle,s] = leadingEdge(obj)
            % The LE is the point of the contour farthest from the TE
            % midpoint (as in XFOIL). The contour is splined in its arc
            % length s and the maximum distance is searched around the
            % farthest data point, so the LE does not have to be a data
            % point (e.g. a cambered NACA airfoil).
            s     = [0; cumsum(hypot(diff(obj.x),diff(obj.z)))];
            i     = obj.iLE();
            xte   = (obj.x(1) + obj.x(end))/2;
            zte   = (obj.z(1) + obj.z(end))/2;
            % data point already on the LE (contour normal to the TE direction
            % there, e.g. any airfoil built by CST2Airfoil): keep it exactly
            tang  = [obj.x(i+1) - obj.x(i-1), obj.z(i+1) - obj.z(i-1)];
            radv  = [obj.x(i) - xte, obj.z(i) - zte];
            if abs(tang*radv') <= 1e-9*norm(tang)*norm(radv)
                xle = obj.x(i); zle = obj.z(i); sle = s(i);
                return
            end
            ppx   = spline(s,obj.x);
            ppz   = spline(s,obj.z);
            dist2 = @(ss) -((ppval(ppx,ss) - xte).^2 + (ppval(ppz,ss) - zte).^2);
            sle   = fminbnd(dist2,s(max(i-1,1)),s(min(i+1,end)),optimset('TolX',1e-10*s(end)));
            xle   = ppval(ppx,sle);
            zle   = ppval(ppz,sle);
        end
        % upper and lower sides, both from the LE to the TE (LE point in both)
        function [xu,zu,xl,zl] = split(obj)
            [xle,zle,sle,s] = obj.leadingEdge();
            tol = 1e-9*s(end);                     % data points on the LE are not repeated
            up  = s < sle - tol;
            lo  = s > sle + tol;
            xu  = [xle; flip(obj.x(up))];   zu = [zle; flip(obj.z(up))];
            xl  = [xle; obj.x(lo)];         zl = [zle; obj.z(lo)];
        end
        % normalized airfoil: LE at (0,0), TE midpoint at (1,0)
        function [AFn,c,twist,xle,zle] = normalize(obj)
            % The physical section is   [x;z] = [xle;zle] + c*R(twist)*[x_n;z_n]
            % with R the nose-up rotation about the LE (rotateSection).
            %   c     - chord (LE to TE midpoint distance)
            %   twist - rad, positive nose up (TE below the LE)
            %   xle, zle - leading edge point (see leadingEdge)
            [xle,zle] = obj.leadingEdge();
            xte   = (obj.x(1) + obj.x(end))/2;
            zte   = (obj.z(1) + obj.z(end))/2;
            c     = hypot(xte - xle, zte - zle);
            twist = atan2(-(zte - zle), xte - xle);
            [xn,zn] = rotateSection((obj.x - xle)/c,(obj.z - zle)/c,-twist,0);
            AFn   = obj;
            AFn.x = xn;
            AFn.z = zn;
        end
        % same airfoil on a new chordwise distribution psi (x/c, both sides)
        function [AFr,zu_r,zl_r] = resample(obj,psi)
            % The sides are interpolated (pchip) in sqrt(x): z ~ sqrt(x) at
            % a round nose, so z is smooth in that variable.
            arguments
                obj
                psi (:,1) double {mustBeInRange(psi,0,1)} = cosineSpacing(101)
            end
            [xu,zu,xl,zl] = obj.split();
            zu_r = interpSide(xu,zu,psi);
            zl_r = interpSide(xl,zl,psi);
            AFr  = Airfoil_class([flip(psi); psi(2:end)],[flip(zu_r); zl_r(2:end)],name=obj.name);
        end
        % maximum thickness (t/c if the airfoil is normalized) and its x
        function [tmax,x_tmax] = thickness(obj)
            psi         = cosineSpacing(201);
            [~,zu,zl]   = obj.resample(psi);
            [tmax,i]    = max(zu - zl);
            x_tmax      = psi(i);
        end
        % plot
        function h = plot(obj,ax,varargin)
            % h = AF.plot()   or   AF.plot(ax,'r-','LineWidth',2)
            if nargin < 2 || isempty(ax); ax = gca; end
            h = plot(ax,obj.x,obj.z,varargin{:});
            axis(ax,'equal');
        end
        % write a Selig file (name in the first line)
        function write(obj,filename)
            fid = fopen(filename,'w');
            if fid < 0; error('\nAirfoil_class.write: cannot open %s\n',filename); end
            fprintf(fid,'%s\n',obj.name);
            fprintf(fid,'%12.8f %12.8f\n',[obj.x obj.z]');
            fclose(fid);
        end
    end
    % ----------------------------------------------------
    methods (Static)
        % read a Selig or Lednicer coordinate file
        function obj = fromFile(filename)
            % Selig:    name line, then x z from TE over the upper side to the LE and back
            % Lednicer: name line, "Nu Nl" line, upper side LE->TE, lower side LE->TE
            lines = strtrim(readlines(filename));
            lines = lines(lines ~= "");
            name  = "";
            data  = zeros(0,2);
            for i = 1:numel(lines)
                v = sscanf(lines(i),'%f');
                if numel(v) >= 2
                    data(end+1,:) = v(1:2)'; %#ok<AGROW>
                elseif i == 1
                    name = lines(i);
                end
            end
            if data(1,1) > 1.5 && data(1,2) > 1.5        % Lednicer: first row = number of points per side
                nu    = data(1,1);
                upper = data(2:1+nu,:);
                lower = data(2+nu:end,:);
                data  = [flip(upper,1); lower(2:end,:)];
            end
            obj = Airfoil_class(data(:,1),data(:,2),name=name);
        end
        % NACA 4 digit airfoil
        function obj = naca4(code,N,options)
            % AF = Airfoil_class.naca4("2412")            % 101 points per side
            % AF = Airfoil_class.naca4("0012",201,closedTE=false)
            arguments
                code (1,1) string
                N    (1,1) double {mustBeInteger,mustBeGreaterThan(N,4)} = 101
                options.closedTE (1,1) logical = true
            end
            d = char(code);
            if numel(d) ~= 4 || any(~isstrprop(d,'digit'))
                error('\nAirfoil_class.naca4: code must have 4 digits (got %s)\n',code);
            end
            m  = str2double(d(1))/100;            % max camber
            p  = str2double(d(2))/10;             % position of max camber
            t  = str2double(d(3:4))/100;          % thickness
            a4 = -0.1015;  if options.closedTE; a4 = -0.1036; end
            xc = cosineSpacing(N);
            yt = 5*t*(0.2969*sqrt(xc) - 0.1260*xc - 0.3516*xc.^2 + 0.2843*xc.^3 + a4*xc.^4);
            yc = zeros(N,1); dyc = zeros(N,1);
            if m > 0 && p > 0
                f       = xc < p;
                yc(f)   = m/p^2*(2*p*xc(f) - xc(f).^2);
                yc(~f)  = m/(1-p)^2*((1 - 2*p) + 2*p*xc(~f) - xc(~f).^2);
                dyc(f)  = 2*m/p^2*(p - xc(f));
                dyc(~f) = 2*m/(1-p)^2*(p - xc(~f));
            end
            th = atan(dyc);
            xu = xc - yt.*sin(th);  zu = yc + yt.*cos(th);
            xl = xc + yt.*sin(th);  zl = yc - yt.*cos(th);
            obj = Airfoil_class([flip(xu); xl(2:end)],[flip(zu); zl(2:end)],name="NACA "+code);
        end
    end
    % ----------------------------------------------------
    methods (Access = private)
        % make sure the points are in Selig order (upper side first)
        function obj = toSelig(obj)
            i = obj.iLE();
            if i == 1 || i == numel(obj.x)
                error('\nAirfoil_class: the points must go around the airfoil from TE to TE (the leading edge is point %d)\n',i);
            end
            % compare the two halves at the same x
            [xa,za] = deal(flip(obj.x(1:i)),flip(obj.z(1:i)));
            [xb,zb] = deal(obj.x(i:end),obj.z(i:end));
            xs = linspace(max(min(xa),min(xb)),min(max(xa),max(xb)),50)';
            if mean(interpSide(xa,za,xs,false)) < mean(interpSide(xb,zb,xs,false))
                obj.x = flip(obj.x);
                obj.z = flip(obj.z);
            end
        end
    end
end

% FUNCTIONs ---------------------------------------------
function zq = interpSide(x,z,xq,useSqrt)
% interpolate one side (from LE to TE) at xq, in sqrt(x) to follow a round nose
if nargin < 4; useSqrt = true; end
x0      = x(1);                                  % LE abscissa (0 for a normalized airfoil)
s       = x - x0;
sq      = xq - x0;
if useSqrt
    s  = sqrt(max(s,0));
    sq = sqrt(max(sq,0));
end
[s,ia]  = unique(s);                            % sorted, duplicates removed
zq      = interp1(s,z(ia),sq,'pchip','extrap');
end
