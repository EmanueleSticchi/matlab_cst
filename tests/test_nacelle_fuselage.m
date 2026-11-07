function tests = test_nacelle_fuselage
% TEST_NACELLE_FUSELAGE  Nacelle and fuselage fit / build.
tests = functiontests(localfunctions);
end

function setupOnce(~)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','src'));
end

function test_nacelleRoundTrip(tc)
% CST -> nacelle -> CST gives back every parameter
CST0       = struct();
CST0.Bu    = [0.17 0.34 0.18 0.33; 0.18 0.30 0.20 0.30; 0.20 0.28 0.22 0.25];
CST0.Bl    = [-0.17 0.04 -0.10 0.04; -0.18 0.02 -0.12 0.03; -0.20 0.00 -0.13 0.01];
CST0.Rle   = [1.0; 1.05; 1.1];
CST0.c     = [3.0; 2.9; 2.8];
CST0.Xle   = [0; -0.05; -0.15];
CST0.twist = deg2rad([0; -1; -2]);
Nac = CST2Nacelle(CST0,phi=linspace(0,pi,9));
CST = Nacelle2CST(Nac,nx=3,ny=2,ny_plan=2);
verifyEqual(tc,CST.Bu,CST0.Bu,'AbsTol',1e-8);
verifyEqual(tc,[CST.Rle CST.c CST.Xle CST.twist],[CST0.Rle CST0.c CST0.Xle CST0.twist],'AbsTol',1e-8);
end

function test_nacelleSmoothAtSymmetryPlane(tc)
% the profile quantities are even in phi: no ridge on crown and keel
CST0 = struct('Bu',[0.2 0.3 0.2 0.3; 0.2 0.3 0.2 0.3],'Bl',[-0.2 0 -0.1 0; -0.2 0 -0.1 0], ...
              'Rle',[1; 1.2],'c',3);
d    = 1e-4;
Nac  = CST2Nacelle(CST0,phi=[0 d pi-d pi]);
verifyLessThan(tc,abs(Nac.Rle(2) - Nac.Rle(1))/d,1e-3);   % dRle/dphi = 0 at the crown
verifyLessThan(tc,abs(Nac.Rle(4) - Nac.Rle(3))/d,1e-3);   % ... and at the keel
end

function test_nacelleAxisymmetric(tc)
Nac = Nacelle_class(0,1,3,0,Airfoil_class.naca4("4412"));
CST = Nacelle2CST(Nac,nx=6);
verifyEqual(tc,[CST.ny CST.ny_plan],[0 0]);
zones = CST2Nacelle(CST).surface(41,25);
r = hypot(zones(1).Y,zones(1).Z);
verifyEqual(tc,max(std(r,0,1)),0,'AbsTol',1e-12);         % same radius all around
end

function test_fuselageRoundTrip(tc)
% CST -> fuselage -> CST gives back every parameter
CST0 = struct('L',30,'N1',0.5,'N2',1,'Ahu',[0.10 0.08 0.09 0.06]','dzhu',0.002, ...
              'Ahl',[0.09 0.07 0.08 0.05]','dzhl',0.001,'Aw',[0.10 0.09 0.08 0.07]','dzw',0.002, ...
              'Bzm',[0 -0.01 0.02]','BNU',[0.5 0.6]','BNL',[0.4 0.5]');
F   = CST2Fuselage(CST0,npsi=41);
CST = Fuselage2CST(F,n=3,nz=2,nN=1);
verifyEqual(tc,[CST.Ahu CST.Ahl CST.Aw],[CST0.Ahu CST0.Ahl CST0.Aw],'AbsTol',1e-9);
verifyEqual(tc,[CST.dzhu CST.dzhl CST.dzw],[CST0.dzhu CST0.dzhl CST0.dzw],'AbsTol',1e-10);
verifyEqual(tc,CST.Bzm,CST0.Bzm,'AbsTol',1e-9);
verifyEqual(tc,[CST.BNU CST.BNL],[CST0.BNU CST0.BNL],'AbsTol',1e-9);
end

function test_fuselageFromSections(tc)
% elliptic and boxy sections: the lobe exponents are recovered
x  = linspace(0,10,5)';
th = linspace(0,2*pi,201);
[Ys,Zs] = deal(cell(5,1));
for k = 1:5
    Ys{k} = 1.5*cos(th);
    up    = sin(th) >= 0;
    Zs{k} = zeros(size(th));
    Zs{k}(up)  = 0.2 + 1.0*(sin(th(up)).^2).^0.5;     % NU = 0.5
    Zs{k}(~up) = 0.2 - 0.8*(sin(th(~up)).^2).^0.2;    % NL = 0.2
end
F = Fuselage_class.fromSections(x,Ys,Zs);
verifyEqual(tc,F.NU,0.5*ones(5,1),'AbsTol',1e-3);
verifyEqual(tc,F.NL,0.2*ones(5,1),'AbsTol',1e-3);
verifyEqual(tc,[F.w F.zm F.zc F.zk],repmat([1.5 0.2 1.2 -0.6],5,1),'AbsTol',1e-9);
end
