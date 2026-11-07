function tests = test_wing
% TEST_WING  Wing3D_class, Wing2CST and CST2Wing.
tests = functiontests(localfunctions);
end

function setupOnce(~)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','src'));
end

function CST = crankedCST()
% a cranked wing written by hand: linear in each of the two segments
CST            = struct();
CST.b          = 15;
CST.eta_breaks = 0.4;
CST.Bu         = [0.20 0.25 0.18 0.22; 0.18 0.22 0.16 0.20; 0.14 0.20 0.12 0.18];
CST.Bl         = [-0.20 -0.10 -0.15 -0.08; -0.18 -0.08 -0.12 -0.05; -0.14 -0.03 -0.10 -0.02];
CST.c          = [6; 4; 1.5];
CST.Xle        = [0; 2; 7];
CST.Zle        = [0; 0.3; 1.0];
CST.twist      = deg2rad([3; 1; -1]);
end

function test_roundTrip(tc)
% CST -> wing -> CST gives back every parameter
CST0 = crankedCST();
W    = CST2Wing(CST0,eta=linspace(0,1,15));
CST  = Wing2CST(W,nx=3,ny=1,ny_plan=1,eta_breaks=0.4);
verifyEqual(tc,CST.Bu,CST0.Bu,'AbsTol',1e-8);
verifyEqual(tc,CST.Bl,CST0.Bl,'AbsTol',1e-8);
verifyEqual(tc,[CST.c CST.Xle CST.Zle CST.twist],[CST0.c CST0.Xle CST0.Zle CST0.twist],'AbsTol',1e-8);
end

function test_planformGeometry(tc)
% sweep, area and aspect ratio of a straight tapered wing
AF = Airfoil_class.naca4("0012");
W  = Wing3D_class([0 1]',10,[4 2]',[0 0]',AF,sweep=deg2rad(30));
verifyEqual(tc,W.S,30,'AbsTol',1e-12);
verifyEqual(tc,W.AR,2*100/30,'AbsTol',1e-12);
verifyEqual(tc,rad2deg(W.sweepLE),30,'AbsTol',1e-10);
verifyEqual(tc,W.Xle(2),10*tan(deg2rad(30)),'AbsTol',1e-12);
end

function test_fromSections(tc)
% sections cut from a known wing give back its planform
W  = CST2Wing(crankedCST(),eta=linspace(0,1,9),npsi=121);
XZ = cell(numel(W.eta),1);
for k = 1:numel(W.eta)
    [Xu,Zu,Xl,Zl] = W.section3D(k,cosineSpacing(121));
    XZ{k} = [flip(Xu') flip(Zu'); Xl(2:end)' Zl(2:end)'];
end
W2 = Wing3D_class.fromSections(W.Yle,XZ,eta_breaks=0.4);
verifyEqual(tc,W2.c,W.c,'AbsTol',1e-9);
verifyEqual(tc,W2.twist,W.twist,'AbsTol',1e-9);
verifyEqual(tc,[W2.Xle W2.Zle],[W.Xle W.Zle],'AbsTol',1e-9);
end

function test_breakNeeded(tc)
% without the break a single straight panel cannot follow the kink
W   = CST2Wing(crankedCST(),eta=linspace(0,1,15));
C1  = Wing2CST(W,nx=3,ny=1,ny_plan=1,eta_breaks=0.4);
C0  = Wing2CST(W,nx=3,ny=1,ny_plan=1,eta_breaks=[]);
verifyLessThan(tc,max(abs(C1.err.c)),1e-8);
verifyGreaterThan(tc,max(abs(C0.err.c)),0.1);
end

function test_normalizeSectionsKeepsGeometry(tc)
% moving the airfoil normalization into the planform keeps the 3D section
AF = Airfoil_class.naca4("4415");                   % geometric LE not at x = 0
W  = Wing3D_class([0 1]',10,[3 2]',deg2rad([4 -2])',AF,Xle=[0 1]',Zle=[0 0.5]');
Wn = W.normalizeSections();
for k = 1:2
    A = W.Airfoil{k};  B = Wn.Airfoil{k};
    [x1,z1] = rotateSection(A.x*W.c(k), A.z*W.c(k), W.twist(k), 0.25*W.c(k));
    [x2,z2] = rotateSection(B.x*Wn.c(k),B.z*Wn.c(k),Wn.twist(k),0.25*Wn.c(k));
    verifyEqual(tc,[x1 + W.Xle(k), z1 + W.Zle(k)],[x2 + Wn.Xle(k), z2 + Wn.Zle(k)],'AbsTol',1e-12);
end
end
