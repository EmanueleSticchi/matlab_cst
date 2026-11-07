function tests = test_core
% TEST_CORE  Class function, Bernstein basis, CST curve and spanwise basis.
tests = functiontests(localfunctions);
end

function setupOnce(~)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','src'));
end

function test_partitionOfUnity(tc)
% the Bernstein polynomials of any order sum to 1 (Kulfan 2008, Fig. 6)
t = linspace(0,1,37)';
for n = 0:8
    verifyEqual(tc,sum(bernsteinBasis(t,n),2),ones(size(t)),'AbsTol',1e-12);
end
end

function test_classFunctionEllipse(tc)
% N1 = N2 = 0.5 and S = 1 is a circle of diameter 1 (Kulfan 2008, Eq. 6)
psi  = linspace(0,1,51)';
zeta = classFunction(psi,0.5,0.5);
verifyEqual(tc,(psi - 0.5).^2 + zeta.^2,0.25*ones(size(psi)),'AbsTol',1e-12);
end

function test_leadingEdgeRadius(tc)
% S(0) = A(1) = sqrt(2 r_LE/c): the curve near psi = 0 is a circle of radius A(1)^2/2
A    = [0.2 0.1 0.15]';
rLE  = A(1)^2/2;
psi  = 1e-6;
zeta = cstCurve(psi,A,0,0.5,1);
verifyEqual(tc,zeta,sqrt(2*rLE*psi),'RelTol',1e-4);
end

function test_fitCstCurveExact(tc)
% a CST curve is recovered exactly (the fit is linear)
psi  = cosineSpacing(60);
A    = [0.18 0.12 0.2 0.1 0.15]';
zeta = cstCurve(psi,A,0.003,0.5,1);
[Af,dz] = fitCstCurve(psi,zeta,4,0.5,1);
verifyEqual(tc,Af,A,'AbsTol',1e-10);
verifyEqual(tc,dz,0.003,'AbsTol',1e-12);
end

function test_spanwiseBasisContinuity(tc)
% the piecewise basis is continuous at the breaks and interpolates the ends
breaks = [0.3 0.7];
ny     = 2;
[B,eta_ctrl] = spanwiseBasis([0; 0.3-1e-12; 0.3; 0.3+1e-12; 1],breaks,ny);
verifyEqual(tc,size(B,2),3*ny + 1);
verifyEqual(tc,B(2,:),B(4,:),'AbsTol',1e-9);
verifyEqual(tc,B(1,1),1);                 % root = first control value
verifyEqual(tc,B(end,end),1);             % tip  = last control value
verifyEqual(tc,eta_ctrl([1 3 5 7]),[0 0.3 0.7 1]','AbsTol',1e-12);
verifyEqual(tc,sum(spanwiseBasis(linspace(0,1,41)',breaks,ny),2),ones(41,1),'AbsTol',1e-12);
end

function test_fitSpanwiseTooFewStations(tc)
% a clear error when a segment cannot determine the order
verifyError(tc,@() fitSpanwise([0 0.5 1]',(1:3)',0.5,2),?MException);
end
