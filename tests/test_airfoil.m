function tests = test_airfoil
% TEST_AIRFOIL  Airfoil_class, Airfoil2CST and CST2Airfoil.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','src'));
tc.TestData.data = fullfile(fileparts(mfilename('fullpath')),'..','examples','data','naca0012.dat');
end

function test_roundTrip(tc)
% CST -> airfoil -> CST gives back the coefficients
CST0 = struct('Au',[0.2 0.15 0.18 0.12 0.1]','Al',[-0.2 -0.1 -0.05 0.02 0.05]','dzu',0.002,'dzl',-0.002);
AF   = CST2Airfoil(CST0,151);
CST  = Airfoil2CST(AF,n=4);
verifyEqual(tc,CST.Au,CST0.Au,'AbsTol',1e-10);
verifyEqual(tc,CST.Al,CST0.Al,'AbsTol',1e-10);
verifyEqual(tc,[CST.dzu CST.dzl],[CST0.dzu CST0.dzl],'AbsTol',1e-12);
end

function test_normalizeRotatedScaled(tc)
% a scaled, rotated, moved airfoil is normalized back and fitted exactly
CST0 = struct('Au',[0.2 0.15 0.18 0.12 0.1]','Al',[-0.2 -0.1 -0.05 0.02 0.05]');
AF   = CST2Airfoil(CST0,101);
[x,z] = rotateSection(2*AF.x,2*AF.z,deg2rad(5),0.5);
AFp  = Airfoil_class(x + 3,z - 1);
[~,c,twist,xle,zle] = AFp.normalize();
verifyEqual(tc,[c rad2deg(twist)],[2 5],'AbsTol',1e-10);
[xl0,zl0] = rotateSection(0,0,deg2rad(5),0.5);
verifyEqual(tc,[xle zle],[xl0 + 3, zl0 - 1],'AbsTol',1e-10);
CST  = Airfoil2CST(AFp,n=4);
verifyEqual(tc,CST.Au,CST0.Au(:),'AbsTol',1e-9);
end

function test_fileAndOrder(tc)
% Selig file read; reversed point order is put back in Selig order
AF  = Airfoil_class(tc.TestData.data);
verifyGreaterThan(tc,numel(AF.x),100);
AFr = Airfoil_class(flip(AF.x),flip(AF.z));
verifyEqual(tc,AFr.x,AF.x);
verifyEqual(tc,AF.thickness(),0.12,'AbsTol',1e-3);
end

function test_fitAccuracy(tc)
% Kulfan 2008: from order ~6 the fit is within wind tunnel tolerance
AF  = Airfoil_class(tc.TestData.data);
CST = Airfoil2CST(AF,n=6);
verifyLessThan(tc,CST.err.max,2e-4);
verifyEqual(tc,CST.Al(1),-CST.Au(1));        % LE continuity
CST = Airfoil2CST(Airfoil_class.naca4("2412"),n=8);
verifyLessThan(tc,CST.err.rms,2e-4);
end

function test_lednicer(tc)
% Lednicer file: same airfoil as Selig
AF  = Airfoil_class.naca4("2412",41);
[xu,zu,xl,zl] = AF.split();
f   = [tempname '.dat'];
fid = fopen(f,'w');
fprintf(fid,'NACA 2412 lednicer\n%d. %d.\n\n',numel(xu),numel(xl));
fprintf(fid,'%.8f %.8f\n',[xu zu]');
fprintf(fid,'\n');
fprintf(fid,'%.8f %.8f\n',[xl zl]');
fclose(fid);
AFl = Airfoil_class(f);
delete(f);
verifyEqual(tc,AFl.thickness(),AF.thickness(),'AbsTol',1e-6);
end
