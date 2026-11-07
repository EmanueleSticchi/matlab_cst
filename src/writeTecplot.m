function writeTecplot(filename,zones,options)
% WRITETECPLOT  Write surface zones to an ASCII Tecplot file.
%
%   writeTecplot('wing.dat',Wing.surface())
%   writeTecplot('wing.dat',[Wing.surface(); Fuselage.surface()],title="aircraft")
%
% One ordered zone (I x J, POINT format) for each element of zones, with
% the variables X, Y, Z. The file opens in Tecplot, ParaView and FieldView.
%
% INPUTS
%   filename - output file
%   zones    - struct array with fields name, X, Y, Z (same size matrices),
%              e.g. from the surface method of Wing3D_class, Nacelle_class,
%              Fuselage_class
%   title    - (name-value) dataset title
%
% See also: plotSurface
arguments
    filename (1,1) string
    zones    struct
    options.title (1,1) string = "matlab_cst"
end
fid = fopen(filename,'w');
if fid < 0; error('\nwriteTecplot: cannot open %s\n',filename); end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'TITLE = "%s"\n',options.title);
fprintf(fid,'VARIABLES = "X" "Y" "Z"\n');
for k = 1:numel(zones)
    [nj,ni] = size(zones(k).X);           % I runs along the columns (chordwise), J along the rows
    fprintf(fid,'ZONE T="%s", I=%d, J=%d, F=POINT\n',zones(k).name,ni,nj);
    X = zones(k).X'; Y = zones(k).Y'; Z = zones(k).Z';
    fprintf(fid,'%16.9e %16.9e %16.9e\n',[X(:) Y(:) Z(:)]');
end
end
