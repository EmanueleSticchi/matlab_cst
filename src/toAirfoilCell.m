function AFs = toAirfoilCell(Airfoil,N)
% TOAIRFOILCELL  N x 1 cell of Airfoil_class from the accepted section inputs.
%
%   AFs = toAirfoilCell(Airfoil,N)
%
% Airfoil can be
%   - one Airfoil_class, one N x 2 array [x z] or one file name (also in a
%     1 x 1 cell): the same airfoil is used in all the N sections
%   - a cell with N elements, each one of the above (NUIST style:
%     Airfoil{i}(:,1) = x/c, Airfoil{i}(:,2) = z/c)
%
% Used by the constructors of Wing3D_class and Nacelle_class.
if ~iscell(Airfoil)
    Airfoil = repmat({Airfoil},N,1);
elseif isscalar(Airfoil)
    Airfoil = repmat(Airfoil,N,1);
end
if numel(Airfoil) ~= N
    error('\nThe sections need %d airfoils (got %d)\n',N,numel(Airfoil));
end
AFs = cell(N,1);
for i = 1:N
    if isa(Airfoil{i},'Airfoil_class')
        AFs{i} = Airfoil{i};
    else
        AFs{i} = Airfoil_class(Airfoil{i});
    end
end
end
