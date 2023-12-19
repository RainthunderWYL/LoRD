function [RMarker, RGridCount] = Tool_Cluster_RegionGrow(xi,yi,zi,Gx,Gy,Gz)
% Clustering grids into connected segmentations by the region growing method. 
% Here, two grids are connected if they are vertices of one cell.
%
% [RMarker, RGridCount] = Tool_Cluster_RegionGrow(xi,yi,zi,Gx,Gy,Gz)
%
% Inputs: 
%   xi, yi, and zi are 1D arrays of grid coordinates; 
%   Gx, Gy, and Gz are mesh coordinates on three directions
% Outputs: 
%   RMarker is a 1D array of natural number, labeling the region
%     markers of inputed grids. RMarker is of the same size as xi. 
%   RGridCount has two columns and N rows. N is the number of recognized 
%     segmentations. The first column is region marker, while the second 
%     column stores the grid numbers of every segments.

dx = mean(diff(Gx));
dy = mean(diff(Gy));
dz = mean(diff(Gz));
XLIM = [min(xi(:)),max(xi(:))];
YLIM = [min(yi(:)),max(yi(:))];
ZLIM = [min(zi(:)),max(zi(:))];

nx = round(diff(XLIM)/dx)+1;
ny = round(diff(YLIM)/dy)+1;
nz = round(diff(ZLIM)/dz)+1;

xi_int = round((xi-XLIM(1))/dx)+1;
yi_int = round((yi-YLIM(1))/dy)+1;
zi_int = round((zi-ZLIM(1))/dz)+1;

% Map to 1D grid
infunc_sub2ind = @(y,x,z) y + (x - 1)*ny + (z - 1)*nx*ny;
% Jp_int = infunc_sub2ind(yi_int,xi_int,zi_int);
Jp_int = sub2ind([ny,nx,nz],yi_int,xi_int,zi_int);

% Gen neighboor M
for k = -1:1:1
    for i = -1:1:1
        for j = -1:1:1
            J = j+2 + 3*(i+2-1) + 9*(k+2-1);
            Mnb(J,:) = [i,j,k];
        end
    end
end
Mnb(14,:) = [];

% Initialization
X_seeds = Jp_int;
NRegion = 0;
RMarker = zeros(length(X_seeds),1);
RGridCount = [];

count_seeds = length(X_seeds);
while count_seeds ~= 0
    % First row is a new seed
    ijk0 = X_seeds(1);
    NRegion = NRegion + 1;

    % set region marker
    RMarker(Jp_int == ijk0) = NRegion;

    % Remove this seed
    X_seeds(X_seeds == ijk0) = [];

    % Start growth of region
    ijk = ijk0;
    while isempty(ijk) == 0
        ijk_tmp = [];
        for iedge = 1:length(ijk)
            % Get index of its 26 neighboors
            [iy,ix,iz] = ind2sub([ny,nx,nz],ijk(iedge));
            M0 = repmat([iy,ix,iz],26,1);
            M1 = M0 + Mnb;
            Jnb = infunc_sub2ind(M1(:,1),M1(:,2),M1(:,3));
            [IsEdge,Loc] = ismember(Jnb,X_seeds);

            if sum(IsEdge) ~= 0
                % Save new edges to cache
                ijk_tmp_i = Jnb(IsEdge);
                ijk_tmp = [ijk_tmp;ijk_tmp_i];

                % Remove them in X_seeds
                X_seeds(Loc(IsEdge)) = [];

                % set region marker
                [Lia,Loc_Jp_int] = ismember(ijk_tmp_i,Jp_int);
                RMarker(Loc_Jp_int(Lia)) = NRegion;
            end
        end
        ijk = ijk_tmp;
    end
    count_seeds_new = length(X_seeds);
    RGridCount = [RGridCount; count_seeds - count_seeds_new];
    count_seeds = count_seeds_new;
end

CountRegions = (1:length(RGridCount))';
[RGridCount,MI_tmp] = sort(RGridCount);
RGridCount = flipud([CountRegions(MI_tmp) RGridCount]);

end