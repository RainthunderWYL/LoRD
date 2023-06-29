function B = func_MapSortedIndex(A,I)
[ny,nx] = size(I);
[i,j] = meshgrid(1:nx,1:ny);
II = reshape(ny*(i-1)+I,[numel(A),1]);
B = reshape(A,[numel(A) 1]);
B = reshape(B(II),[ny,nx]);
end