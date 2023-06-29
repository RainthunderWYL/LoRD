function NPCells = func_APNP_CellWithNullPoints2D( Bx,By )

dims = size(Bx);
% Locate cells include Null points
iix = 1:dims(2)-1;
iiy = 1:dims(1)-1;
[IX,IY] = meshgrid(iix,iiy);

% Bx
Bx_Cell = zeros(dims(1)-1,dims(2)-1,4);
Bx_Cell(:,:,1) = Bx(iiy,iix); % 
Bx_Cell(:,:,2) = Bx(iiy,iix+1); %
Bx_Cell(:,:,3) = Bx(iiy+1,iix); % 
Bx_Cell(:,:,4) = Bx(iiy+1,iix+1); %

b_BxNull = (max(Bx_Cell,[],3).*min(Bx_Cell,[],3))<=0;
clear Bx_Cell;

% By
By_Cell = zeros(dims(1)-1,dims(2)-1,4);
By_Cell(:,:,1) = By(iiy,iix); % 
By_Cell(:,:,2) = By(iiy,iix+1); %
By_Cell(:,:,3) = By(iiy+1,iix); % 
By_Cell(:,:,4) = By(iiy+1,iix+1); %

b_ByNull = (max(By_Cell,[],3).*min(By_Cell,[],3))<=0;
clear By_Cell;

b_CellWithNP = b_BxNull & b_ByNull;
clear b_BxNull b_ByNull;

% 
ix_NPs = IX(b_CellWithNP);
iy_NPs = IY(b_CellWithNP);
clear IX IY;
% 
ix_NPs = reshape(ix_NPs,numel(ix_NPs),1);
iy_NPs = reshape(iy_NPs,numel(iy_NPs),1);

NPCells = [iy_NPs, ix_NPs];
end

