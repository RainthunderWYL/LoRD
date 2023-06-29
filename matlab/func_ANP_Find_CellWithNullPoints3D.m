function NPCells = func_ANP_Find_CellWithNullPoints3D( Bx,By,Bz )

dims = size(Bx);
% Locate cells include Null points
iix = 1:dims(2)-1;
iiy = 1:dims(1)-1;
iiz = 1:dims(3)-1;
[IX,IY,IZ] = meshgrid(iix,iiy,iiz);

% Bx
Bx_Cell = zeros(dims(1)-1,dims(2)-1,dims(3)-1,8);
Bx_Cell(:,:,:,1) = Bx(iiy,iix,iiz); % 
Bx_Cell(:,:,:,2) = Bx(iiy,iix+1,iiz); %
Bx_Cell(:,:,:,3) = Bx(iiy+1,iix,iiz); % 
Bx_Cell(:,:,:,4) = Bx(iiy+1,iix+1,iiz); %
Bx_Cell(:,:,:,5) = Bx(iiy,iix,iiz+1); % 
Bx_Cell(:,:,:,6) = Bx(iiy,iix+1,iiz+1); %
Bx_Cell(:,:,:,7) = Bx(iiy+1,iix,iiz+1); % 
Bx_Cell(:,:,:,8) = Bx(iiy+1,iix+1,iiz+1); %

b_BxNull = (max(Bx_Cell,[],4).*min(Bx_Cell,[],4))<=0;
clear Bx_Cell;

% By
By_Cell = zeros(dims(1)-1,dims(2)-1,dims(3)-1,8);
By_Cell(:,:,:,1) = By(iiy,iix,iiz); % 
By_Cell(:,:,:,2) = By(iiy,iix+1,iiz); %
By_Cell(:,:,:,3) = By(iiy+1,iix,iiz); % 
By_Cell(:,:,:,4) = By(iiy+1,iix+1,iiz); %
By_Cell(:,:,:,5) = By(iiy,iix,iiz+1); % 
By_Cell(:,:,:,6) = By(iiy,iix+1,iiz+1); %
By_Cell(:,:,:,7) = By(iiy+1,iix,iiz+1); % 
By_Cell(:,:,:,8) = By(iiy+1,iix+1,iiz+1); %

b_ByNull = (max(By_Cell,[],4).*min(By_Cell,[],4))<=0;
clear By_Cell;

% Bz
Bz_Cell = zeros(dims(1)-1,dims(2)-1,dims(3)-1,8);
Bz_Cell(:,:,:,1) = Bz(iiy,iix,iiz); % 
Bz_Cell(:,:,:,2) = Bz(iiy,iix+1,iiz); %
Bz_Cell(:,:,:,3) = Bz(iiy+1,iix,iiz); % 
Bz_Cell(:,:,:,4) = Bz(iiy+1,iix+1,iiz); %
Bz_Cell(:,:,:,5) = Bz(iiy,iix,iiz+1); % 
Bz_Cell(:,:,:,6) = Bz(iiy,iix+1,iiz+1); %
Bz_Cell(:,:,:,7) = Bz(iiy+1,iix,iiz+1); % 
Bz_Cell(:,:,:,8) = Bz(iiy+1,iix+1,iiz+1); %

b_BzNull = (max(Bz_Cell,[],4).*min(Bz_Cell,[],4))<=0;
clear Bz_Cell;

b_CellWithNP = b_BxNull & b_ByNull & b_BzNull;
clear b_BxNull b_ByNull b_BzNull;

% 
ix_NPs = IX(b_CellWithNP);
iy_NPs = IY(b_CellWithNP);
iz_NPs = IZ(b_CellWithNP);
clear IX IY IZ;
% 
ix_NPs = reshape(ix_NPs,numel(ix_NPs),1);
iy_NPs = reshape(iy_NPs,numel(iy_NPs),1);
iz_NPs = reshape(iz_NPs,numel(iz_NPs),1);

NPCells = [iy_NPs, ix_NPs, iz_NPs];
end

