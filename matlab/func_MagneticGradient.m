function DB = func_MagneticGradient(B1,B2,B3,x1,x2,x3,Method)
% Only for uniform Cartesian mesh
% Method:
% 0 (Central difference);
% 1 (Suppose B is on cell center, rebuild B on face and calculate gradients)

[nx2,nx1,nx3] = size(B1);
dx1 = mean(diff(x1));
dx2 = mean(diff(x2));
dx3 = mean(diff(x3));

if Method == 0
    [DB{1},DB{4},DB{7}] = gradient(B1,dx1,dx2,dx3);
    [DB{2},DB{5},DB{8}] = gradient(B2,dx1,dx2,dx3);
    [DB{3},DB{6},DB{9}] = gradient(B3,dx1,dx2,dx3);
elseif Method == 1
    % Rebuild B on cell faces
    B1_f = zeros(nx2,nx1+1,nx3);
    B1_f(:,2:end-1,:) = 0.5*(B1(:,1:end-1,:) + B1(:,2:end,:));
    B1_f(:,1,:) = 2*B1(:,1,:)-B1_f(:,2,:);
    B1_f(:,end,:) = 2*B1(:,end,:)-B1_f(:,end-1,:);

    B2_f = zeros(nx2+1,nx1,nx3);
    B2_f(2:end-1,:,:) = 0.5*(B2(1:end-1,:,:) + B2(2:end,:,:));
    B2_f(1,:,:) = 2*B2(1,:,:)-B2_f(2,:,:);
    B2_f(end,:,:) = 2*B2(end,:,:)-B2_f(end-1,:,:);

    B3_f = zeros(nx2,nx1,nx3+1);
    B3_f(:,:,2:end-1) = 0.5*(B3(:,:,1:end-1) + B3(:,:,2:end));
    B3_f(:,:,1) = 2*B3(:,:,1)-B3_f(:,:,2);
    B3_f(:,:,end) = 2*B3(:,:,end)-B3_f(:,:,end-1);

    DB{1} = diff(B1_f,1,2)/dx1;
    DB{5} = diff(B2_f,1,1)/dx2;
    DB{9} = diff(B3_f,1,3)/dx3;
    DivB = DB{1} + DB{5} + DB{9};

else
    error('func_MagneticGradient: Wrong type of Method (0 or 1)!!!');
end
end