function NPInfo=func_ANP_CleanAdjacentPoints(NPInfo_All,x1,x2,x3)
DL = [mean(diff(x1)),mean(diff(x2)),mean(diff(x3))];
dl = mean(DL);

% Clean Near points
num_NPs = size(NPInfo_All,1);
i = 1;
while i < num_NPs
    r_ref = NPInfo_All(i,1:3);

    x = NPInfo_All(:,1);
    y = NPInfo_All(:,2);
    z = NPInfo_All(:,3);

    r_i = sqrt((x-r_ref(1)).^2 + (y-r_ref(2)).^2 + (z-r_ref(3)).^2);
    b = r_i<=dl;
    b(i) = false;
    NPInfo_All(b,:) = [];
    i = i+1;
    num_NPs = size(NPInfo_All,1);
end
NPInfo = NPInfo_All;
end