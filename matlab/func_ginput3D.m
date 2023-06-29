function Cor3D = func_ginput3D(n)
    datacursormode on;
    dcmObj = datacursormode(gcf);

    set(dcmObj,'SnapToDataVertex','off','enable', 'on','DisplayStyle','window');

    Cor3D = zeros(n,3);
    for i = 1:n
        disp('Click to get a coordinate');
        waitforbuttonpress;
        point = dcmObj.getCursorInfo;
        Cor3D(i,1) = point.Position(1);
        Cor3D(i,2) = point.Position(2);
        Cor3D(i,3) = point.Position(3);
    end
end