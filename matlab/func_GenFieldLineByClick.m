function [hfl, cors] = func_GenFieldLineByClick(axes1,X,Y,Z,B1,B2,B3,len)

cors = ginput3D(1);
xs = cors(1);
ys = cors(2);
zs = cors(3);

hfl = DrawSymmStreamLine3D(X,Y,Z,B1,B2,B3,xs,ys,zs,[1,len],[1,0,0],'-',1,axes1);

end