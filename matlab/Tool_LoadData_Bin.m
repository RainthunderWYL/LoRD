function [F,x1,x2,x3] = Tool_LoadData_Bin(filename,n1,n2,n3,precision,varargin)
% fB#: file name of B#
% n#: number of grids of x#
% precision: double, float32, ... (help fwrite for details)
% Domain: [x1,x2,y1,y2,z1,z2]

fp = fopen(filename,'r');
F = permute(reshape(fread(fp,precision),[n1,n2,n3]),[2,1,3]);
fclose(fp);

if ~isempty(varargin)
    domain = varargin{1};
    x1 = linspace(domain(1),domain(2),n1);
    x2 = linspace(domain(3),domain(4),n2);
    x3 = linspace(domain(5),domain(6),n3);
else
    x1 = 0:n1-1;
    x2 = 0:n2-1;
    x3 = 0:n3-1;
end
end