function [F,x1,x2,x3] = Tool_LoadData_Bin(filename,dim,precision,varargin)

fp = fopen(filename,'r');
F = permute(reshape(fread(fp,precision),dim),[2,1,3]);
fclose(fp);

if ~isempty(varargin)
    domain = varargin{1};
    x1 = linspace(domain(1),domain(2),dim(1));
    x2 = linspace(domain(3),domain(4),dim(2));
    x3 = linspace(domain(5),domain(6),dim(3));
else
    x1 = 0:dim(1)-1;
    x2 = 0:dim(2)-1;
    x3 = 0:dim(3)-1;
end
end
