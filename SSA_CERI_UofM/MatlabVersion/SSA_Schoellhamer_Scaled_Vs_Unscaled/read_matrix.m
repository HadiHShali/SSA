function [data datasize]=read_matrix(filename)
%reads NUMERICAL file with data in N columns (freeform)
    data=load(filename);
    datasize=size(data);
end