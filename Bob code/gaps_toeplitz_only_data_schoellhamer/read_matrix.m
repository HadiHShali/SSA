function [data]=read_matrix(filename)
%reads NUMERICAL file with data in N columns (freeform)
    data=load(filename);
end