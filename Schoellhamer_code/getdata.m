function x_in_cont=getdata()

%read a trend and jumps removed gappy hector file from Hadi and return with NaNs in gaps
%0=ew,1=NS,2=ud, ordered but not uniformly spaced, column vector
filename_root='C:\Users\GeodesyLab\OneDrive - The University of Memphis\Desktop\SSA\Schoellhamer code';
sitename='1LSU';
compname='E';
compnum='1';
filename='1LSU_E/1LSU_0_SSA.dat'; %(all real) input data has gaps
full_filename=strcat(filename_root,"/",sitename,"_",compnum,"_SSA.dat");

input_data_discont=load(full_filename);

x_in_discont=input_data_discont(:,3); %discontinuous in time data vector, cols: fracyr, mjd, position (1 comp/file)
x_in_discont_indx=input_data_discont(:,1)-input_data_discont(1,1)+1; %this has indices of days with data
N=input_data_discont(end,1)-input_data_discont(1,1)+1; %num days from first to last, >= to number days with data
x_in_cont=nan(N,1); %to make a copy of input vector that is continuous, fill with NaNs, will remain in holes
x_in_cont(x_in_discont_indx)=x_in_discont; %move in data, discont_in_indx has index
end

