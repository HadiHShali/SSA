clear
close all

ts=getdata_bob;
m=500; %set embedding dimension to 3.5 years (rounded to integer),
%gives number of columns in trajectory matrix and number RCs we will obtain
f=0;

[xmean,xstd,xk,ro,lam,xk_]=vssanan_bob(ts,m,f);
xf=rc_bob(xk,ro);
xf_=rc_bob(xk_,ro);

nan_indx=find(isnan(ts))';
RC=(sum(xf)*xstd)+xmean;
RC_=(sum(xf_)*xstd)+xmean;
RC(nan_indx)=nan;
RC_(nan_indx)=nan;

max(abs(ts-RC'))
max(abs(ts-RC_'))

plot(ts,'b')
hold on
plot(RC,'g')
plot(RC_,'r--')
legend('Original Timeseries', 'Scaled Reconstructed Signal', 'Unscaled Reconstructed Signal')
