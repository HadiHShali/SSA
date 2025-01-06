
function antidiagmean=myantidiagmean(InMat) %takes longer, don't use
    tamano=size(InMat);
    tmp=NaN(sum(tamano)-1,tamano(2));
    tmp(1:tamano(1),1:tamano(2))=InMat;
    for colnum=2:tamano(2)
        tmp(:,colnum)=circshift(tmp(:,colnum),colnum-1);
    end
    antidiagmean=mean(tmp,2,"omitnan");
end