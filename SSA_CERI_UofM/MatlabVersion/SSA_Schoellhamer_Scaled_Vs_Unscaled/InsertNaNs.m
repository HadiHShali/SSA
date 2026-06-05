function XNaN=InsertNaNs(X,nanindx)
    XNaN=X;
    XNaN(nanindx,:)=NaN;
end