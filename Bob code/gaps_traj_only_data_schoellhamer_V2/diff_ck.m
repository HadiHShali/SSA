function diff_ck(a1,a2,vari)
    del=max(abs(a1(:)-a2(:)));
    disp(strcat("diff check ",vari," ",num2str(del)));
end