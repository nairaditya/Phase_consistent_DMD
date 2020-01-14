function out = edgeDistance(L)
Lt = ones(size(L,1)+2,size(L,2)+2);
Lt(2:end-1,2:end-1) = double(L);
Lt(2:end-1,[1,end]) = double(L(:,[1,end]));
Lt([1,end],2:end-1) = double(L([1,end],:));

Lt = colfilt(Lt,[3 3],'sliding',@sum);
Lt(Lt==4) = 6;
Lt(Lt==2) = 3;
Lt(Lt==1) = 0;
Lt(round(Lt)~=0) = Lt(round(Lt)~=0)./3-1;

Lt(:,[1,end]) = Lt(:,[2,end-1]);
Lt([1,end],:) = Lt([2,end-1],:);


out = zeros(size(Lt));
outLast = out;
cont = 1;
i = 0;
while cont
    
    i = i + 1;
    out(Lt == 1) = i;
    Lt(Lt == 1) = 0;
    Lt(Lt>1) = 1;
    Lt = colfilt(Lt,[3 3],'sliding',@sum);
    Lt(Lt==4) = 6;
    Lt(Lt==2) = 3;
    Lt(Lt==1) = 0;
    Lt(round(Lt)~=0) = Lt(round(Lt)~=0)./3-1;

    
    Lt(:,[1,end]) = Lt(:,[2,end-1]);
    Lt([1,end],:) = Lt([2,end-1],:);
    
    if all(out(:)==outLast(:))
        cont = 0;
    else
        outLast = out;
    end
    
end

out = out(2:end-1,2:end-1);

end