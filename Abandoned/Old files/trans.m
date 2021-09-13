function TB=trans(r1,r2,gammax,nalt,nvar,IB)
m=size(gammax);
TB=zeros(m+1, 1+nvar+2*m+nalt);
for i=1:length(IB)
    TB(i,IB(i))=1;
    if IB(i)>1+nvar+m
        TB(i,IB(i)-m-nalt)=-1;
    else
        TB(i,IB(i)+m+nalt)=-1;
        TB(end,IB(i)+m+nalt)=1;
    end
    TB(:,r1)=gammax;
    r20=r2>0;
    TB(:,r2(r20))=
end


end