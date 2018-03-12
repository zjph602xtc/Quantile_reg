function [ estimate,j, A,TB,IB ] = altermyqr3(xnew,A,TB,IB,c,nvar,y)


m=length(y);
xnew(y<0)=-xnew(y<0);
bi=inv(A(:,IB));
bix=bi*xnew;
if (any(IB==3))
    flag=-1;
    TB(:,4)=[bix; 0];
    %      TB(:,4)=[bix; -c(IB)'*bix];
    TB=TB(:,[1:(end-1) 3 end]);
    IB(IB==3)=3+2*nvar+2*m;
    %     TB(:,3)=[-bix; 0];
    TB(:,3)=[-bix; 0];
    %     TB(end,[1:3 5:(end-2) end])=-TB(find(IB==[3+2*nvar+2*m]),[1:3 5:(end-2) end]);
    TB(end,[1:(end-2) end])=-TB(find(IB==[3+2*nvar+2*m]),[1:(end-2) end]);
    %         TB(end,[1:end])=-TB(find(IB==[3+2*nvar+2*m]),[1:end]);
    %         TB(end,[1:end])=-TB(find(IB==[3+2*nvar+2*m]),[1:end]);
    A(:,3:4)=[-xnew xnew ];
else
    flag=1;
    TB(:,3)=[bix; 0];
    %     TB(:,3)=[bix; -c(IB)'*bix];
    TB=TB(:,[1:(end-1) 4 end]);
    IB(IB==4)=3+2*nvar+2*m;
    TB(:,4)=[-bix; 0];
    TB(end,[1:(end-2) end])=-TB(find(IB==[3+2*nvar+2*m]),[1:(end-2) end]);
    %         TB(end,[1:end])=-TB(find(IB==[3+2*nvar+2*m]),[1:end]);
    A(:,3:4)=[xnew -xnew ];
end

% TB(end,(end-1))=TB(1:(end-1),3+2*nvar+2*m)'*bi*xnew;
c=[c ;0];
cib=c(IB);
TB=[TB; c(1:end-1)'-cib'*bi*A 0 0];



j=1;
while (j<5000)
    % start iteration
    % step 2
    r=TB(end-1,1:(end-1));
    if (all(r>=-1e-7))
        %         sprintf('done')
        break
    end
    % step 3
    [~,t]=min(r);
    %     t=find(r<0,1);
    % step 4
    yy=TB(1:(end-2),t);
    if (all(yy<=0))
        sprintf('unboundness')
        break
    end
    % step 5
    b=TB(1:(end-2),end);
    k=b./yy;
    k(yy<=1e-7)=inf;
    [~,k]=min(k);
    
    % step 6
    ee=-TB(:,t)./TB(k,t);
    ee(k)=1/TB(k,t)-1;
    
    IB(k)=t;
    % ID=setdiff(1:(2*(1+nvar)+2*m),IB);
    % step 7
    %     TB=E * TB;
    TB=TB+ee*TB(k,:);
    
    
    %     TB(abs(TB)<1e-7)=0; % too time consuming
    j=j+1;
    if (j==4800); sprintf('not converge stage 1'); end;
end

TB=TB([1:(end-2) end],[1:(end-2) end]);

while (j<5000)
    % start iteration
    % step 2
    r=TB(end,1:(end-1));
    if (all(r>=-1e-7))
        %         sprintf('done')
        break
    end
    % step 3
    [~,t]=min(r);
    % step 4
    yy=TB(1:(end-1),t);
    %     if (all(yy<=0))
    %         sprintf('unboundness')
    %         break
    %     end
    % step 5
    b=TB(1:(end-1),end);
    k=b./yy;
    k(yy<=1e-7)=inf;
    [~,k]=min(k);
    
    % step 6
    ee=-TB(:,t)./TB(k,t);
    ee(k)=1/TB(k,t)-1;
    
    IB(k)=t;
    % ID=setdiff(1:(2*(1+nvar)+2*m),IB);
    % step 7
    %     TB=E * TB;
    TB=TB+ee*TB(k,:);
    
    
    %     TB(abs(TB)<1e-7)=0; % too time consuming
    j=j+1;
    if (j==4800); sprintf('not converge stage 2'); end;
end


b=TB(1:(end-1),end);
estimate=zeros((nvar+1),1);
for i=1:2*(nvar+1)
    if isempty(b(IB==i))
        continue
    end
    if (mod(i,2)==1)
        estimate((i+1)/2)=b(IB==i);
    else
        estimate(i/2)=-b(IB==i);
    end
end
estimate(2)=flag*estimate(2);
% A=A(:,1:(end-1));
% c=c(1:(end-1));
end


