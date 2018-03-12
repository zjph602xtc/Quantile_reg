function [ estimate,j, A,TB,IB,c,nvar ] = altermyqr(xnew,A,TB,IB,c,nvar,y)
% m=length(y);
% This function might not be right under some circumstances.

xnew(y<0)=-xnew(y<0);
b=y;
b(b<0)=-b(b<0);
A(:,3:4)=[xnew -xnew];
B=A(:,IB);
Bi=inv(B);
Ba=Bi*A;
Bb=Bi*b;
% TB=[inv(B) zeros(m,1); -c(IB)'/B 1] * [A b; c' 0];
TB=[Ba Bb;c'-c(IB)'*Ba c(IB)'*Bb];

% TB=sparse(TB);
% !!!!!! This code shows only one changed var is correct!!!
j=1;
while (any(TB(1:(end-1),end)<0))
    [~,k]=min(TB(1:(end-1),end));
    r=TB(k,1:(end-1));
    [~,t]=min(r);
    IB(k)=t;
    ee=-TB(:,t)./TB(k,t);
    ee(k)=1/TB(k,t)-1;
    TB=TB+ee*TB(k,:); 
    j=j+1;
end


while (1)
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
end


