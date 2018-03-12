function [ estimate,j ] = multiaddmyqr(xnew,A,TB,IB,c,nvar,y)
m=length(y);
[~,nvarnew]=size(xnew); 
xnew(y<0,:)=-xnew(y<0,:);

aix = zeros(m, nvarnew*2);
aix(:,1:2:end) = xnew;
aix(:,2:2:end) = -xnew;


BinvX=A(:,IB)\aix;

NewTB=[BinvX; -c(IB)'*BinvX];
TB=[TB(:,1:(end-1)) NewTB TB(:,end)];


j=1;
while (1)
    % start iteration
    % step 2
    r=TB(end,1:(end-1));
    if (all(r>=-1e-7))
%         sprintf('done')
        break
    end
    % step 3
    % t=find(r==min(r(ID)),1);
    [~,t]=min(r);
    % step 4
    yy=TB(1:(end-1),t);
    if (all(yy<=0))
        sprintf('unboundness')
        break
    end
    % step 5
    b=TB(1:(end-1),end);
    k=b./yy;
    k(yy<=1e-7)=inf;
    [~,k]=min(k);
    % k=find(k==min(k(k>=0)),1);
    
    % step 6
    ee=-TB(:,t)./TB(k,t);
    ee(k)=1/TB(k,t)-1;
    
    IB(k)=t;
    % ID=setdiff(1:(2*(1+nvar)+2*m),IB);
    % step 7
    TB=TB+ee*TB(k,:);
%     TB(abs(TB)<1e-7)=0; % too time consuming
    j=j+1;
end

clear estimate;

b=TB(1:(end-1),end);
estimate=[];
estimate1=[];
[~,II]=size(TB);
betaindex=[1:2:2*(nvar+1) II-2*nvarnew:2:II-1];
for i = 1:nvar+nvarnew+1
    if any(IB==betaindex(i))
        estimate = [estimate; b(IB==betaindex(i))];
        estimate1 = [estimate1; 0];
    end
    if any(IB==betaindex(i)+1)
        estimate = [estimate; 0];
        estimate1 = [estimate1; b(IB==betaindex(i)+1)];
    end
end

estimate=estimate-estimate1;
end

