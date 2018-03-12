function [ estimate,j ] = addmyqr(xnew,TB,IB,iABI,c,nvar,y)
m=length(y);

xnew(y<0)=-xnew(y<0);
BinvX=iABI*xnew;
NewTB=[BinvX -BinvX; -c(IB)'*BinvX c(IB)'*BinvX];
TB=[TB(:,1:(end-1)) NewTB TB(:,end)];
% !!!!!! This code shows only one added var is correct!!!

j=1; nnzz=[];
while (1)
    % start iteration
    % step 2
    r=TB(end,1:(end-1));
    if (all(r>=-1e-12))
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
    
%         nnzz = [nnzz nnz(TB)];
%     if (j~=1 && nnzz(j)-nnzz(j-1)>1800)
%         TB(abs(TB)<1e-10)=0;
%     end
    
    % step 6
    ee=-TB(:,t)./TB(k,t);
    ee(k)=1/TB(k,t)-1;
    TB=TB+ee*TB(k,:);
    
%     E=eye(m+1,m+1);
%     E(:,k)=[-TB(:,t)./TB(k,t)];
%     E(k,k)=1/TB(k,t);
%     TB=E * TB;
    

        
%     assignin('base','nnzz',nnz(TB));
%     assignin('base','nnzztrue',nnz(round(TB,8)));
%         evalin('base','nn=[nn [nnzz;nnzztrue]];');
    
    
%     TB=round(TB,6);
% TB=mtimesx(E,TB,'SPEED');
%    TB1=mtimesx(E,TB,'BLAS'); 
%       TB1=mtimesx(E,TB,'LOOPSOMP'); 
%                            TB1=E*TB;
   
    IB(k)=t;
    % ID=setdiff(1:(2*(1+nvar)+2*m),IB);
    % step 7

%     TB(abs(TB)<1e-7)=0; % too time consuming
    j=j+1;
end

clear estimate;
[~,nvarnew]=size(xnew); 
b=TB(1:(end-1),end);
estimate=zeros((nvar+nvarnew+1),1);
[~,II]=size(TB);
for i=[1:2*(nvar+1) II-2*nvarnew:II-1]
    if isempty(b(IB==i))
        continue
    end
    if (mod(i,2)==1)
        estimate((i+1)/2)=b(IB==i);
    else
        estimate(i/2)=-b(IB==i);
    end
end
estimate(estimate==0)=[];
end

