function invA = invA(Aib,IB,m,y,nvar,nalt,invAzero)
% x adds ones,  if y<0 x=-x

ys=y<0;

R=find(IB <= (1+nvar) | (IB > 1+nvar+m & IB <= 1+nvar+m+nalt));
% IBR=IB(R);
% IBRs = IBR>1+nvar+m;
% IBR(IBRs) = IBR(IRBs)-m-nalt;
IB(R)=[];

IBs = IB> 1+nvar+m;
R=[R' my_setdiff(1:m,R)];

IB(IBs)=IB(IBs)-m-nalt;
IB=IB-1-nvar;
ysIB=ys(IB);

negsign =[false(1+nvar,1); xor(IBs,ysIB)];
IB = [my_setdiff(1:m,IB) IB'];

%
%
% Aib=A(:,IB);
% temp = Aib(IB,R);

%
% x=[ones(m,1) x];
% x(y<0,:)=-x(y<0,:);

x12=Aib(IB(1:1+nvar),R(1:1+nvar));
x21=Aib(IB(2+nvar:end),R(1:1+nvar));

b112=inv(x12);
invAtmp=[[b112; -x21*b112] invAzero];
invAtmp(negsign,:)=-invAtmp(negsign,:);
[~,IBsinv]=sort(IB);
[~,Rsinv]=sort(R);
invA=invAtmp(Rsinv,IBsinv);


end
