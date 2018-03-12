% test whether the quantreg is right
clear all
n=1000;
beta=[0.2 0.3];
gamma=[0 1];
x=random('norm',0,1,n,1);
z=random('norm',4,1,n,1);

y=1+beta(1)*x+beta(2)*z+(1+gamma(1)*x+gamma(2)*z).*normrnd(0,1,n,1);
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('sample.mat')
tau=0.8;

Z=arrayfun(@(i)sprintf('Z%d',i),1:40,'unif',0);
% X=arrayfun(@(i)sprintf('X%d',i),1:3394,'unif',0);
X=arrayfun(@(i)sprintf('X%d',i),1:290,'unif',0);
y=samplemat{:,'Y'};

tic
nn=0;
for i=1:50
    [~,n ] = myqr(samplemat{:,[Z X(i)]},y,tau);
    nn=nn+n';
end
toc
% nn=45946  43.9 s  N=50  tau=0.8

%%%%%%%% compare
% ????
[ A,TB,IB,c,nvar ] = mypreqr(samplemat{:,Z},y,tau);
TB=full(TB);
tic
nn=0;estall=[];
for i=1:50
    [estimate,n ] = addmyqr(samplemat{:,X(i)},A,TB,IB,c,nvar,y);
    nn=nn+n;
        estall=[estall estimate];
end
toc
estall=[estall([1 end 2:end-1 ],:)];
% nn=1264  1.7 s  N=50  tau=0.8
%%%%%%%%%%%
% ??? ???
[ A,TB,IB,c,nvar ] = mypreqr(samplemat{:,['X60' Z]},y,tau);

tic
nn=0;
for i=1:50
    [~,n ] = altermyqr(samplemat{:,X(i)},A,TB,IB,c,nvar,samplemat{:,'Y'});
    nn=nn+n;
end
toc
% nn=1749  2.4 s  N=50  tau=0.8
%%%%%%%%%%%%%%%%%%%
% ?? ???
startx='X1';
sortedx = sortx( samplemat(:,X),startx );
[ A,~,IB,c,nvar ] = mypreqr(samplemat{:,[startx Z]},y,tau);
profile on
tic
nn=0;nall=[];
estall_1=[];
for i=1:355
%     [ estimate,n, A,~,IB,~,~] = altermyqr(samplemat{:,sortedx(i)},A,'',IB,c,nvar,y);
       [ estimate,n, A,~,IB,~,~] = altermyqr(samplemat{:,X(i)},A,'',IB,c,nvar,y);
%      [ estimate,n, A,~,IB,~,~] = cccaltermyqr_mex(samplemat{:,sortedx(i)},A,IB,c,nvar,y);
    nn=nn+n;
    estall_1=[estall_1 estimate];
    nall=[nall n];
end
toc
profile viewer
% nn=714  0.72 s  N=50  tau=0.8

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%

startx='X51';
sortedx = sortx( samplemat(:,X),startx );
[ A,TB,IB,c,nvar ] = mypreqr(samplemat{:,[startx Z]},y,tau);
TB=full(TB);
% profile on
tic
nn=0;nall=[];
estmy=[];
for i=1:20
    [ estimate,n, A,TB,IB ] = altermyqr3(samplemat{:,sortedx(i)},A,TB,IB,c,nvar,y);
    if any(IB==641); break; end;
    nn=nn+n;
    nall=[nall n];
    estmy=[estmy estimate];
end
toc
% profile viewer
% nn=579  0.70 s  N=50  tau=0.8
