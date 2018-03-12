addpath(genpath('E:\dataarr'));
listfile=ls('E:\dataarr\hou');
listfile=listfile(3:end,:);
rng(123)
sampleindex=datasample(1:23948,5,'Replace',false);
listfile=listfile(sampleindex,:);
cov=readtable('cov.txt','FileType','text','ReadVariableNames',1,'ReadRowNames',1);

i = 3
while 1
    dat=readtable(strtrim(listfile(i,:)),'FileType','text');
    if width(dat)~=1
        break
    end
    sprintf('Change index')
    sampleindex(i)=datasample(setdiff(1:23948,sampleindex),1);
end

%dat=dat(:,1:200);
y=dat{:,1};
dat=[dat(:,2:end) cov];
%%%%%%%%%%%%%%%%%%%%%%%
%1) change one x at each time
tau=0.5
[ A,TB,IB,c,nvar ] = mypreqr(dat{:,[1 end-39:end]},y,tau); % fit all z's + one x
TB=full(TB);

tic
iter=0;
for i=2:500 % only calculated 50 x's of 290, to save time
    [~,n,A,TB,IB] = altermyqr3(dat{:,i},A,TB,IB,c,nvar,y);
    iter=iter+n;
    if mod(i,50)==0
    sprintf('iter = %.f,  i = %.f',iter,i)
    end
end
toc

%2) fit all z's first, then add one x at each time
[ A,TB,IB,c,nvar ] = mypreqr(dat{:,[end-39:end]},y,tau); % fit all z's
iABI=inv(A(:,IB));
TB(abs(TB)<1e-7)=0;iABI(abs(iABI)<1e-7)=0;
TB=full(TB);

tic
iteradd=0;
est=zeros(42,500);
for i=1:500 % only calculated 50 x's of 290, to save time
    [est(:,i),n ] = addmyqr(dat{:,i},TB,IB,iABI,c,nvar,y);
    iteradd=iteradd+n;
    if mod(i,50)==0
    sprintf('iter = %.f,  i = %.f',iteradd,i)
    end
end
toc

%%%%%%%%%%%%%%%%%%
% 3) one by one 
tic
iteradd=0;
est=zeros(42,500);
for i=1:500 % only calculated 50 x's of 290, to save time
    [est(:,i),n ] = myqr3(dat{:,[4 end-39:end]},y,tau);
    iteradd=iteradd+n;
    if mod(i,50)==0
    sprintf('iter = %.f,  i = %.f',iteradd,i)
    end
end
toc


