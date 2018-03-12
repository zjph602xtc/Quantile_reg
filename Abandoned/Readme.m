%%%% Change 1 predictor at each regression %%%%
% addmyqr.m, altermyqr3.m are designed for this situation

% read sample data. There are 40 covariates that do not change.
load('sample.mat')
tau=0.8;
Z=arrayfun(@(i)sprintf('Z%d',i),1:40,'unif',0);
X=arrayfun(@(i)sprintf('X%d',i),1:290,'unif',0);
y=samplemat{:,'Y'};

%1) the general way to do quantile regression
tic
nn=0;
for i=1:50 % only calculated 50 x's of 290, to save time
    [~,n ] = myqr(samplemat{:,[Z X(i)]},y,tau);
    nn=nn+n';
end
toc

%2) fit all z's first, then add one x at each time
[ A,TB,IB,c,nvar ] = mypreqr(samplemat{:,Z},y,tau); % fit all z's
TB=full(TB);

tic
nn=0;
for i=1:50 % only calculated 50 x's of 290, to save time
    [~,n ] = addmyqr(samplemat{:,X(i)},A,TB,IB,c,nvar,y);
    nn=nn+n;
end
toc

%3) change one x at each time
startx='X51'; % choose a starting variable x
sortedx = sortx( samplemat(:,X),startx ); % I sorted x's here. We could use the location info to sort data in the real data set.
[ A,TB,IB,c,nvar ] = mypreqr(samplemat{:,[startx Z]},y,tau); % fit all z's + one x
TB=full(TB);

tic
nn=0;
for i=1:50 % only calculated 50 x's of 290, to save time
    [~,n,A,TB,IB] = altermyqr3(samplemat{:,sortedx(i)},A,TB,IB,c,nvar,y);
    nn=nn+n;
end
toc


%%%% Change m predictors at each regression %%%%
% multiaddmyqr.m, multialtermyqr.m are designed for this situation

% read sample data. There 10 z's.
n = 500;
p = 10;
x = normrnd(0,1,n,1);
z = normrnd(0,1,n,p);
bet0 = @(tau, alp00 , alp01, x, k){
        alp00 ./ p .* k + (alp01 + 0.5 .* x) .* log(tau);
    }
y = [];
for i = 1:n
    u = unifrnd(0,1,1,1);
    betas = [];
    for k = 1:p
        betas = [betas, cell2mat(bet0(u, -3.3,  0.86, x(i), k))];
    end
    y = [y; z(i,:)*betas'];
end
x = sort(x);
xx=[];
for j = (p + 2):(n - p - 2)
    xx = [xx 1 * (x <= x(j))];
end
tau=0.8;
% The model is y= x*beta_1 + z*beta_2 + x*z*beta_3.
% We change x and x*z at each regression.

%1)  the general way to do quantile regression
tic
nn=0;
for i=1:50
    inter = bsxfun(@times, z, xx(:,i));
    [~,n ] = myqr([z xx(:,i) inter],y,tau);
    nn=nn+n';
end
toc

%2) fit all z's first, then add x and x*z at each time
[ A,TB,IB,c,nvar ] = mypreqr(z,y,tau); % fit all z's
TB=full(TB);

tic
nn=0;
for i=1:50
    inter = bsxfun(@times, z, xx(:,i));
    [~,n ] = multiaddmyqr([xx(:,i) inter],A,TB,IB,c,nvar,y);
    nn=nn+n;
end
toc

%3) change x and x*z at each time
startx=1; % choose a starting variable x
inter = bsxfun(@times, z, xx(:,startx));
[ A,TB,IB,c,nvar ] = mypreqr([x(:,startx) inter z],y,tau); % fit all startx & startx*z & z 
TB=full(TB);

tic
nn=0;
for i=1:50
    inter = bsxfun(@times, z, xx(:,i));
    [~,n,A,TB,IB] = multialtermyqr([xx(:,i) inter],A,TB,IB,c,nvar,y);
    nn=nn+n;
end
toc
