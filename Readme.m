% The resutlt 'est' is the coefficients for intercept fisrt, then all x's.

% read sample data. There are 18 covariates that do not change.
load('sample.mat')
tau=0.8;

%1) the general way to do quantile regression (interior points)
est = [];
pvalue = [];
tic
for i=1:500 % only calculated 500 x's of 10000, to save time
    [beta, p]=qr_standard([z x(:,i)], y, tau,'test','wald','method','interior');
    est = [est beta];
    pvalue = [pvalue p];
end
toc

%2) the general way to do quantile regression (standard linear programming)
est = [];
pvalue = [];
tic
for i=1:500 % only calculated 500 x's of 10000, to save time
    [beta, p] = qr_standard([z x(:,i)], y, tau,'test','kernel');
    est = [est beta];
    pvalue = [pvalue p];
end
toc

%3) fit all z's first, then add one x at each time
tic
[est, pvalue] = qr_add(z,x,y,tau,1,'test','kernel');  % for all 10000 x's
toc


%4) change one x at each time
tic
[est, pvalue] = qr_alter(z,x(:,1:500),y,tau,1);  % for 500 x's
toc


%%%% Change 2 predictors at each regression %%%%
%5) fit all z's first, then add two x's at each time
tic
[est, pvalue] = qr_add(z,x,y,tau,2,'test','kernel');  % for all 5000 x's pairs
toc


%6) change two x's at each time
tic
[est, pvalue] = qr_alter(z,x(:,1:500),y,tau,2);  % for 250 x's pairs
toc