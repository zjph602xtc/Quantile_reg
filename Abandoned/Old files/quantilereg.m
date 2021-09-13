function [beta tstats VCboot itrat PseudoR2 betaboot]=quantilereg(y,x,p)

% This Mfile estimates quantile regression based on weighted least squares.
 
% Inputs:
%  y, Dependent variable
%  x, matrix of independent variables

% Outputs:
%  beta, estimated Coefficients.
%  tstats, T- students of the coefficients.
%  VCboot, Variance-Covariance of coefficients by bootstrapping method.
%  itrat, number of iterations for convergence to roots.
%  PseudoR2, in quatile regression another definition of R2 is used namely 
%  PseudoR2.
%  betaboot, estimated coefficients by bootstrapping method.


% This code can be used for quantile regression estimation as whole,and LAD
% regression as special case of it, when one sets tau=0.5.
 
 
% Copyright(c) Shapour Mohammadi, University of Tehran, 2008
% shmohammadi@gmail.com
 
%Ref:
% 1-Birkes, D. and Y. Dodge(1993). Alternative Methods of Regression,
% John Wiley and Sons. 
% 2-Green,W. H. (2008). Econometric Analysis. Sixth Edition. International
% Student Edition.
% 3-LeSage, J. P.(1999),Applied Econometrics Using MATLAB, 
 
% Keywords: Least Absolute Deviation(LAD) Regression, Quantile Regression,
% Regression, Robust Estimation.
%__________________________________________________________________________


% tic
ry=length(y);
[rx cx]=size(x);
x=[ones(rx,1) x];
cx=cx+1;
%______________Finding first estimates by solving the system_______________

% Some lines of this section is based on a code written by 
% James P. Lesage in Applied Econometrics Using MATLAB(1999).PP. 73-4.
itrat=0;
xstar=x;
diff=1;
beta=ones(cx,1);
while itrat<1000 && diff>1e-6
    itrat=itrat+1;
    beta0=beta;
beta=inv(xstar'*x)*xstar'*y;
resid=y-x*beta;
resid(abs(resid)<.000001)=.000001;
resid(resid<0)=p*resid(resid<0);
resid(resid>0)=(1-p)*resid(resid>0);
resid=abs(resid);
z=[];
for i=1:cx 
    z0 = x(:,i)./resid;
    z=[z z0];
end

xstar=z;
beta1=beta;
diff=max(abs(beta1-beta0));

end

e=y-x*beta;

%_______estimating variances based on Green 2008(quantile regression)______

iqre=iqr(e);
if p==0.5
  h=0.9*std(e)/(ry^0.2);
else
  h=0.9*min(std(e),iqre/1.34)/(ry^0.2);
end
u=(e/h);
fhat0=(1/(ry*h))*(sum(exp(-u)./((1+exp(-u)).^2)));
D(ry,ry)=0;
DIAGON=diag(D);
DIAGON(e>0)=(p/fhat0)^2;
DIAGON(e<=0)=((1-p)/fhat0)^2;
D=diag(DIAGON);
VCQ=(x'*x)^(-1)*x'*D*x*(x'*x)^(-1);


%____________________Standarad errores and t-stats_________________________

tstats=beta./diag(VCQ).^0.5;
stderrors=diag(VCQ).^0.5;
PValues=2*(1-tcdf(abs(tstats),ry-cx));

%______________________________ Quasi R square_____________________________

ef=y-x*beta;
ef(ef<0)=(1-p)*ef(ef<0);
ef(ef>0)=p*ef(ef>0);
ef=abs(ef);

ered=y-quantile(y,p);
ered(ered<0)=(1-p)*ered(ered<0);
ered(ered>0)=p*ered(ered>0);
ered=abs(ered);

PseudoR2=1-sum(ef)/sum(ered);



%__________________Bootstrap standard deviation (Green 2008)_______________

betaboot=zeros(cx);
for ii=1:100
[bootm estar]=bootstrp(1,@mean,e);

ystar=x*beta+e(estar);

itratstar=0;
xstarstar=x;
diffstar=1;
betastar=ones(cx,1);
while itratstar<1000 && diffstar>1e-6
    itratstar=itratstar+1;
    betastar0=betastar;
betastar=inv(xstarstar'*x)*xstarstar'*ystar;

residstar=ystar-x*betastar;
residstar(abs(residstar)<.000001)=.000001;
residstar(residstar<0)=p*residstar(residstar<0);
residstar(residstar>0)=(1-p)*residstar(residstar>0);
residstar=abs(residstar);
zstar=[];
for i=1:cx 
    zstar0 = x(:,i)./residstar;
    zstar=[zstar zstar0];
end
xstarstar=zstar;
beta1star=betastar;
diffstar=max(abs(beta1star-betastar0));
end

betaboot=[betaboot + (betastar-beta)*(betastar-beta)'];
end
VCboot=(1/100)*betaboot;

tstatsboot=beta./diag(VCboot).^0.5;
stderrorsboot=diag(VCboot).^0.5;
PValuesboot=2*(1-tcdf(abs(tstatsboot),ry-cx));

%_______________________________Display Results____________________________

% disp(' ')
% disp('  Results of Quantile Regression       ')
% disp('___________________________________________________________________')
% disp([blanks(4) 'Coef.'  blanks(5) 'SE.Ker'  blanks(4) 't.Ker' blanks(5)...
%  'P.Ker' blanks(5) 'SE.Boot' blanks(2) 't.Boot'  blanks(4) 'P.Boot'])
% disp('___________________________________________________________________')
% disp(  [ beta     ,     stderrors   ,       tstats,    PValues,...
%                         stderrorsboot,  tstatsboot,    PValuesboot] )
% disp('___________________________________________________________________')
% disp(  [blanks(4)  'Pseudo R2' ] );disp(                 PseudoR2)
% disp('___________________________________________________________________')
% toc

%____________________________________END___________________________________


