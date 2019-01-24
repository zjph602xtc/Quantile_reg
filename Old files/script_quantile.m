n=1000;
beta=[0.2 0.3];
gamma=[0 1];
tau=[1:4]/5;

coe=[];
for i=1:150
    x=random('norm',0,1,[1,n]);
    z=random('norm',4,1,[1,n]);
    y=1+beta(1)*x+beta(2)*z+(1+gamma(1)*x+gamma(2)*z).*random('chi2',2,[1,n]);
    dat=table(y',x',z','VariableNames',{'y' 'x' 'z'});
    coe=[coe mycompositerq(dat,tau,'x','y',{'z'})];
    i
end

mean(coe,2,'omitnan')
norminv(tau)+1

chi2inv(tau,2)+1


