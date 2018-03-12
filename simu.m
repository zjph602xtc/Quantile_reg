beta=[0.2, 0.3];
tau5=[1:5]/6;
tau9=[1:9]/10;
tau15=[1:15]/16;
n=600;
gamma=[0 1];

[p,slopeall]=simulation(n,gamma,tau5,'bino','norm');
xlswrite('simu-8_4.xlsx',p',1,'BP7')

[p,slopeall]=simulation(n,gamma,tau9,'bino','norm');
xlswrite('simu-8_4.xlsx',p',1,'CH7')

[p,slopeall]=simulation(n,gamma,tau5,'norm','norm');
xlswrite('simu-8_4.xlsx',p',1,'BY7')

[p,slopeall]=simulation(n,gamma,tau9,'norm','norm');
xlswrite('simu-8_4.xlsx',p',1,'CQ7')

[p,slopeall]=simulation(n,gamma,tau5,'bino','chi2');
xlswrite('simu-8_4.xlsx',p',1,'BQ7')

[p,slopeall]=simulation(n,gamma,tau9,'bino','chi2');
xlswrite('simu-8_4.xlsx',p',1,'CI7')

[p,slopeall]=simulation(n,gamma,tau5,'norm','chi2');
xlswrite('simu-8_4.xlsx',p',1,'BZ7')

[p,slopeall]=simulation(n,gamma,tau9,'norm','chi2');
xlswrite('simu-8_4.xlsx',p',1,'CR7')

[p,slopeall]=simulation(n,gamma,tau5,'bino','cauchy');
xlswrite('simu-8_4.xlsx',p',1,'BS7')

[p,slopeall]=simulation(n,gamma,tau9,'bino','cauchy');
xlswrite('simu-8_4.xlsx',p',1,'CK7')

[p,slopeall]=simulation(n,gamma,tau5,'norm','cauchy');
xlswrite('simu-8_4.xlsx',p',1,'CB7')

[p,slopeall]=simulation(n,gamma,tau9,'norm','cauchy');
xlswrite('simu-8_4.xlsx',p',1,'CT7')

%%
[p,slopeall]=simulation(n,gamma,tau15,'bino','norm');
xlswrite('simu-8_4.xlsx',p',1,'CZ7')

[p,slopeall]=simulation(n,gamma,tau15,'bino','t');
xlswrite('simu-8_4.xlsx',p',1,'DB7')

[p,slopeall]=simulation(n,gamma,tau15,'norm','norm');
xlswrite('simu-8_4.xlsx',p',1,'DI7')

[p,slopeall]=simulation(n,gamma,tau15,'norm','t');
xlswrite('simu-8_4.xlsx',p',1,'DK7')

%%
mean(p<0.05)
plot(0:0.01:1,arrayfun(@(x)mean(pvalue<x),0:0.01:1),'r')
line(0:0.01:1,arrayfun(@(x)mean(p<x),0:0.01:1),'color','b')
line([0 1],[0 1])

%%
beta=[0.2, 0.3];
tau=[1:5]/6;
n=600;
gamma=[0 1];

pvalue=[];
slope_all=[];
estimate=[];
ltau=length(tau);

for i=1:250
    x=random('bino',2,0.3,n,1);
    err=random('norm',0,1,n,1);
    z=random('norm',4,1,n,1);
    y=1+beta(1)*x+beta(2)*z+(1+gamma(1)*x+gamma(2)*z).*err;
    dat=table(y,x,z,'VariableNames',{'y' 'x' 'z'});
    
    res=mycompositerq(dat,tau,'x','y',{'z'});
    resi=arrayfun(@(i) y - res(1)*x- res((i-1)*2+3)*z-res((i-1)*2+2),1:ltau,'unif',0);
    ranks=cellfun(@(u, tau) tau-1*(u<0),resi,num2cell(tau),'unif',0);
    Sn_i_tau=bsxfun(@times,cell2mat(ranks),x);
    slope=arrayfun(@(n)sum((tau-mean(tau)).*(Sn_i_tau(n,:)-mean(Sn_i_tau(n,:))))/sum((tau-mean(tau)).^2),1:n);
    est=sum(slope); % estimate sum(beta_1)
    
    %     vv=reshape(datasample(slope,1200*n),n,1200);
    %     slopeSd=std(sum(vv));
    %%% save
    estimate=[estimate est];
    %     slope_all=[slope_all;slope];
    %     pvalue=[pvalue 2*normcdf(-abs(est)./slopeSd)];
   %% %%boot
   bootslope=zeros(1,200);
   for k=1:200
       if k==1
           datsample=dat;
       else
           datsample=datasample(dat,n);
           x=datsample{:,'x'};
           y=datsample{:,'y'};
           z=datsample{:,'z'};
       end
       res=mycompositerq(datsample,tau,'x','y',{'z'});
       resi=arrayfun(@(i) y - res(1)*x- res((i-1)*2+3)*z-res((i-1)*2+2),1:ltau,'unif',0);
       ranks=cellfun(@(u, tau) tau-1*(u<0),resi,num2cell(tau),'unif',0);
       ranks=cell2mat(ranks);
       Sn_i_tau=bsxfun(@times,ranks,x);
       %   Sn_i_tau=ranks.*repmat(x,1,ltau);
       slope=arrayfun(@(n)sum((tau-mean(tau)).*(Sn_i_tau(n,:)-mean(Sn_i_tau(n,:))))/sum((tau-mean(tau)).^2),1:n);
       bootslope(k)=sum(slope);
       if mod(k,20)==0
           fprintf('This is inner %.0f.\n',k)
       end
   end
    
    i
end

%% check for overall
plot(0:0.01:1,arrayfun(@(t)mean(pvalue_t<t),0:0.01:1))
line([0 1],[0 1],'lines','--')
line(0:0.01:1,arrayfun(@(t)mean(pvalue<t),0:0.01:1),'color','k')
line(0:0.01:1,arrayfun(@(t)mean(pvaluewild<t),0:0.01:1),'color','g')
%% check for slope estimation
plot(sum(slope_t_all,2),'.')
hold on
plot(sum(slope_all,2),'k.')

mean(sum(slope_t_all,2))
mean(sum(slope_all,2))

%% check for std
plot(slopeSd_t_all,'.')
hold on
plot(slopeSd_all,'k.')

mean(slopeSd_t_all)
mean(slopeSd_all)

%% correct
std(estimate)
pright=normcdf(-abs(estimate)/std(estimate))*2;

plot(0:0.01:1,arrayfun(@(t)mean(pvalue_t<t),0:0.01:1))
line([0 1],[0 1],'lines','--')
line(0:0.01:1,arrayfun(@(t)mean(pright<t),0:0.01:1),'color','k')


%% one sample
pvalue=[];pvalue_t=[];
slope_all=[];slope_t_all=[];
estimate=[];estimate_t=[];
slopeSd_all=[];slopeSd_t_all=[];
ltau=length(tau);
x=random('bino',2,0.3,n,1);
err=random('norm',0,1,n,1);
z=random('norm',4,1,n,1);
y=1+beta(1)*x+beta(2)*z+(1+gamma(1)*x+gamma(2)*z).*err;
dat1=table(y,x,z,'VariableNames',{'y' 'x' 'z'});
for i=1:250
    dat=datasample(dat1,n);
    res=mycompositerq(dat,tau,'x','y',{'z'});
    resi=arrayfun(@(i) y - res(1)*x- res((i-1)*2+3)*z-res((i-1)*2+2),1:ltau,'unif',0);
    ranks=cellfun(@(u, tau) tau-1*(u<0),resi,num2cell(tau),'unif',0);
    Sn_i_tau=cell2mat(ranks).*repmat(x,1,ltau);
    slope=arrayfun(@(n)sum((tau-mean(tau)).*(Sn_i_tau(n,:)-mean(Sn_i_tau(n,:))))/sum((tau-mean(tau)).^2),1:n);
    est=sum(slope); % estimate sum(beta_1)
    vv=reshape(datasample(slope,1200*n),n,1200);
    slopeSd=std(sum(vv));
    %%% save
    slopeSd_all=[slopeSd_all slopeSd];
    estimate=[estimate est];
    slope_all=[slope_all;slope];
    pvalue=[pvalue 2*normcdf(-abs(est)./slopeSd)];
    
    % ture value
    resi_t=arrayfun(@(taui )y-1-beta(1)*x-beta(2)*z-(1+gamma(1)*x+gamma(2)*z).*icdf('norm',taui,0,1),tau,'unif',0);
    ranks_t=cellfun(@(u, tau) tau-1*(u<0), resi_t,num2cell(tau),'unif',0);
    Sn_i_tau_t=cell2mat(ranks_t).*repmat(x,1,ltau);
    slope_t=arrayfun(@(n)sum((tau-mean(tau)).*(Sn_i_tau_t(n,:)-mean(Sn_i_tau_t(n,:))))/sum((tau-mean(tau)).^2),1:n);
    est_t=sum(slope_t); % estimate sum(beta_1_t)
    vv_t=reshape(datasample(slope_t,1200*n),n,1200);
    slopeSd_t=std(sum(vv_t));
    %%% save
    slopeSd_t_all=[slopeSd_t_all slopeSd_t];
    estimate_t=[estimate_t est_t];
    slope_t_all=[slope_t_all;slope_t];
    pvalue_t=[pvalue_t 2*normcdf(-abs(est_t)./slopeSd_t)];
    
    i
end

%% 09/26 newly added
tau9=[1:9]/10;
gamma=[1 1];
[p,slopeall]=simulation(n,gamma,tau9,'bino','norm');
mean(p<0.05)

[p,slopeall]=simulation(n,gamma,tau9,'bino','chi2');
mean(p<0.05)

[p,slopeall]=simulation(n,gamma,tau9,'bino','cauchy');
mean(p<0.05)
%
[p,slopeall]=simulation(n,gamma,tau19,'bino','norm');
mean(p<0.05)

[p,slopeall]=simulation(n,gamma,tau19,'bino','chi2');
mean(p<0.05)

[p,slopeall]=simulation(n,gamma,tau19,'bino','cauchy');
mean(p<0.05)

%%
[p,slopeall]=simulation(n,gamma,tau5,'norm','norm');
mean(p<0.05)

[p,slopeall]=simulation(n,gamma,tau5,'norm','chi2');
mean(p<0.05)

[p,slopeall]=simulation(n,gamma,tau5,'norm','cauchy');
mean(p<0.05)

[p,slopeall]=simulation(n,gamma,tau9,'norm','norm');
mean(p<0.05)

[p,slopeall]=simulation(n,gamma,tau9,'norm','chi2');
mean(p<0.05)

[p,slopeall]=simulation(n,gamma,tau9,'norm','cauchy');
mean(p<0.05)
%
[p,slopeall]=simulation(n,gamma,tau19,'norm','norm');
mean(p<0.05)

[p,slopeall]=simulation(n,gamma,tau19,'norm','chi2');
mean(p<0.05)

[p,slopeall]=simulation(n,gamma,tau19,'norm','cauchy');
mean(p<0.05)