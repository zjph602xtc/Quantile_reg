function [pvalue,allslope] = simulation(n,gamma,tau,xsit,esit)
pvalue=[];
ltau=length(tau);
beta=evalin('base','beta');
allslope=[];
for i=1:500
    fprintf('\n ******************\n This is OUTER %.0f.\n ******************\n',i)
    switch xsit
        case 'bino'
            x=random('bino',2,0.3,n,1);
        case 'norm'
            x=random('norm',0,1,n,1);
        otherwise
            warning('Wrong xsit')
    end
    z=random('norm',4,1,n,1);
    switch esit
        case 'norm'
            err=random('norm',0,1,n,1);
        case 'chi2'
            err=random('chi2',2,n,1);
        case 't'
            err=random('t',2,n,1);
        case 'cauchy'
            err=random('t',1,n,1);
        otherwise
            warning('Wrong esit')
    end
    
    y=1+beta(1)*x+beta(2)*z+(1+gamma(1)*x+gamma(2)*z).*err;
    dat=table(y,x,z,'VariableNames',{'y' 'x' 'z'});
    
    if i<=30
        bootslope=zeros(1,200);
        parfor k=1:200
            if k==1
                datsample=dat;
            else
                datsample=datasample(dat,n);
            end
            x=datsample{:,'x'};
            y=datsample{:,'y'};
            z=datsample{:,'z'};
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
        est=bootslope(1);
        allslope=[allslope;bootslope];
        pvalue=[pvalue 2*normcdf(-abs(est)./std(bootslope))];
    else
        res=mycompositerq(dat,tau,'x','y',{'z'});
        resi=arrayfun(@(i) y - res(1)*x- res((i-1)*2+3)*z-res((i-1)*2+2),1:ltau,'unif',0);
        ranks=cellfun(@(u, tau) tau-1*(u<0),resi,num2cell(tau),'unif',0);
        ranks=cell2mat(ranks);
        Sn_i_tau=bsxfun(@times,ranks,x);
        %   Sn_i_tau=ranks.*repmat(x,1,ltau);
        slope=arrayfun(@(n)sum((tau-mean(tau)).*(Sn_i_tau(n,:)-mean(Sn_i_tau(n,:))))/sum((tau-mean(tau)).^2),1:n);
        pvalue=[pvalue 2*normcdf(-abs(sum(slope))./mean(std(allslope,[],2)))];
    end
    % ture value
    %     resi_t=arrayfun(@(taui )y-1-beta(1)*x-beta(2)*z-(1+gamma(1)*x+gamma(2)*z).*icdf('bino',taui,1,0.3,tau,'unif',0);
    %     ranks_t=cellfun(@(u, tau) tau-1*(u<0), resi_t,num2cell(tau),'unif',0);
    %     Sn_i_tau_t=cell2mat(ranks_t).*repmat(x,1,ltau);
    %     slope_t=arrayfun(@(n)sum((tau-mean(tau)).*(Sn_i_tau_t(n,:)-mean(Sn_i_tau_t(n,:))))/sum((tau-mean(tau)).^2),1:n);
    %     est_t=sum(slope_t); % estimate sum(beta_1_t)
    %     vv_t=reshape(datasample(slope_t,1200*n),n,1200);
    %     slopeSd_t=std(sum(vv_t));
    %     %%% save
    %     slopeSd_t_all=[slopeSd_t_all slopeSd_t];
    %     estimate_t=[estimate_t est_t];
    %     slope_t_all=[slope_t_all;slope_t];
    %     pvalue_t=[pvalue_t 2*normcdf(-abs(est_t)./slopeSd_t)];
      
end

end

