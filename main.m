n=500;
ab_size_vec=[1000:500:5000];
alph=0.05;
k_max=7;
coverage_BC_vec=zeros(length(ab_size_vec),1);
coverage_mine_vec=zeros(length(ab_size_vec),1);
log_vol_mine_vec=zeros(length(ab_size_vec),1);
log_vol_BC_vec=zeros(length(ab_size_vec),1);


for ab_size_ind=1:length(ab_size_vec)
    ab_size=ab_size_vec(ab_size_ind);
    %p=[1:1:ab_size].^-1.01;
    p=ones(ab_size,1);
    p=p'/sum(p);


    CI=CI_for_non_frequent_symbols(n,k_max,alph*(1-n/(ab_size*(k_max+1))));
        
        
    num_of_exp=10^3;   
    log_vol_naive_exp=zeros(num_of_exp,1);
    log_vol_mine_exp=zeros(num_of_exp,1);
    ok_naive=0;
    ok_mine=0;
    
    for exp_ind=1:num_of_exp
        [ab_size exp_ind]
        r = mnrnd(n,p);
        
        [~,pci_naive] = binofit(r,n,alph/ab_size);
        log_vol_naive_exp(exp_ind)=sum(log(pci_naive(:,2)-pci_naive(:,1)));
            
        if sum((p>=pci_naive(:,1)).*(p<=pci_naive(:,2)))==ab_size
            ok_naive=ok_naive+1;
        end
       
        
        pci_mine=pci_naive;
        for ind=0:k_max
            pci_mine(r==ind,:)=repmat(CI(ind+1,:),length(pci_mine(r==ind,1)),1);
        end
        log_vol_mine_exp(exp_ind)=sum(log(pci_mine(:,2)-pci_mine(:,1)));
            
        if sum((p>=pci_mine(:,1)).*(p<=pci_mine(:,2)))==ab_size
            ok_mine=ok_mine+1;
        end
       
        
    end
    coverage_BC_vec(ab_size_ind,1)=ok_naive/num_of_exp;
    coverage_mine_vec(ab_size_ind,1)=ok_mine/num_of_exp;
    log_vol_BC_vec(ab_size_ind,1)=mean(log_vol_naive_exp);
    log_vol_mine_vec(ab_size_ind,1)=mean(log_vol_mine_exp);
        
        
end






BC=log_vol_BC_vec./ab_size_vec';
kmax7=log_vol_mine_vec./ab_size_vec';

figure(1)
hold on
plot(ab_size_vec,BC,'r','linewidth',2)
plot(ab_size_vec,kmax7,'b','linewidth',2)
set(gca,'FontSize', 13)
legend('BC','k max=7')
xlabel('ab size','fontsize',20,'interpreter','latex')
ylabel('log-volume','fontsize',20,'interpreter','latex')
box on
hold off






