function CI=CI_for_non_frequent_symbols(n,k_max,alph)
    k_min=0;
    k_vec=[k_min:1:k_max];
    CI_marginals=zeros(length(k_vec),2);
    r_vec=zeros(length(k_vec),1);
    m_vec=zeros(length(k_vec),1);
    for k_ind=1:length(k_vec)
        k=k_vec(k_ind);
        
        if k>=2
            [b,r_opt] = upper_edge(n,k,alph/2);
            [a,m_opt] = lower_edge(n,k,alph/2);
        else
            [b,r_opt] = upper_edge(n,k,alph);
            a=0;
            m_opt=-3;
        end
        CI_marginals(k_ind,1)=a;
        CI_marginals(k_ind,2)=b;
        r_vec(k_ind)=r_opt;
        m_vec(k_ind)=m_opt;
    end
    
    
    c_res=10^-3;
    c_vec=[1:c_res:10];
    t_res=10^-6;
    t=[t_res:t_res:1-t_res];
    ind=1;
    done=0;
    while ~done
        
        c=c_vec(ind);
        current_CI_marginals=[mean(CI_marginals,2)-c*0.5*(CI_marginals(:,2)-CI_marginals(:,1)),mean(CI_marginals,2)+c*0.5*(CI_marginals(:,2)-CI_marginals(:,1))];
        current_CI_marginals(current_CI_marginals<0)=0;
        obj=zeros(1,length(t));
        for k_ind=1:length(k_vec)
            k=k_vec(k_ind);
            %current_obj=nchoosek(n,k)*(current_CI_marginals(k_ind,2).^-r_vec(k_ind).*t.^(r_vec(k_ind)+k-1).*(1-t).^(n-k)+current_CI_marginals(k_ind,1).^(-m_vec(k_ind)+k-1).*t.^m_vec(k_ind).*(1-t).^(n-k));
            current_obj=binopdf(k,n,t).*(current_CI_marginals(k_ind,2).^-r_vec(k_ind).*t.^(r_vec(k_ind)-1)+current_CI_marginals(k_ind,1).^(-m_vec(k_ind)+k-1).*t.^(m_vec(k_ind)-k));

            obj=obj+current_obj;
        end 
        max_obj=max(obj);
        if max_obj<alph
            done=1;
        else
            ind=ind+1;
            if ind>length(c_vec)
                bug=1;
            end
        end
        [k_max ind max_obj]
    
    end
    c_opt=c_vec(ind);
    CI=current_CI_marginals;






end