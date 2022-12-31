function [a,m_opt] = lower_edge(n,k,alph)
    res=10^-5;
    m=[res:res:k-1-res];
    t=m./(m+n-k);
    %obj=(alph./(nchoosek(n,k).*t.^m.*(1-t).^(n-k))).^(1./(k-1-m));
    obj=(alph./(binopdf(k,n,t).*t.^(m-k))).^(1./(k-1-m));
    [a,ind]=max(obj);
    m_opt=m(ind);

end