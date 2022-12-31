function [b,r_opt] = upper_edge(n,k,alph)
    res=10^-3;
    r=[1:res:20];
    t=(r+k-1)./(r+n-1);
    obj=(binopdf(k,n,t).*t.^(r-1)/alph).^(1./r);
    [b,ind]=min(obj);
    r_opt=r(ind);

end