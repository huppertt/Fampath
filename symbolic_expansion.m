
syms c0 c1 c2 hbo hbr W lambda rho a b r
assume(c0,'real')
assume(c1,'real')
assume(c2,'real')
assume(hbo,'real')
assume(hbr,'real')
assume(W,'real')
assume(lambda,'real')
assume(rho,'real')
assume(a,'real')
assume(b,'real')
assume(r,'real')

assume(c0,'positive')
assume(c1,'positive')
assume(c2,'positive')
assume(hbo,'positive')
assume(hbr,'positive')
assume(W,'positive')
assume(lambda,'positive')
assume(rho,'positive')
assume(a,'positive')
assume(b,'positive')
assume(r,'positive')

mua = c0.*W+c1.*hbo+c2.*hbr;
mus = a.*lambda.^(-b);
mue = (3.*mua.*(mua + mus)).^(1/2); % effective attenuation coefficient
D = 1./(3.*(mua+mus));
z0 = 1./(mua+mus);
zb = ((1+0.493)./(1-0.493)).*2.*D; %0.493 is Reff at n=1.4
x0 = z0.^2 + rho.^2;
x1 = (z0+2.*zb).^2 + rho.^2; 
refl = simplify((1/4*pi).*((z0.*(mue+1./sqrt(x0)).*exp(-mue.*sqrt(x0))./x0) +...
    ((z0 + 2*zb).*(mue + 1./sqrt(x1)).*exp(-mue.*sqrt(x1))./x1)));

s = simplify(r/(refl*100000)); 

ds_dhbo=diff(s,hbo);
ds_dhbr=diff(s,hbr);
ds_dW=diff(s,W);
ds_da=diff(s,a);
ds_db=diff(s,b);

ds_dx = [ds_dhbo ds_dhbr ds_dW ds_da ds_db]

ds_dx=simplify(ds_dx)


