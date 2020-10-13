function dxdt=diff_eq(t,x,x_sym,f)

f=subs(f,x_sym,sym(x));
f=sym(f);

dxdt=eval(f);


end