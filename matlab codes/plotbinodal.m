x1=0.5;
N1=4;
N2=200;
z
function output=d2Spinodal(x1,N1,N2,zdW,T)
k=0.001985875;
chi=zdW/(k*T);
output=1./(N1*x1)+1./(N2*(1-x1))-2*chi;
output=output*(k*T);
end