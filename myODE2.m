function dydt = myODE2(t,y,A1,B1,D1,E1,F1,A2,B2,D2,E2,F2,t_,alpha,beta,delta,c1,c2,k)

X1=y(1);
X2=y(2);

A1_=interp1(t_,A1,t);
B1_=interp1(t_,B1,t);
D1_=interp1(t_,D1,t);
E1_=interp1(t_,E1,t);
F1_=interp1(t_,F1,t);
A2_=interp1(t_,A2,t);
B2_=interp1(t_,B2,t);
D2_=interp1(t_,D2,t);
E2_=interp1(t_,E2,t);
F2_=interp1(t_,F2,t);

P=A1_-A2_-B1_+B2_+2*D1_*X1-2*D2_*X2-2*E1_*X2+2*E2_*X1-F1_*X1+F1_*X2-F2_*X1+F2_*X2;

dydt=zeros(2,1);

dydt(1)=alpha-beta*(c1-c2)/3+k*(2*X1+X2)/3+beta*P/3-delta*X1;
dydt(2)=alpha-beta*(c2-c1)/3+k*(2*X2+X1)/3-beta*P/3-delta*X2;

end