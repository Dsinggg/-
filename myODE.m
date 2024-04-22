function dydt = myODE(t,y,alpha,beta,delta,c1,c2,k)

A1=y(1);
B1=y(2);
D1=y(3);
E1=y(4);
F1=y(5);
A2=y(6);
B2=y(7);
D2=y(8);
E2=y(9);
F2=y(10);

dydt=zeros(10,1);

dydt(1)= (4*D1*alpha)/3 - A1*delta + (4*E2*alpha)/3 + (4*F1*alpha)/3 - (2*F2*alpha)/3 + (4*A1*k)/9 - (4*A2*k)/9 + (5*B1*k)/9 + (4*B2*k)/9 - (4*c1*k)/9 + (4*c2*k)/9 + (4*A1*D1*beta)/9 - (4*A2*D1*beta)/9 - (4*B1*D1*beta)/9 + (4*A1*E2*beta)/9 + (4*B2*D1*beta)/9 - (4*A2*E2*beta)/9 - (2*A1*F1*beta)/9 - (2*A1*F2*beta)/9 + (2*A2*F1*beta)/9 - (4*B1*E2*beta)/9 + (2*A2*F2*beta)/9 + (4*B2*E2*beta)/9 + (2*B1*F1*beta)/9 + (2*B1*F2*beta)/9 - (2*B2*F1*beta)/9 - (2*B2*F2*beta)/9 - (4*D1*beta*c1)/9 + (4*D1*beta*c2)/9 - (4*E2*beta*c1)/9 + (4*E2*beta*c2)/9 + (2*F1*beta*c1)/9 - (2*F1*beta*c2)/9 + (2*F2*beta*c1)/9 - (2*F2*beta*c2)/9 + (4*alpha*k)/(3*beta);
dydt(2)= (8*E1*alpha)/3 - B1*delta - (4*D2*alpha)/3 + (2*F1*alpha)/3 + (2*F2*alpha)/3 + (2*A1*k)/9 - (2*A2*k)/9 + (7*B1*k)/9 + (2*B2*k)/9 - (2*c1*k)/9 + (2*c2*k)/9 - (4*A1*D2*beta)/9 + (4*A2*D2*beta)/9 - (4*A1*E1*beta)/9 + (4*A2*E1*beta)/9 + (4*B1*D2*beta)/9 - (4*B2*D2*beta)/9 + (2*A1*F1*beta)/9 + (4*B1*E1*beta)/9 + (2*A1*F2*beta)/9 - (2*A2*F1*beta)/9 - (4*B2*E1*beta)/9 - (2*A2*F2*beta)/9 - (2*B1*F1*beta)/9 - (2*B1*F2*beta)/9 + (2*B2*F1*beta)/9 + (2*B2*F2*beta)/9 + (4*D2*beta*c1)/9 - (4*D2*beta*c2)/9 + (4*E1*beta*c1)/9 - (4*E1*beta*c2)/9 - (2*F1*beta*c1)/9 + (2*F1*beta*c2)/9 - (2*F2*beta*c1)/9 + (2*F2*beta*c2)/9 + (2*alpha*k)/(3*beta);
dydt(3)= (4*beta*D1^2)/9 + (8*beta*D1*E2)/9 - (4*beta*D1*F1)/9 - (4*beta*D1*F2)/9 + (8*D1*k)/9 - 2*delta*D1 + (4*beta*E2^2)/9 - (4*beta*E2*F1)/9 - (4*beta*E2*F2)/9 + (8*E2*k)/9 + (beta*F1^2)/9 + (2*beta*F1*F2)/9 + (5*F1*k)/9 + (beta*F2^2)/9 - (4*F2*k)/9 + (4*k^2)/(9*beta);
dydt(4)= (4*beta*D2^2)/9 + (8*beta*D2*E1)/9 - (4*beta*D2*F1)/9 - (4*beta*D2*F2)/9 - (4*D2*k)/9 + (4*beta*E1^2)/9 - (4*beta*E1*F1)/9 - (4*beta*E1*F2)/9 + (14*E1*k)/9 - 2*delta*E1 + (beta*F1^2)/9 + (2*beta*F1*F2)/9 + (2*F1*k)/9 + (beta*F2^2)/9 + (2*F2*k)/9 + k^2/(9*beta);
dydt(5)= (4*D1*k)/9 - 2*F1*delta - (8*D2*k)/9 + (10*E1*k)/9 + (4*E2*k)/9 + (11*F1*k)/9 + (2*F2*k)/9 - (2*F1^2*beta)/9 - (2*F2^2*beta)/9 + (4*k^2)/(9*beta) - (8*D1*D2*beta)/9 - (8*D1*E1*beta)/9 - (8*D2*E2*beta)/9 + (4*D1*F1*beta)/9 + (4*D1*F2*beta)/9 + (4*D2*F1*beta)/9 - (8*E1*E2*beta)/9 + (4*D2*F2*beta)/9 + (4*E1*F1*beta)/9 + (4*E1*F2*beta)/9 + (4*E2*F1*beta)/9 + (4*E2*F2*beta)/9 - (4*F1*F2*beta)/9;
% (beta*A1^2)/9 - (2*beta*A1*A2)/9 - (2*beta*A1*B1)/9 + (2*beta*A1*B2)/9 + (2*A1*alpha)/3 - (2*beta*A1*c1)/9 + (2*beta*A1*c2)/9 + (beta*A2^2)/9 + (2*beta*A2*B1)/9 - (2*beta*A2*B2)/9 - (2*A2*alpha)/3 + (2*beta*A2*c1)/9 - (2*beta*A2*c2)/9 + (beta*B1^2)/9 - (2*beta*B1*B2)/9 + (4*B1*alpha)/3 + (2*beta*B1*c1)/9 - (2*beta*B1*c2)/9 + (beta*B2^2)/9 + (2*B2*alpha)/3 - (2*beta*B2*c1)/9 + (2*beta*B2*c2)/9 + alpha^2/beta - (2*alpha*c1)/3 + (2*alpha*c2)/3 + (beta*c1^2)/9 - (2*beta*c1*c2)/9 + (beta*c2^2)/9 + D1*s1^2 + E1*s2^2 

dydt(7)= (8*E2*alpha)/3 - B2*delta - (4*D1*alpha)/3 + (2*F1*alpha)/3 + (2*F2*alpha)/3 - (2*A1*k)/9 + (2*A2*k)/9 + (2*B1*k)/9 + (7*B2*k)/9 + (2*c1*k)/9 - (2*c2*k)/9 + (4*A1*D1*beta)/9 - (4*A2*D1*beta)/9 - (4*B1*D1*beta)/9 + (4*A1*E2*beta)/9 + (4*B2*D1*beta)/9 - (4*A2*E2*beta)/9 - (2*A1*F1*beta)/9 - (2*A1*F2*beta)/9 + (2*A2*F1*beta)/9 - (4*B1*E2*beta)/9 + (2*A2*F2*beta)/9 + (4*B2*E2*beta)/9 + (2*B1*F1*beta)/9 + (2*B1*F2*beta)/9 - (2*B2*F1*beta)/9 - (2*B2*F2*beta)/9 - (4*D1*beta*c1)/9 + (4*D1*beta*c2)/9 - (4*E2*beta*c1)/9 + (4*E2*beta*c2)/9 + (2*F1*beta*c1)/9 - (2*F1*beta*c2)/9 + (2*F2*beta*c1)/9 - (2*F2*beta*c2)/9 + (2*alpha*k)/(3*beta);
dydt(6)= (4*D2*alpha)/3 - A2*delta + (4*E1*alpha)/3 - (2*F1*alpha)/3 + (4*F2*alpha)/3 - (4*A1*k)/9 + (4*A2*k)/9 + (4*B1*k)/9 + (5*B2*k)/9 + (4*c1*k)/9 - (4*c2*k)/9 - (4*A1*D2*beta)/9 + (4*A2*D2*beta)/9 - (4*A1*E1*beta)/9 + (4*A2*E1*beta)/9 + (4*B1*D2*beta)/9 - (4*B2*D2*beta)/9 + (2*A1*F1*beta)/9 + (4*B1*E1*beta)/9 + (2*A1*F2*beta)/9 - (2*A2*F1*beta)/9 - (4*B2*E1*beta)/9 - (2*A2*F2*beta)/9 - (2*B1*F1*beta)/9 - (2*B1*F2*beta)/9 + (2*B2*F1*beta)/9 + (2*B2*F2*beta)/9 + (4*D2*beta*c1)/9 - (4*D2*beta*c2)/9 + (4*E1*beta*c1)/9 - (4*E1*beta*c2)/9 - (2*F1*beta*c1)/9 + (2*F1*beta*c2)/9 - (2*F2*beta*c1)/9 + (2*F2*beta*c2)/9 + (4*alpha*k)/(3*beta);
dydt(9)= (4*beta*D1^2)/9 + (8*beta*D1*E2)/9 - (4*beta*D1*F1)/9 - (4*beta*D1*F2)/9 - (4*D1*k)/9 + (4*beta*E2^2)/9 - (4*beta*E2*F1)/9 - (4*beta*E2*F2)/9 + (14*E2*k)/9 - 2*delta*E2 + (beta*F1^2)/9 + (2*beta*F1*F2)/9 + (2*F1*k)/9 + (beta*F2^2)/9 + (2*F2*k)/9 + k^2/(9*beta);
dydt(8)= (4*beta*D2^2)/9 + (8*beta*D2*E1)/9 - (4*beta*D2*F1)/9 - (4*beta*D2*F2)/9 + (8*D2*k)/9 - 2*delta*D2 + (4*beta*E1^2)/9 - (4*beta*E1*F1)/9 - (4*beta*E1*F2)/9 + (8*E1*k)/9 + (beta*F1^2)/9 + (2*beta*F1*F2)/9 - (4*F1*k)/9 + (beta*F2^2)/9 + (5*F2*k)/9 + (4*k^2)/(9*beta);
dydt(10)= (4*D2*k)/9 - (8*D1*k)/9 - 2*F2*delta + (4*E1*k)/9 + (10*E2*k)/9 + (2*F1*k)/9 + (11*F2*k)/9 - (2*F1^2*beta)/9 - (2*F2^2*beta)/9 + (4*k^2)/(9*beta) - (8*D1*D2*beta)/9 - (8*D1*E1*beta)/9 - (8*D2*E2*beta)/9 + (4*D1*F1*beta)/9 + (4*D1*F2*beta)/9 + (4*D2*F1*beta)/9 - (8*E1*E2*beta)/9 + (4*D2*F2*beta)/9 + (4*E1*F1*beta)/9 + (4*E1*F2*beta)/9 + (4*E2*F1*beta)/9 + (4*E2*F2*beta)/9 - (4*F1*F2*beta)/9;
% (beta*A1^2)/9 - (2*beta*A1*A2)/9 - (2*beta*A1*B1)/9 + (2*beta*A1*B2)/9 - (2*A1*alpha)/3 - (2*beta*A1*c1)/9 + (2*beta*A1*c2)/9 + (beta*A2^2)/9 + (2*beta*A2*B1)/9 - (2*beta*A2*B2)/9 + (2*A2*alpha)/3 + (2*beta*A2*c1)/9 - (2*beta*A2*c2)/9 + (beta*B1^2)/9 - (2*beta*B1*B2)/9 + (2*B1*alpha)/3 + (2*beta*B1*c1)/9 - (2*beta*B1*c2)/9 + (beta*B2^2)/9 + (4*B2*alpha)/3 - (2*beta*B2*c1)/9 + (2*beta*B2*c2)/9 + alpha^2/beta + (2*alpha*c1)/3 - (2*alpha*c2)/3 + (beta*c1^2)/9 - (2*beta*c1*c2)/9 + (beta*c2^2)/9 + E2*s1^2 + D2*s2^2

end