syms alpha beta p1 p2 k X1 X2 s1 s2 delta c1 c2

syms V1 V1_1 V1_2 V1_11 V1_22 W1 A1 B1 D1 E1 F1

syms V2 V2_1 V2_2 V2_11 V2_22 W2 A2 B2 D2 E2 F2

% solve(diff((alpha-beta*(p1-p2)+k*X1)*(p1-c1)+V1_1*(alpha-beta*(p1-p2)+(k-delta)*X1)+V1_2*(alpha-beta*(p2-p1)+(k-delta)*X2)+V1_11*s1^2/2+V1_22*s2^2/2,p1)==0,p1)

solve([diff((alpha-beta*(p1-p2)+k*X1)*(p1-c1)+V1_1*(alpha-beta*(p1-p2)+(k-delta)*X1)+V1_2*(alpha-beta*(p2-p1)+(k-delta)*X2)+V1_11*s1^2/2+V1_22*s2^2/2,p1)==0;...
       diff((alpha-beta*(p2-p1)+k*X2)*(p2-c2)+V2_1*(alpha-beta*(p1-p2)+(k-delta)*X1)+V2_2*(alpha-beta*(p2-p1)+(k-delta)*X2)+V2_11*s1^2/2+V2_22*s2^2/2,p2)==0],[p1,p2])

V1=W1+A1*X1+B1*X2+D1*X1^2+E1*X2^2+F1*X1*X2;

V1_1=A1+2*D1*X1+F1*X2;
V1_2=B1+2*E1*X2+F1*X1;
V1_11=2*D1;
V1_22=2*E1;

V2=W2+A2*X2+B2*X1+D2*X2^2+E2*X1^2+F2*X1*X2;

V2_2=A2+2*D2*X2+F2*X1;
V2_1=B2+2*E2*X1+F2*X2;
V2_22=2*D2;
V2_11=2*E2;

p1 = (3*alpha - 2*V1_1*beta + 2*V1_2*beta + V2_1*beta - V2_2*beta + 2*X1*k + X2*k + 2*beta*c1 + beta*c2)/(3*beta);
p2 = (3*alpha - V1_1*beta + V1_2*beta + 2*V2_1*beta - 2*V2_2*beta + X1*k + 2*X2*k + beta*c1 + 2*beta*c2)/(3*beta);

RHS1=(alpha-beta*(p1-p2)+k*X1)*(p1-c1)+V1_1*(alpha-beta*(p1-p2)+(k-delta)*X1)+V1_2*(alpha-beta*(p2-p1)+(k-delta)*X2)+V1_11*s1^2/2+V1_22*s2^2/2;

[coef1,var]=coeffs(RHS1,[X1,X2]);

coef1_simp=simplify(coef1)

RHS2=(alpha-beta*(p2-p1)+k*X2)*(p2-c2)+V2_1*(alpha-beta*(p1-p2)+(k-delta)*X1)+V2_2*(alpha-beta*(p2-p1)+(k-delta)*X2)+V2_11*s1^2/2+V2_22*s2^2/2;

[coef2,var]=coeffs(RHS2,[X1,X2]);

coef2_simp=simplify(coef2)

var