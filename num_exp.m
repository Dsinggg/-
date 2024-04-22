% 绘制的图像: (1) X (2) p (3) q (4) pi (5) S
    
syms t A1(t) A2(t) B1(t) B2(t) D1(t) D2(t) E1(t) E2(t) F1(t) F2(t)

rng(12345)

alpha=20;   
beta=8;
c1=5;
c2=6;
k=0.1;
delta=0.05;
X10=0;
X20=50;
s1=0;
s2=0;
T=10;

[t_, y] = ode45(@(t,y) myODE(t,y,alpha,beta,delta,c1,c2,k), [T, 0], zeros(10,1));

A1=y(:,1);
B1=y(:,2);
D1=y(:,3);
E1=y(:,4);
F1=y(:,5);
A2=y(:,6);
B2=y(:,7);
D2=y(:,8);
E2=y(:,9);
F2=y(:,10);

[t__, x] = ode45(@(t,y) myODE2(t,y,A1,B1,D1,E1,F1,A2,B2,D2,E2,F2,t_,alpha,beta,delta,c1,c2,k), [0, T], [X10,X20]);

X1=x(:,1);
X2=x(:,2);

timeline=0:0.1:T;
dt=0.1;

dW1 = sqrt(dt) * randn(1, length(timeline));
W1 = cumsum(dW1);
dW2 = sqrt(dt) * randn(1, length(timeline));
W2 = cumsum(dW2);

X1__=zeros(length(timeline),1);
X2__=zeros(length(timeline),1);
p1=zeros(length(timeline),1);
p2=zeros(length(timeline),1);
q1=zeros(length(timeline),1);
q2=zeros(length(timeline),1);

index=1;

for t=timeline

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
    X1_=interp1(t__,X1,t)+s1*W1(index);
    X2_=interp1(t__,X2,t)+s2*W2(index);

    V1_1=A1_+2*D1_*X1_+F1_*X2_;
    V1_2=B1_+2*E1_*X2_+F1_*X1_;

    V2_1=B2_+2*E2_*X1_+F2_*X2_;
    V2_2=A2_+2*D2_*X2_+F2_*X1_;
    
    p1_i=(2*c1+c2)/3+(3*alpha+k*(2*X1_+X2_))/(3*beta)-(2*V1_1-2*V1_2+V2_2-V2_1)/3;
    p2_i=(2*c2+c1)/3+(3*alpha+k*(2*X2_+X1_))/(3*beta)-(2*V2_2-2*V2_1+V1_1-V1_2)/3;

    X1__(index)=X1_;
    X2__(index)=X2_;

    p1(index)=p1_i;
    p2(index)=p2_i;

    q1(index)=alpha-beta*(p1_i-p2_i)+k*X1_;
    q2(index)=alpha-beta*(p2_i-p1_i)+k*X2_;
   
    index=index+1;

end

plot(timeline,X1__,'-')
hold on
plot(timeline,X2__,'-')
figure

plot(timeline,p1)
hold on
plot(timeline,p2)
color_p1 = get(gca, 'ColorOrder');
color_p2 = get(gca, 'ColorOrder');
line([0 T], [c1 c1], 'Color', color_p1(1, :), 'LineStyle', '--');
line([0 T], [c2 c2], 'Color', color_p2(2, :), 'LineStyle', '--');
figure

plot(timeline,q1)
hold on
plot(timeline,q2)
figure

pi1=(p1-c1).*q1;
pi2=(p2-c2).*q2;
plot(timeline,pi1)
hold on
plot(timeline,pi2)
figure

S1 = cumtrapz(timeline, pi1);
S2 = cumtrapz(timeline, pi2);
plot(timeline,S1)
hold on
plot(timeline,S2)
figure

plot(timeline,p2-p1,'-')
