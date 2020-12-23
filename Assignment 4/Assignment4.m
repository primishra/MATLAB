%Q1-A
clear all
clc
syms x;
f=(x^2);
q=int(f,x,0,5);
fprintf('the integral is:%f \n',double(q));


%%
%Q1-B
clear all
clc
tspan=[0,1];
x0=0;
[t,x]=ode45(@(t,x) t^2-x , tspan , x0);
fprintf('      t        numerical solution \n\n');
disp([t, x]);
 

%%
clear all
clc
%Q1-C
tspan=[0,1];
x0=0;
[t,x]=ode15s(@(t,x) t^2-x , tspan , x0);
fprintf('      t        numerical solution \n \n');
disp([t,x]);




%%
%Q2
clear all
clc
syms x;
% A.
f=sin(x^2);
p=int(f,x,0,1);
fprintf('the integral is:%f \n',double(p));


% B.
g=exp(-(x^2));
q=int(g,x,0,inf);
fprintf('the integral is:%f \n',double(q));



%%
%Q3
clear all
clc
syms x;
g=1/(1+x^2);
integral=int(g,x,0,6); 
res=double(integral);
% range of definite integral
a=0;
b=6;
n=input('enter no. of sub-interval:');
f=matlabFunction(g); % converts symbolic expression to funtion handle
h=abs(b-a)/n; % calculating the value of x

% Trapezoidal Rule
fprintf('\n trapeoidal rule\n');

sum1=f(a)+f(b); % sum of first and last term
for i=1:n-1
    sum1=sum1+2*f(a+(i*h)); % sum of middle terms
    
end
res1=sum1*(h/2); % INTEGRAL USING TRAPEZOIDAL RULE
fprintf(' approximate value of integral is :%f \n',res1);

err1= abs(res-res1); % difference between the actual & approximated value

fprintf('difference between the actual & approximated value:%f \n',err1);

% Simpson’s One third rule
fprintf("\n Simpson's One third rule \n");
sum3=f(a)+f(b);

for i=1:n-1
    if mod(i,2)==0
        sum3=sum3+2*f(a+(i*h)); 
    else
       sum3=sum3+4*f(a+(i*h)); 
    end
    
end
res3=sum3*(h/3);
fprintf(' approximate value of integral is :%f \n',res3);
err3=abs(res-res3);
fprintf('difference between the actual & approximated value:%f \n',err3);

% simpson's three eighth rule
fprintf(" \n simpson's three eighth rule \n ");
sum3=f(a)+f(b);
for i=1:n-1
    if mod(i,3)==0
        sum3=sum3+2*f(a+(i*h)); 
    else
       sum3=sum3+3*f(a+(i*h)); 
    end
    
end
res3=sum3*((3*h)/8);
fprintf(' approximate value of integral is :%f \n',res3);
err3=abs(res-res3);
fprintf('difference between the actual & approximated value:%f \n',err3);





%%
%Q4
clear all
clc
tspan=[0,2];
x0=0;
sol=ode45(@(t,x) t^2-x , tspan , x0);
p=linspace(0,2,10); 
%  x1 gives solution using ode45 at points contained in p
x1=deval(sol,p);
disp('Numerical solution:');
disp('t:');
disp(p);
disp('corrosponding numerical solution:');
disp(x1);


syms x(t);
ode=diff(x,t)==t^2-x;
cond=x(0)==0;
xsol(t)=dsolve(ode,cond);
temp=subs(xsol,p);
% x2 gives exact solution points contained in p
x2=double(temp);
disp('Exact solution:');
disp('t:');
disp(p);
disp('corrosponding  solution:');
disp(x2);

plot(p,x1,'+k',p,x2,'or');
title('comparision of numerical solution and exact solution ');

xlabel('t');
ylabel('solution');
legend('numerical solution','exact solution');


 

%%
% Q5
clear all
clc

alpha=1; beta=0.05; mu=0.02; eta=0.5;
%initial population of both  prey & predator is 10 
y0=[10;10];
% for prey: alpha*x(1)-beta*x(1)*x(2)
% for predator: mu*x(1)*x(2)-eta*x(2)
%  x(1) denotes the number of prey.
% x(2) refer to the number of predators.
f=@(t,x)[alpha*x(1)-beta*x(1)*x(2);mu*x(1)*x(2)-eta*x(2)];
[t,y]=ode45(f,[0 50],y0);
% [t,y] consists of a column of times  (t) and a matrix of populations (y)
% 1st column of y corresponds with x(1)
% and the 2nd column corresponds  with  x(2) 

y1=y(:,1); % prey
y2=y(:,2); % predator 
% solutions as x vs t and y vs t in a single graph.
figure
plot(t,y1,'g',t,y2,'r');
xlabel('time');
ylabel('population ');
legend('prey','predator');
grid on;
title('solutions as x vs t and y vs t in a single graph.');

%%
% Q5.C

alpha=1; beta=0.05;  eta=0.5;
%initial population of both  prey & predator is 10 
y0=[10;10];
% for prey: alpha*x(1)-beta*x(1)*x(2)
% for predator: mu*x(1)*x(2)-eta*x(2)
%  x(1) denotes the number of prey.
% x(2) refer to the number of predators.

i=2;

for mu=0.01:0.005:0.025
f=@(t,x)[alpha*x(1)-beta*x(1)*x(2);mu*x(1)*x(2)-eta*x(2)];

[t,y]=ode45(f,[0 20],y0);
% [t,y] consists of a column of times  (t) and a matrix of populations (y)
% 1st column of y corresponds with x(1)
% and the 2nd column corresponds  with  x(2) 

y1=y(:,1); % prey
y2=y(:,2); % predator 
% solutions as x vs t and y vs t in a single graph.
subplot(3,3,i);
i=i+2;
p=plot(t,y1,'g',t,y2,'r');

xlabel('time')
ylabel('population ');
suptitle('solutions as x vs t and y vs t in a single graph');
title([' mu= ' num2str(mu)]);

end

subplot(3,3,3);
plot(0,0,'g' ,0,0,'r');
axis off;
legend('prey','predator');
