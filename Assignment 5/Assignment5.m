%{
Question 1. Evaluate the working of the following functions with suitable examples (a) plot, (b) subplot (c) plot3 (d) bar (e) meshgrid, surfplot Write a program on the following and display the output:
}%

%Create Line Plot- Question-1a
x = 0:pi/100:2*pi;
y = sin(x);
plot(x,y)
%Plot Multiple Lines
x = linspace(-2*pi,2*pi);
y1 = sin(x);
y2 = cos(x);
figure
plot(x,y1,x,y2)
%Create Line Plot From Matrix
Y = magic(4)
figure
plot(Y)
%Specify Line Style
x = 0:pi/100:2*pi;
y1 = sin(x);
y2 = sin(x-0.25);
y3 = sin(x-0.5);

figure
plot(x,y1,x,y2,'--',x,y3,':')
%Specify Line Style, Color, and Marker
x = linspace(0,10);
y = sin(x);
plot(x,y,'-o','MarkerIndices',1:5:length(y))
%Specify Line Width, Marker Size, and Marker Color
x = -pi:pi/10:pi;
y = tan(sin(x)) - sin(tan(x));

figure
plot(x,y,'--gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
%Plot Circle
r = 2;
xc = 4;
yc = 3;

theta = linspace(0,2*pi);
x = r*cos(theta) + xc;
y = r*sin(theta) + yc;
plot(x,y)
axis equal

%Question 1-b

%Upper and Lower Subplots
subplot(2,1,1);
x = linspace(0,10);
y1 = sin(x);
plot(x,y1)

subplot(2,1,2); 
y2 = sin(5*x);
plot(x,y2)
%Quadrant of Subplots

subplot(2,2,1)
x = linspace(0,10);
y1 = sin(x);
plot(x,y1)
title('Subplot 1: sin(x)')

subplot(2,2,2)
y2 = sin(2*x);
plot(x,y2)
title('Subplot 2: sin(2x)')

subplot(2,2,3)
y3 = sin(4*x);
plot(x,y3)
title('Subplot 3: sin(4x)')

subplot(2,2,4)
y4 = sin(8*x);
plot(x,y4)
title('Subplot 4: sin(8x)')
%Subplots with Different Sizes
subplot(2,2,1);
x = linspace(-3.8,3.8);
y_cos = cos(x);
plot(x,y_cos);
title('Subplot 1: Cosine')

subplot(2,2,2);
y_poly = 1 - x.^2./2 + x.^4./24;
plot(x,y_poly,'g');
title('Subplot 2: Polynomial')

subplot(2,2,[3,4]);
plot(x,y_cos,'b',x,y_poly,'g');
title('Subplot 3 and 4: Both')

%Question - 1-c
%Plot 3-D Helix
t = 0:pi/50:10*pi;
st = sin(t);
ct = cos(t);
plot3(st,ct,t)
%Plot Multiple Lines
t = 0:pi/500:pi;
xt1 = sin(t).*cos(10*t);
yt1 = sin(t).*sin(10*t);
zt1 = cos(t);

xt2 = sin(t).*cos(12*t);
yt2 = sin(t).*sin(12*t);
zt2 = cos(t);
figure
plot3(xt1,yt1,zt1,xt2,yt2,zt2)
%Specify Equally-Spaced Tick Units and Axis Labels
t = 0:pi/500:40*pi;
xt = (3 + cos(sqrt(32)*t)).*cos(t);
yt = sin(sqrt(32) * t);
zt = (3 + cos(sqrt(32)*t)).*sin(t);
figure
plot3(xt,yt,zt)
axis equal
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')

%Question-1-D
%Create Bar Graph
y = [75 91 105 123.5 131 150 179 203 226 249 281.5];
bar(y)
%Create Bar Graph with Categorical Data
c = categorical({'apples','pears','oranges'});
prices = [1.23 0.99 2.3];
bar(c,prices)
%Display Stacked Bars
y = [2 2 3; 2 5 6; 2 8 9; 2 11 12];
figure
bar(y,'stacked')

%Question-1-e
%2-D Grid
x = 1:3;
y = 1:5;
[X,Y] = meshgrid(x,y)
X.^2 + Y.^2
%Plot Surface
x = -2:0.25:2;
y = x;
[X,Y] = meshgrid(x);
F = X.*exp(-X.^2-Y.^2);
surf(X,Y,F)
surf(x,y,x.*exp(-x.^2-(y').^2))
%%





%{
Question 2. . Using Newton–Raphson’s method, find the real root of the equation f(x) = x ? e ?x correct to four decimal places. Compare the result with Bisection and Regula Falsi method.
}%

f=@(x) (x-exp(-x));
display(f);
fderv=@(x) (1+exp(-x));
for i=-1000:1000
    j=i+1;
    if f(i)*f(j)<0
        a2=i;
        b2=j;
    end
end
lhl=a2;
rhl=b2;
fprintf('The roots lie in the inteval [%d,%d]\n',lhl,rhl);
fprintf('***NEWTON RAPHSON***\n');
x0=((a2+b2)/2);
fprintf('Taking initial roots to be %f\n',x0);
count1=1;
while (abs(f(x0))>0.0001)
    fprintf('iteration %d\n',count1);
    y=x0-(f(x0)/fderv(x0));
    fprintf('root %d: %f\n',count1,y);
    count1=count1+1;
    x0=y;
end
fprintf('the root of the equation is : %f\n',x0);
fprintf('Newton Raphson Method took %d iterations\n',count1-1);
fprintf('***BISECTION METHOD***\n');
a1=lhl;
b1=rhl;
c1=(a1+b1)/2;
count2=1;
while abs(f(c1))>0.0001
    fprintf('iteration %d\n',count2);
    c1=(a1+b1)/2;
    if f(c1)*f(a1)<0
        b1=c1;
    elseif f(c1)*f(b1)<0
        a1=c1;
    end
    fprintf('The roots lie in the inteval [%d,%d]\n',a1,b1);
    count2=count2+1;
end
fprintf('the root of the equation is : %f\n',c1);
fprintf('Bisection Method took %d iterations\n',count2-1);
fprintf('***REGULA FALSI METHOD***\n');
a2=lhl;
b2=rhl;
c2=(a2*f(b2)-b2*f(a2))/(f(b2)-f(a2));
count3=1;
while abs(f(c2))>0.0001
    fprintf('iteration %d\n',count3);
    c2=(a2*f(b2)-b2*f(a2))/(f(b2)-f(a2));
    if f(c2)*f(a2)<0
        b2=c2;
    elseif f(c2)*f(b2)<0
        a2=c2;
    end
    fprintf('The roots lie in the inteval [%d,%d]\n',a2,b2);
    count3=count3+1;
end
fprintf('the root of the equation is : %f\n',c2);
fprintf('Regula Falsi Method took %d iterations\n',count3-1);
%%
Output


% 
% The roots lie in the inteval [0,1]
% ***NEWTON RAPHSON***
% Taking initial roots to be 0.500000
% iteration 1
% root 1: 0.566311
% iteration 2
% root 2: 0.567143
% the root of the equation is : 0.567143
% Newton Raphson Method took 2 iterations
% ***BISECTION METHOD***
% iteration 1
% The roots lie in the inteval [5.000000e-01,1]
% iteration 2
% The roots lie in the inteval [5.000000e-01,7.500000e-01]
% iteration 3
% The roots lie in the inteval [5.000000e-01,6.250000e-01]
% iteration 4
% The roots lie in the inteval [5.625000e-01,6.250000e-01]
% iteration 5
% The roots lie in the inteval [5.625000e-01,5.937500e-01]
% iteration 6
% The roots lie in the inteval [5.625000e-01,5.781250e-01]
% iteration 7
% The roots lie in the inteval [5.625000e-01,5.703125e-01]
% iteration 8
% The roots lie in the inteval [5.664063e-01,5.703125e-01]
% iteration 9
% The roots lie in the inteval [5.664063e-01,5.683594e-01]
% iteration 10
% The roots lie in the inteval [5.664063e-01,5.673828e-01]
% iteration 11
% The roots lie in the inteval [5.668945e-01,5.673828e-01]
% iteration 12
% The roots lie in the inteval [5.671387e-01,5.673828e-01]
% the root of the equation is : 0.567139
% Bisection Method took 12 iterations
% ***REGULA FALSI METHOD***
% iteration 1
% The roots lie in the inteval [0,6.126998e-01]
% iteration 2
% The roots lie in the inteval [0,5.721814e-01]
% iteration 3
% The roots lie in the inteval [0,5.677032e-01]
% iteration 4
% The roots lie in the inteval [0,5.672056e-01]
% the root of the equation is : 0.567206
% Regula Falsi Method took 4 iterations 


%%
%{
Question 3. Find the double root of the equation f(x) = x 3?7x 2+16x?12 using Newton–Raphson’s method correct to three decimal places. Compare the obtained result with Generalized Newton’s method (or modified Newton–Raphson method).
}%

syms x;
f=x^3-7*x^2+16*x-12;; %Enter the Function here
g=diff(f); %The Derivative of the Function
n=input('Enter the number of decimal places:');
epsilon = 5*10^-(n+1)
x0 = input('Enter the intial approximation:');
for i=1:100
     f0=vpa(subs(f,x,x0)); %Calculating the value of function at x0
     f0_der=vpa(subs(g,x,x0)); %Calculating the value of function derivative at x0
  y=x0-f0/f0_der; % The Formula
err=abs(y-x0);
if err<epsilon %checking the amount of error at each iteration
break
end
x0=y;
end
y = y - rem(y,10^-n); %Displaying upto required decimal places
fprintf("This is from NewtonsRaphson Method")
fprintf('The Root is : %f \n',y);
fprintf('No. of Iterations : %d\n',i);

fprintf("######################################\n")
syms x;
f=x^3-7*x^2+16*x-12; %Enter the Function here
g=diff(f); %The Derivative of the Function
n=input('Enter the number of decimal places:');
epsilon = 5*10^-(n+1);
u=f/g;
v=diff(u);
x0 = input('Enter the intial approximation:');
for i=1:100
     f0=vpa(subs(u,x,x0)); %Calculating the value of function at x0
     f0_der=vpa(subs(v,x,x0)); %Calculating the value of function derivative at x0
     
  y=x0-f0/f0_der; % The Formula
err=abs(y-x0);
if err<epsilon %checking the amount of error at each iteration
break
end
x0=y;
end
y = y - rem(y,10^-n); %Displaying upto required decimal places
fprintf("This from Modified Newton's raphson")
fprintf('The Root is : %f \n',y);
fprintf('No. of Iterations : %d\n',i);

%{
	Output:

Ques3
Enter the number of decimal places:3

epsilon =

   5.0000e-04

Enter the intial approximation:1
This is from NewtonsRaphson MethodThe Root is : 1.999000 
No. of Iterations : 12
######################################
Enter the number of decimal places:3
Enter the intial approximation:1
This from Modified Newton's raphsonThe Root is : 2.000000 
No. of Iterations : 4
%}



%{
Question 4. Evaluate the positive root of the equation lying near 0.5 of the equation cos(x)?xex = 0 correct to four decimal places using Newton’s method. Compare the obtained result with Regula–Falsi method.
}%

% f-It is the function to be solved 
% a=First end; b=Second end 
% tol-It is the tolerance limit(ans correct upto given decimal point )
%  kmax -Maximum number of iteration
% Regula falsi
%#######################
%%

f=@(x) (cos(x)-x*exp(x));

a=0;b=2; 
kmax=10000;
double tol;
tol=input("Enter the value of tolearnace ");
a_o=a;b_o=b;
ya=f(a);
yb=f(b);
if sign(ya)==sign(yb)
    error("function have same sign at a and b");
end
%%%disp("      Step      a         b        m         ya       yb         ym     ")
for k=1:kmax
    m=b -(f(b)*(b-a))/(f(b)-f(a));

    ym=f(m);
    iter=k;
    %out=[iter,a,b,m,f(a),f(b),f(m)];
    %disp( out);
    if(abs(ym)<tol)
        disp("Root  has reached tolerance limit ");
        break;
    end
    if(sign(ym)==sign(ya))
        a=m;
    else 
        b=m;
    end
    if(iter>kmax)
        disp("Maximum number of iteration crossed,tolrance limit not reached")
    end

end

fprintf("Regual falsi\n");
fprintf("No of iterartions=%d\n",iter);
%fprintf("a=%d,b=%d\n",a_o,b_o);
disp("Answer")
disp("m")
disp(m)
disp("value at m")
disp(num2str(ym,'%.10f\n'))
disp("##############")
disp("Newton's Raphson")
syms x;
f=(cos(x)-x*exp(x)); %Enter the Function here
g=diff(f); %The Derivative of the Function
n=input('Enter the number of decimal places:');
epsilon = 10^-(n+1);
x0 = input('Enter the intial approximation:');
for i=1:100
     f0=vpa(subs(f,x,x0)); %Calculating the value of function at x0
     f0_der=vpa(subs(g,x,x0)); %Calculating the value of function derivative at x0
  y=x0-f0/f0_der; % The Formula
err=abs(y-x0);
if err<epsilon %checking the amount of error at each iteration
break
end
x0=y;
end
y = y - rem(y,10^-n); %Displaying upto required decimal places
fprintf('The Root is : %f \n',y);
fprintf('No. of Iterations : %d\n',i);

%{

Output:
Ques4
Enter the value of tolearnace 0.0001
Root  has reached tolerance limit 
Regual falsi
No of iterartions=29
Answer
m
    0.5177

value at m
0.0000784489
##############
Newton's Raphson
Enter the number of decimal places:4
Enter the intial approximation:1
The Root is : 0.517700 
No. of Iterations : 5%}




