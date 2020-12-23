%Q1.a
fn= @round;
x=pi;
y=5;
y=feval(fn, x, y)


%Q1.b
p = [4 3 1];
x = [2 4 6];
y = polyval(p,x)


%Q1.c
p = [2 3 -6];
r = roots(p)
poly(r)


%Q1.d
syms a b c x
eqn = 3*x^2 + 2*x + 6 == 0
S = solve(eqn,x)
Sv =vpasolve(eqn,x)


%Q1.e
p= 3*x^4 + 4*x^2 + 6;
coeffs(p)


%Q1.f
syms x
M = [x x^3; x^2 x^4];
f(x) = M


%Q1.g
fun = @sin;
x0 = 1;
x = fzero(fun,x0)





%%
%Q2
p=[6 -25 31 0 -31 25 6];
r= roots(p)





%%

%Q3.
syms x;

% Input Section
y = input('Enter non-linear equations: ');
a = input('Enter first guess: ');
b = input('Enter second guess: ');
fprintf('Tolerable error');
e = 10^(-4)

a1=a;
b1=b;

a2=a;
b2=b;

n=0;
n1=0;

% Finding Functional Value
fa = eval(subs(y,x,a));
fb = eval(subs(y,x,b));

% Implementing Bisection Method
fprintf('\n Criteria using value of the function: ');
if fa*fb > 0 
    disp('Given initial values do not bracket the root.');
else
    c = (a+b)/2;
    fc = eval(subs(y,x,c));
    fprintf('\n\n\tn\t\ta\t\t\tb\t\t\tc\t\t\tf(c)\n');
    while abs(fc)>e
        fprintf('%f\t%f\t%f\t%f\t%f\n',n,a,b,c,fc);
        if fa*fc< 0
            b =c;
            n=n+1;
        else
            a =c;
            n=n+1;
        end
        c = (a+b)/2;
        fc = eval(subs(y,x,c));
    end
    fprintf('\nRoot is: %f\n', c);
end

% Finding Functional Value
fa1 = eval(subs(y,x,a1));
fb1 = eval(subs(y,x,b1));


fprintf('\n Criteria using difference of roots: ');
if fa1*fb1 > 0 
    disp('Given initial values do not bracket the root.');
else
    c1 = (a1+b1)/2;
    fc1 = eval(subs(y,x,c1));
    fprintf('\n\n\tn1\t\ta1\t\t\tb1\t\t\tc1\t\t\tf(c1)\n');
    while abs(b1-a1)>e
        fprintf('%f\t%f\t%f\t%f\t%f\n',n1,a1,b1,c1,fc1);
        if fa1*fc1< 0
            b1 =c1;
            n1=n1+1;
        else
            a1 =c1;
            n1=n1+1;
        end
        c1 = (a1+b1)/2;
        fc1 = eval(subs(y,x,c1));
    end
    fprintf('\nRoot is: %f\n', c1);
end

fprintf('\nThe root of the equation calculated using inbuilt func: ');
p=vpasolve(y,x,[a2 b2]);
p1=double(p)

err= abs(((c-p1)/p1) *100);
fprintf('\nRelative error using Function value as termination criteria: %f ', err);

err1= abs(((c1- p1)/p1) *100);
fprintf('\nRelative error using difference of roots: %f ', err1);





%%
%Q4

% Setting x as symbolic variable
syms x;

% Input Section
y = input('Enter non-linear equations: ');
a = input('Enter first guess: ');
b = input('Enter second guess: ');

% Finding Functional Value
fa = eval(subs(y,x,a));
fb = eval(subs(y,x,b));

a1=a;
b1=b;

a2=a;
b2=b

n=0;
n1=0;

fprintf('\nImplementing Bisection Method');
if fa*fb > 0 
    disp('Given initial values do not bracket the root.');
else
    c = (a+b)/2;
    fc = eval(subs(y,x,c));
    fprintf('\n\n\tn\t\ta\t\t\tb\t\t\tc\t\t\tf(c)\n');
    for j=0:10
        fprintf('%f\t%f\t%f\t%f\t%f\n',n,a,b,c,fc);
        if fa*fc< 0
            b =c;
            n=n+1;
        else
            a =c;
            n=n+1;
        end
        c = (a+b)/2;
        fc = eval(subs(y,x,c));
    end
    fprintf('\nRoot(bisection) is: %f\n', c);
end



% Finding Functional Value
fa1 = eval(subs(y,x,a1));
fb1 = eval(subs(y,x,b1));

fprintf('\nImplementing Regula Falsi Method');
if fa1*fb1 > 0 
    disp('Given initial values do not bracket the root.');
else
    c1 = a1 - (a1-b1) * fa1/(fa1-fb1);
    fc1 = eval(subs(y,x,c1));
    fprintf('\n\n\tn1\t\ta1\t\t\tb1\t\t\tc1\t\t\tf(c1)\n');
    for i=0:10
        fprintf('%f\t%f\t%f\t%f\t%f\n',n1,a1,b1,c1,fc1);
        if fa1*fc1< 0
            b1 =c1;
            fb1 = eval(subs(y,x,b1));
        else
            a1=c1;
            fa1 = eval(subs(y,x,a1));
        end
        n1=n1+1;
        c1 = a1 - (a1-b1) * fa1/(fa1-fb1);
        fc1 = eval(subs(y,x,c1));
    end
    fprintf('\nRoot (Regula Falsi) is: %f\n', c1);
end

fprintf('\nThe root of the equation calculated using inbuilt func: ');
p=vpasolve(y,x,[a2 b2]);
p1=double(p)

err= abs(((c-p1)/p1) *100);
fprintf('\nRelative error using Bisection method: %f ', err);

err1= abs(((c1- p1)/p1) *100);
fprintf('\nRelative error using Regula Falsi method: %f ', err1);







%%
%Q5
syms x;

% Input Section
y = input('Enter non-linear equations: ');
a=0;
b=0;
fprintf('Tolerable error');
e = 10^(-3)


%Finding suitable interval
for k=-100:100
    fk= eval(subs(y,x,k));
    fk1= eval(subs(y,x,k+1));
    if fk*fk1<0
        fprintf('\nThe values are %f and %f', k, k+1);
        a=k;
        b=k+1;
        break
    end
end

% Finding Functional Value
fa = eval(subs(y,x,a));
fb = eval(subs(y,x,b));
 
a1=a;
b1=b;

a2=a;
b2=b

n=0;
n1=0;

fprintf('\nImplementing Bisection Method');
if fa*fb > 0 
    disp('Given initial values do not bracket the root.');
else
    c = (a+b)/2;
    fc = eval(subs(y,x,c));
    fprintf('\n\n\tn\t\ta\t\t\tb\t\t\tc\t\t\tf(c)\n');
    while abs(a-b)>e
        fprintf('%f\t%f\t%f\t%f\t%f\n',n,a,b,c,fc);
        if fa*fc< 0
            b =c;
            n=n+1;
        else
            a =c;
            n=n+1;
        end
        c = (a+b)/2;
        fc = eval(subs(y,x,c));
    end
    fprintf('\nRoot(bisection) is: %f\n', c);
end



% Finding Functional Value
fa1 = eval(subs(y,x,a1));
fb1 = eval(subs(y,x,b1));

fprintf('\nImplementing Regula Falsi Method');
if fa1*fb1 > 0 
    disp('Given initial values do not bracket the root.');
else
    c1 = a1 - (a1-b1) * fa1/(fa1-fb1);
    fc1 = eval(subs(y,x,c1));
    fprintf('\n\n\tn1\t\ta1\t\t\tb1\t\t\tc1\t\t\tf(c1)\n');
    while abs(fc1)>e
        fprintf('%f\t%f\t%f\t%f\t%f\n',n1,a1,b1,c1,fc1);
        if fa1*fc1< 0
            b1 =c1;
            fb1 = eval(subs(y,x,b1));
        else
            a1=c1;
            fa1 = eval(subs(y,x,a1));
        end
        n1=n1+1;
        c1 = a1 - (a1-b1) * fa1/(fa1-fb1);
        fc1 = eval(subs(y,x,c1));
    end
    fprintf('\nRoot (Regula Falsi) is: %f\n', c1);
end

fprintf('\nThe root of the equation calculated using inbuilt func: ');
p=vpasolve(y,x,[a2 b2]);
p1=double(p)

err= abs(((c-p1)/p1) *100);
fprintf('\nRelative error using Bisection method: %f ', err);

err1= abs(((c1- p1)/p1) *100);
fprintf('\nRelative error using Regula Falsi method: %f ', err1);