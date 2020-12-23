% Q1
a = 3;
b = 5;
c = -3;
x = b - a/(b+((c+a)/c*a));
disp (x)
fprintf x

%%
%Q2
a = 3;
b = 5;
c = -3;
x = b - a/(b+((c+a)/c*a));
disp (x)
fprintf x

%%
% Q3
n = input('Enter value of n = ');
m = input('Enter value of m = ');
ncm = factorial(n)/(factorial(m)*(factorial(n-m)));
f=1; 
g=1;
k=1;
%Finding n factorial, f = n!
for i=1:n
    f=f*i;
end
%g = m!
for i=1:m
    g=g*i;
end
%k = (n-m)!
for i=1:(n-m)
    k=k*i;
end
ans = f/(g*k);
fprintf('value of nCm (using function): %d\n',ncm);
fprintf('value of nCm: (using loops) %d\n',ans);


%%
% Q4
a = input('Enter value of a = ');
r = input('Enter value of r = ');
n = input('Enter value of n = ');
sum = 0;
if n>100000
    if abs(r)>=1
        fprintf('Series is Divergent.\n');
    else
        sum = a/(1-r);
        fprintf('Series is Convergent.\n');
        fprintf('Sum = %d\n',sum);
    end
else
    for i=1:n
        sum = sum + power(a,i);
    end
    fprintf('Sum of Series = %d\n', sum);
end


%%
% Q5
%clc
%clear all
k = 1000000;
sum = 0;
for n=0:k
    sum=sum+exp(-n);
end

F1 = 1/(1-(1/exp(1)));
fprintf('Calculated Sum = %d\n', sum);
fprintf('Actual Sum = %d\n', F1);



%%
% Q6
n = input ('Enter value of n: ');
p = 1;
for i=1:n
    p = p*(1+(2/i));
end
fprintf('Value of product = %d\n',p);


%%
% Q7
fibonacci = [0 1];
n=10;
sum=0;
for i = 1:n-2
    fibonacci = [fibonacci fibonacci(end)+fibonacci(end-1)];
end
for i = 1:n
    sum = sum + fibonacci(i);
end
disp(fibonacci);
disp (sum);

%%
% Q8
function f = Question8(x)
    if x<0
        f=0;
    elseif x>=0 && x<=1
        f = x;
    elseif x>=1 && x<=2
        f = 2-x;
    elseif x>2
        f=0;
    end
end


%%
% Q9
x=1;
while (x*x*x)<2000
    fprintf('%d, ',(x*x*x));
    x = x+1;
end

%%
% Q10
clc
M = input ('Enter Month: ','s');
Y = input ('Enter Year: ');
%M = {'January' 'February' 'March' 'April' 'May' 'June' 'July' 'August' 'September' 'October' 'November' 'December'};
switch (M)
    case {'January', 'March', 'May', 'July', 'August', 'October', 'December'}
        fprintf('%s - %d\n',M,31);
    case {'April' 'June' 'September' 'November'}
        fprintf('%s - %d\n',M,30);
    case 'February'
        if (mod(Y,400)==0)
            fprintf('%s - 29\n',M);
        elseif (mod(Y,100)==0)
            fprintf('%s - 28\n',M);
        elseif (mod(Y,4)==0)
            fprintf('%s - 29\n',M);
        else
            fprintf('%s - 28\n',M);
        end
end
