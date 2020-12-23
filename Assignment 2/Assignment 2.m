% Q1
%%
clc
d = [1 2 3 4];
D = diag(d);
disp (D);

%%
clc
A = [1:4;5:8;9:12];
x = diag(A);
disp (A);
disp (x);




%%
%Q2
clc;
L1 = linspace(1, 20);
disp(L1);

%%
clc;
L1 = linspace(1, 20, 5);
disp(L1);

%%
clc;
L1 = logspace(2, 4);
disp(L1);

%%
clc;
L1 = logspace(1, 4, 8);
disp(L1);

%%
clc;
L1 = 2:10;
disp(L1);

%%
clc;
L1 = 2:1.5:20;
disp(L1);





%%
%Q3
clc
A1 = ones(6,4); %populate with ones A1 = 6x4
disp (A1);
A2 = ones(5); %populate with ones A2 = 5x5
disp (A2);

%%
clc
B = zeros(4,4); %populate with zeros B = 4x4
disp (B);

%%
clc
C = eye(5); %eye-dentity matrix C=5x5
disp (C);

%%
clc
D = rand(8); %random 2x2 Matrix
disp (D);
disp(D(:,3));
disp(D(2,:));
disp(D(2,3));
D(3:6,3:6)=zeros(4,4);
disp(D);




%%
%Q4
format 
clc;
clear all;
A = rand(8);
disp(A);
rm = [];
cm = [];
for i=1:8
    rm(i) = max(A(i,:));
end
disp (rm(:)); %row max
for i=1:8
    cm(i) = max(A(:,i));
end
disp(cm); %column max

if max(cm)>max(rm)
    M = max(cm);
else
    M = max(rm);
end
disp (M);





%%
%Q5
n = input('Enter value of n:');
M = magic(n);
rsum = M(1,:);
csum = M(:,1);
diagsum = 0;
adiagsum = 0;
flag = 0;
for i=1:n-1
    rsum = sum(M(i,:));
    rsum1 = sum(M(i+1,:));
    if rsum == rsum1
       flag = 0;
    else
       flag = flag + 1;
    end
end
for i=1:n-1
    csum = sum(M(:,i));
    csum1 = sum(M(:,i+1));
    if csum == csum1
       flag = 0;
    else
       flag = flag + 1;
    end
end
for i=1:n
    for j=1:n
        if i == j
            diagsum = diagsum + M(i,j);
        end
    end
end
for i=1:n
    for j=1:n
        if i+j == n+1
            adiagsum = adiagsum + M(i,j);
        end
    end
end
if diagsum == adiagsum
   flag = 0;
else
   flag = flag + 1;
end

disp(M);
disp(rsum);
disp(csum);
disp(diagsum);
disp(adiagsum);
if flag ~= 0
    disp('It is not a Magic Matrix');
else
    disp('Magic Matrix Verified!');
end



%%
%Q6
A=rand(3);
I=eye(3);
B=A^(-1);
C=I/A;

if B==C
    disp('true');
else
    disp('false');
end

B=A.^(-1);
C=I./A;

if B==C
    disp('True');
else
    disp('False');
end


%%
%Q7
p = [1 2 3 4 5];
disp((length(p):-1:1).*p);



%%
%Q8
n = input('Enter value of n: ');
A = eye(n);
for i=1:n
    for j=1:n
        if i==j || j==n
           A(i,j) = 1;
        elseif i>j
            A(i,j) = -1;
        else
            A(i,j) = 0;
        end
    end
end
disp(A)




%%
%Q9
clc
%%
n = input('Enter n: ');
v1 = zeros(n,1)-2;
v2 = zeros(n-1,1)+1;
v3 = zeros(1,1)+1;
D1 = diag(v1);
D2 = diag(v2,1);
D3 = diag(v2,-1);
D4 = diag(v3, n-1);
D5 = diag(v3, 1-n);
D = D1+D2+D3+D4+D5;
disp(D);
%%
a = [-2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
D = toeplitz(a);
disp(D);
%%
format rat
a = [1 2 3 4 5 6 7 8];
X = toeplitz(a);
for i=1:8
    for j=1:8
        if i>j
            X(i,j) = 0;
        end
    end
end
disp(X);
%%
format rat
a = [1 1/2 1/3 1/4 1/5 1/6 1/7 1/8];
X = toeplitz(a);
disp(X);

