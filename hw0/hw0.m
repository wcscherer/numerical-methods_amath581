%Answers to 1 Newton-Raphson Method for zeros of f(x) = xsin(3x)-exp(x)

% x(n+1) = x(n) + f(x(n))/f'(x(n)) with x1 = -1.6, find root at x ~ -0.5

ERR = 10^(-6);  %  establishing the error thresh-hold for convergence

f=@(x) x.*sin(3*x)-exp(x);
df=@(x) 3*x*cos(3*x)+sin(3*x)-exp(x);

x1 = [-1.6]; % vector of x values
errorn = [1]; % vector of calculation error values
nitr = 1; % number of iterations

while errorn > ERR
    
    x1(end + 1) = x1(end) - f(x1(end))/df(x1(end));
    
    nitr = nitr + 1;
    
    errorn(end+1) = abs(f(x1(end)));
       
end

A1 = x1;
save A1.dat  A1  -ASCII

% checking the zeros using the bisection method

% use the same f definition as above

a = -0.7; b = -0.4;

p = [(a+b)/2];
errb = abs([f(p)]);
nitrb = 1;

   while errb > ERR
    
   if f(a)*f(p(end))<0 
       b = p(end);
   else
       a = p(end);          
   end
   p(end+1) = (a + b)/2; 
   errb = abs(f(p(end)));
   
   nitrb = nitrb + 1;
   end

A2 = p;
save A2.dat  A2  -ASCII

A3 = [nitr, nitrb];
save A3.dat  A3  -ASCII


%Answers to Questions 2 a-i
%Defining initial constants
A = [1, 2; -1, 1]; B = [2, 0; 0, 2]; C = [2, 0, -3; 0, 0, -1];

D = [1, 2; 2, 3; -1, 0]; x = [1; 0]; y = [0;1]; z = [1; 2; -1];


% 2a: calculate A+B
A4 = A + B;
save A4.dat  A4  -ASCII

% 2b: calculate 3x - 4y
A5 = 3.*x - 4.*y;
save A5.dat  A5  -ASCII

% 2c: calculate Ax
A6 = A*x;
save A6.dat  A6  -ASCII

% 2d: calculate B(x - y)
A7 = B*(x-y);
save A7.dat  A7  -ASCII

% 2e: calculate Dx
A8 = D*x;
save A8.dat  A8  -ASCII

% 2f: calculate Dy + z
A9 = D*y + z;
save A9.dat  A9  -ASCII

% 2g: calculate A*B
A10 = A*B;
save A10.dat  A10  -ASCII

% 2h: calculate B*C
A11 = B*C;
save A11.dat  A11  -ASCII

% 2i: calculate C*D
A12 = C*D;
save A12.dat  A12  -ASCII

