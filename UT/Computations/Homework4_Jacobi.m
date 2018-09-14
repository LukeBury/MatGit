clc
clear

a = [10 2 -1; -3 6 2; 1 1 5]; b = [27;-61.5;-21.5];
n = size(a,1);

% The Jacobi method to solve Ax = b
x_old = zeros(n,1);
x_new = zeros(n,1);
eps = 100;
tol = .05;

while eps > tol
    
    for k=1:n
       
       sum = 0;
       for j = 1:n
            
           if( j ~= k)
                sum = sum + a(k,j)*x_old(j);
           end
           
           x_new(k) = 1/a(k,k) * ( b(k) - sum);
           
       end
    end
    
    eps = norm(x_new-x_old,2) / norm(x_new,2)
    x_old
    x_new
    x_old = x_new;
    
end

x_new