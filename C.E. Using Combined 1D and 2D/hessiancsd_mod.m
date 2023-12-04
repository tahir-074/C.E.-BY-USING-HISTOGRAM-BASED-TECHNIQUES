% sparse tridiagonal hessian computation with complex step differentiation
% modified version of the following submission
% http://www.mathworks.com/matlabcentral/fileexchange/18177
% By Yi Cao at Cranfield University, 02/01/2008
function [A,g] = hessiancsd_mod(fun,x)
    n=numel(x);                         % size of independent
    A=sparse(n,n);                      % allocate memory for the Hessian matrix
    g=sparse(n,1);                      % allocate memory for the gradient matrix
    h=sqrt(eps);                        % differentiation step size
    h2=h*h;                             % square of step size
    %% determine the hessian and gradient
    for k=1:n-1                         % loop for each independent variable 
        x1=x;                           % reference point
        x1(k)=x1(k)+h*1i;                % increment in kth independent variable
        if nargout>1
            v=fun(x1);                  % function call with a comlex step increment
            g(k)=imag(v)/h;             % the kth gradient
        end
        % only tridiagonal
        for l=k:k+1                     % loop for off diagonal Hessian
            x2=x1;                      % reference with kth increment
            x2(l)=x2(l)+h;              % kth + lth (positive real) increment
            u1=fun(x2);                 % function call with a double increment
            x2(l)=x1(l)-h;              % kth + lth (negative real) increment
            u2=fun(x2);                 % function call with a double increment
            A(k,l)=imag(u1-u2)/h2/2;    % Hessian (central + complex step)
            A(l,k)=A(k,l);              % symmetric
        end  
    end
    %% determine the final point A(n,n)
    k=n;
    x1=x;                           % reference point
    x1(k)=x1(k)+h*1i;               % increment in kth independent variable
    if nargout>1
        v=fun(x1);                  % function call with a comlex step increment
        g(k)=imag(v)/h;             % the kth gradient
    end
    l = k;
    x2=x1;                      % reference with kth increment
    x2(l)=x2(l)+h;              % kth + lth (positive real) increment
    u1=fun(x2);                 % function call with a double increment
    x2(l)=x1(l)-h;              % kth + lth (negative real) increment
    u2=fun(x2);                 % function call with a double increment
    A(k,l)=imag(u1-u2)/h2/2;    % Hessian (central + complex step)
    A(l,k)=A(k,l);              % symmetric
end