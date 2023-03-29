function [lambda,x]=Potencia(A,tol)
% notas:
% Sólo se busca el mayor valor propio.
% No puede haber otros valores propios con la misma magnitud que el mayor
% valor propio. (multiplicidad debe ser 1)
% El mayor valor propio es un número real.
% La matriz para la que se determina el valor propio no puede modificarse.

    [F, ~] = size(A);
    x = ones(F,1);
    x0 = zeros(F,1);
    i = 0;
    while norm((x - x0), inf) > tol
        x0 = x;
        x = A*x0;
        lambda = abs(max(x));
        x = x / lambda;
        i = i+1;
        fprintf('iteración número %d\n',i)
        fprintf('autovector:');
        x
        fprintf('valor propio:');
        lambda
    end
end
