function [lambda,x] = Potencia_inv(A,tol)
% notas:
% Sólo se busca el menor valor propio.
% No puede haber otros valores propios con la misma magnitud que el mayor
% valor propio. (multiplicidad debe ser 1)
% El mayor valor propio es un número real.
% La matriz para la que se determina el valor propio no puede modificarse.

    [F, ~] = size(A);
    x = ones(F,1);
    x0 = zeros(F,1);
    i = 0;
    [L, U] = LUCrout(A);
    while norm((x - x0), inf) > tol
        x0 = x;
        y = SubsAdel(L,x0);
        x = SubsAtras(U,y);
        lambda = 1 / min(x);
        x = x * lambda;
        i = i+1;
    end
    % Invertimos el autovalor para obtener el menor autovalor de A
    lambda = 1 / lambda;

end
