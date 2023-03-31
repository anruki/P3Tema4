function [lambda, x] = Potencia_inv(A, tol)
% Función que permite calcular el autovalor de menor valor absoluto y su 
% autovector asociado.
% El algoritmo utilizado es el Método de la potencia inversa.
% REQUERIDO: Se usan la funciones LUCrout, SubsAdel, SubsAtras.
% Inputs: 
%   A = matriz cuadrada con autovalores distintos reales y distintos de cero
%   tol = tolerancia para detener el algoritmo
% Outputs:
%   lambda = autovalor de menor magnitud (solución númerica)
%   x = autovector asociado en forma de vector columna (solución númerica)
% El error se calculará como la norma infinito del vector diferencia entre
% iteraciones.
    % El procedimiento consiste en calcular el autovalor de mayor magnitud 
    % de la matriz inversa, que por propiedades de los autovalores sabemos 
    % que será el inverso del de menor magnitud de la matriz inicial. 
    [f, ~] = size(A);
    % Inicialización del vector estimación
    xest = ones(f, 1);
    % El vector solución "xs" sin normalizar se obtiene de hacer
    % xs = inv(A)*xest, no obstante podemos expresar xs como solución del
    % sistema A*x = xest y resolver el sistema usando el Método LUCrout
    % (Es ineficiente el cálculo de inversa)
    % Nótese que no podemos trabajar con una matriz que tenga autovalor 0
    [L, U] = LUCrout(A);
    y = SubsAdel(L, xest);
    x = SubsAtras(U, y);
    % Normalizamos el vector resultante
    % de forma que su componente (de mayor magnitud) sea 1
    [~, index] = max(abs(x));
    lambda = x(index);
    x = x / lambda;
    i = 1;
    % El proceso se repite hasta que la norma infinito del vector 
    % diferencia entre iteraciones sea menor que la tolerancia exigida
    while norm(x - xest, 'Inf') > tol
        xest = x;
        y = SubsAdel(L, xest);
        x = SubsAtras(U, y);
        [~, index] = max(abs(x));
        lambda = x(index);
        x = x / lambda;
        i = i+1;
    end
    % El de mayor magnitud de la inversa será el inverso del de menor 
    % magnitud de la matriz inicial
    lambda = 1 / lambda;
end