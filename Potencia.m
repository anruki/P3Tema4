function [lambda, x] = Potencia(A, tol)
% Función que permite calcular el autovalor de mayor valor absoluto y su 
% autovector asociado.
% El algoritmo utilizado es el Método de la potencia.
% Inputs:
%   A = matriz cuadrada con autovalores distintos reales
%   tol = tolerancia para detener el algoritmo
% Outputs:
%   lambda = autovalor de mayor magnitud (solución númerica)
%   x = autovector asociado en forma de vector columna (solución númerica)
% El error se calculará como la norma infinito del vector diferencia entre
% iteraciones.
    [f, ~] = size(A);
    % Inicialización del vector estimación
    xest = ones(f, 1);
    % El vector solución sin normalizar se obtiene como
    x = A * xest;
    % Normalizamos el vector resultante
    % de forma que su componente (de mayor magnitud) sea 1
    [~, index] = max(abs(x));
    lambda = x(index);
    x = x / lambda;
    i = 1;
    % Impresión por pantalla de la lista con los resultados
    fprintf( ...
        "Iteración %d: x = [%s] lambda = %.4f\n", ...
        i, num2str(x', " %.4f "), lambda);
    % El proceso se repite hasta que la norma infinito del vector 
    % diferencia entre iteraciones sea menor que la tolerancia exigida
    while norm(x - xest, 'Inf') > tol
        xest = x;
        x = A * xest;
        [~, index] = max(abs(x));
        lambda = x(index);
        x = x / lambda;
        i = i+1;
        % Impresión por pantalla de la lista con los resultados
        fprintf( ...
            "Iteración %d: x = [%s] lambda = %.4f\n", ...
            i, num2str(x', " %.4f "), lambda);
    end
end
