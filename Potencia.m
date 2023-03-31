function [lambda, x] = Potencia(A, tol)
% Función que permite calcular el autovalor de mayor valor absoluto y su 
% autovector asociado dada una matriz de autovalores distintos reales.
% El algoritmo utilizado es el Método de la potencia.
% Inputs: 
% A = matriz cuadrada con autovalores distintos reales
% tol = tolerancia para detener el algoritmo
% Outputs:
% lambda = autovalor de mayor magnitud (solución númerica)
% x = autovector asociado en forma de vector columna (solución númerica)
% El error se calculará como la norma infinito del vector diferencia entre
% iteraciones, deteniéndose el algoritmo cuando este sea menor que la 
% tolerancia exigida.
% Observación: por simplicidad se asume que todos los autovalores son 
% reales y distintos, por tanto no entraremos en el tratamiento de 
% excepciones para estos casos.
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
    % El proceso se repite hasta que la norma infinito del vector 
    % diferencia entre iteraciones sea menor que la tolerancia exigida
    while norm(x - xest, 'Inf') > tol
        xest = x;
        x = A * xest;
        [~, index] = max(abs(x));
        lambda = x(index);
        x = x / lambda;
        % Impresión por pantalla de la lista 
        fprintf('Iteración %d: x = [', i);
        for k = 1:f-1
            fprintf('%.4f, ', x(k));
        end
        fprintf('%.4f], lambda = %.4f\r', x(f), lambda);
        
        i = i+1;
    end
end
