function [lambda,i] = DiagonalGS(A, tol)
% Función que permite calcular autovalores de una matriz mediante
% iteraciones con el método de Gram-Smith, evitando el error.
% REQUERIDO: Se usa la función qrmodgrsch.
% Input:
% A = matriz cuadrada que queremos factorizar
% tol = error máximo calculado
% Outputs:
% lambda = vector columna que contiene los autovalores
% i = número de iteraciones
% El error se calcula como la norma infinito de la diferencia entre los
% vectores columna de autovalores entre iteraciones.
    % El algoritmo consiste en hallar matrices semejantes a la anterior 
    % iteración, convergiendo a una matriz triangular superior que por
    % serlo tiene sus autovalores en la diagonal y que por ser semejante
    % son los autovalores de la inicial.
    lambda_est = diag(A);
    [Q,R] = qrmodgrsch(A);
    % La nueva A será semejante a la A de la anterior iteración por ser Q
    % (que es ortogonal) la matriz de paso.
    A = R*Q;
    lambda = diag(A);
    i = 1;
    % Se repite hasta que la norma infinito de la diferencia entre 
    % iteraciones sea menor que la tolerancia exigida
    while norm(lambda - lambda_est, 'Inf') >= tol
        lambda_est = lambda;
        [Q,R] = qrmodgrsch(A);
        A = R*Q;
        lambda = diag(A);
        i = i + 1;
    end
end
