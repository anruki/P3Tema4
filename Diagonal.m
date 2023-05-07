function [lambda,i] = Diagonal(A, tol)
% Función que permite calcular autovalores de una matriz mediante
% iteraciones con factorización QR.
% REQUERIDO: Se usa la función QRFact.
% Input:
%   A = matriz cuadrada que queremos factorizar
%   tol = error máximo calculado
% Outputs:
%   lambda = vector columna que contiene los autovalores
%   i = número de iteraciones
% El error se calcula como la norma infinito de la diferencia entre los
% vectores columna de autovalores entre iteraciones.
    lambda_est = diag(A);
    [Q,R] = QRFact(A);
    % La nueva A será semejante a la A de la anterior iteración por ser Q
    % (que es ortogonal) la matriz de paso.
    A = R*Q;
    lambda = diag(A);
    i = 1;
    % Se repite hasta que la norma infinito de la diferencia entre 
    % iteraciones sea menor que la tolerancia exigida
    while norm(lambda - lambda_est, 'Inf') >= tol
        lambda_est = lambda;
        [Q,R] = QRFact(A);
        A = R*Q;
        lambda = diag(A);
        i = i + 1;
    end
end
