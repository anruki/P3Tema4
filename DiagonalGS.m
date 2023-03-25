function [lambda, i] = DiagonalGS(A, tol)
% Función que calcula los autovalores de una matriz utilizando la 
% factorización QR con ortogonalización de Gram-Schmidt 
% (hace uso de la función [Q, R] = qrmodgrsch(V)).
% 
% Input:
%       A: matriz de tamaño n*n
%       tol: error máximo permitido calculado como la norma infinita de la resta
%       de la iteración actual menos la iteración anterior
% 
% Output:
%       lambda: vector de tamaño n que contiene los autovalores de A
%       i: número de iteraciones realizadas

% Inicialización
lambda = diag(A);
V = eye(size(A,1));
error = tol+1;
i = 0;

% Iteraciones
while error > tol
    i = i+1;
    [Q,~] = qrmodgrsch(V'*A*V);
    V = V*Q;
    lambda_ant = lambda;
    lambda = diag(V'*A*V);
    error = norm(lambda-lambda_ant, inf);
end
end