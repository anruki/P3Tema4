function [Q, R] = qrmodgrsch(V)
% Función que factoriza una matriz en el producto de una ortogonal por una 
% triangular superior.
% El algoritmo utilizado es el Método de factorización QR por Gram-Schmidt.
% Inputs
%   V = matriz cuadrada
% Outputs:
%   Q = matriz ortogonal (columnas forman una base ortonormal)
%   R = matriz triangular superior
    % Comprobamos que las columnas de V son linealmente independientes
    [m, n] = size(V);
    if rank(V) > n
        disp('Las columnas de la matriz no son linealmente independientes')
        return
    end
    % Inicializar la matriz Q y R
    Q = zeros(m,n);
    R = zeros(n,n);
    % Ortonormalización de Gram-Schmidt modificado
    for j = 1:n
    % Proyectar la j-ésima columna de V sobre el espacio generado por las
    % primeras j-1 columnas de Q
        v = V(:,j);
        for i = 1:j-1
            R(i,j) = Q(:,i)'*V(:,j);
            v = v - R(i,j)*Q(:,i);
        end
        R(j,j) = norm(v);
        Q(:,j) = v/R(j,j);
    end
end
