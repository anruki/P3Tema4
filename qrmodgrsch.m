function [Q, R] = qrmodgrsch(V)
% Comprobar que las columnas de V son linealmente independientes
if rank(V) < size(V,2)
    error('Las columnas de la matriz no son linealmente independientes')
end

% Inicializar la matriz Q y R
[m, n] = size(V);
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