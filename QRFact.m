function [Q, R] = QRFact(A)
% Función que factoriza una matriz en una ortogonal y una triangular
% superior.
% El algoritmo utilizado es el Método de factorización QR for reflexiones 
% de Housholder.
% Inputs:
% A = matriz cuadrada de la que se quiere obtener su factorización
% Ouputs:
% Q = matriz ortogonal (columnas forman una base ortonormal)
% R = matriz triangular superior
    [f, ~] = size(A);
    I = eye(f);
    Q = I;
    for i=1: f-1
        a = A(i:end, i);
        u = a - norm(a) * I(i:end, i);
        % Se cambia el signo con el fin de evitar el mal condicionamiento
        % del problema
        u = (-1)^(a(1)<0) * u / norm(u);
        % Matriz de housholder que refleja respecto del subespacio 
        % generado por ei (ei es un vector de la base canónica)
        Hi = I(i:end, i:end) - 2 * (u * u');
        H = I;
        H(i:end, i:end) = Hi;
        % La matriz ortogonal Q será el producto de matrices de housholder
        Q = Q * H;
        A = H * A;
    end
    R = A;
end
