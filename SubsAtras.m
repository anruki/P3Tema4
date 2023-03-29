function y = SubsAtras(A,b)
    % Función que utiliza el método de sustitución hacia atrás
    % para obtener la solución de un sistema cuya matriz de
    % coeficientes es triangular inferior. 
    % INPUTS:
    %   A = matriz triangular inferior
    %   b = vector columna de términos independientes
    % OUTPUTS
    %   y = vector columna solución del sistema
    Ub=[A,b];
    [f, c] = size(Ub);
    y = zeros(f, 1);
    % Despejamos la última componente de la solución
    y(f) = Ub(f, c) / Ub(f, f);
    for i = f-1:-1:1
        % Se realiza sustitución hacia atrás
        y(i)= (Ub(i, c) - Ub(i, i+1:f)*y(i+1:f)) / Ub(i, i);
    end
end
