function y = SubsAdel(A,b)
    % Función que utiliza el método de sustitución hacia adelante
    % para obtener la solución de un sistema cuya matriz de
    % coeficientes es triangular superior. 
    % INPUTS:
    %   A = matriz triangular superior
    %   b = vector columna de términos independientes
    % OUTPUTS
    %   y = vector columna solución del sistema
    Lb=[A,b];
    [f, c] = size(Lb);
    y = zeros(f, 1);
    % Despejamos la primera componente de la solución
    y(1) = Lb(1, c) / Lb(1, 1);
    for i = 2: f
        % Se realiza sustitución hacia delante
        y(i)= (Lb(i, c) - Lb(i, 1:i-1)*y(1:i-1)) / Lb(i, i);
    end
end
