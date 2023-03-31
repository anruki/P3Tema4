function [lambda, x] = Potencia(A, tol)
    [f, ~] = size(A)
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
    % Vectores auxiliares para lista
    iter=zeros(8,1);
    autoval=zeros(8,1);
    autovect=zeros(8,3);
    % Se guardan en otro vector para posterior lista
    iter(1)=i;
    autoval(1)=lambda;
    xx=x';
    for j=1:3
        autovect(1,j)=xx(1,j);
    end
    % El proceso se repite hasta que la norma infinito del vector 
    % diferencia entre iteraciones sea menor que la tolerancia exigida
    while norm(x - xest, 'Inf') > tol
        i= i+1;
        xest = x;
        x = A * xest;
        [~, index] = max(abs(x));
        lambda = x(index);
        x = x / lambda;
        % Se guardan los valores para posterior lista      
        iter(i)=i;
        autoval(i)=lambda;
        xx=x';
        for j=1:3
            autovect(i,j)=xx(1,j);
        end        
    end
  T=table(iter,num2str(autoval,'%.4f'),num2str(autovect,'   %.4f   '),'VariableNames',["NUM. ITERACION","AUTOVALOR","AUTOVECTOR"]);
  disp(T)
end
