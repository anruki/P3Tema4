function [Q, R] = QRFact(A)
    [f, ~] = size(A);
    I = eye(f);
    Q = I;
    for i=1: f-1
        a = A(i:end, i);
        u = a - norm(a) * I(i:end, i);
        u = (-1)^(a(1)<0) * u / norm(u);
        Hi = I(i:end, i:end) - 2 * (u * u');
        H = I;
        H(i:end, i:end) = Hi;
        Q = Q * H;
        A = H * A;
    end
    R = A;
end