function T = TDMA(A, B, C, D, N)

T = zeros(1,N);
T(1:N) = 300;
% T = T_Old;
P = zeros(1,N);
Q = zeros(1,N);

P(1) = B(1)/A(1);
Q(1) = D(1)/A(1);

for jj = 2:N
    
    P(jj) = B(jj)/(A(jj) - C(jj)*P(jj-1));
    Q(jj) = (D(jj) + C(jj)*Q(jj-1))/(A(jj) - C(jj)*P(jj-1));
    
end

T(N) = Q(N);
for ii = N-1:-1:1
    
    T(ii) = P(ii)*T(ii+1) + Q(ii);
    
end

end