% internal energy
N = 7;
derivativeF = zeros(N);
derivativeF(1:2, 1:2) = [1 -1; -1 1]; % how?

% alpha matrix
alpha = zeros(N);
for i=0:N-1
    alpha = alpha + circshift(derivativeF, [i i]);
end

% creating tri-diagonal branded matrix:
r = [2 -1 zeros(1,N-2)];
alpha2 = toeplitz(r);
alpha2(1, 1) = 2;
alpha2(1, N) = -1;
alpha2(N, 1) = -1;
alpha2(N, N) = 2;
disp(alpha2);

r2 = [6 -4 1 zeros(1,N-3)];
beta = toeplitz(r2); % <-- creates a diagonal branded matrix
% update the corner values
beta(1, 1) =  6;
beta(1, N) = -4;
beta(N, 1) = -4;
beta(N, N) =  6;
beta(1, N-1) = 1;
beta(2, N) = 1;
beta(N-1, 1) = 1;
beta(N, 2) = 1;
disp(beta);
% disp(alpha2+eye(N));
% disp(inv(alpha2+eye(N)));

% X = [1; 2; 3; 4; 5; 6; 7]
% Y = [1; 3; 7; 11; 7; 3; 1]
Y = [2.9; 2.9; 2.9; 2.9; 2.9; 2.9; 2.9];
% Y = [3; 3; 3; 3; 3; 3; 3]
Y_ext = [3; 3; 3; 7; 3; 3; 3];
% X_new = inv(alpha2+eye(N))*X;
% Y_int = inv(alpha2+eye(N))*Y_ext;
% plot(Y, 'b-','DisplayName','Y');
% hold all
% plot(Y_ext, 'r-','DisplayName','Y\_ext');
% plot(Y_int, 'g-','DisplayName','Y\_int');
% xlabel('X');
% ylabel('Y');
% legend
% grid on

% disp(X_new);
% disp(isequal(alpha, alpha2));


