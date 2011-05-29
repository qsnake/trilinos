% very simple MALTAB script that creates a sparse matrix from the gallery,
% then computes 100 matrix-vector products, and reports the CPU time.
clear all
format long
n = 1000;
A = gallery('poisson', n);
size(A)
x = rand(n * n,1);
t1 = cputime;
for i=1:100
  b = A * x;
end
cputime - t1

