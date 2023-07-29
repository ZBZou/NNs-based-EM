function O = soft_orth(x,N)

O = (sum(sum(abs(x'*x))) - trace(x'*x));
O/(N*(N-1));
