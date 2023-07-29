function O = soft_orth(x)

O = (sum(sum(abs(x'*x))) - trace(x'*x));
