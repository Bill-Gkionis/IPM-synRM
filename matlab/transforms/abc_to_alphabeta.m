function alpha_beta = abc_to_alphabeta(a,b,c)

clarke_matrix = [1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2];
abc_vector = [a b c]';

alpha_beta = (2/3)*clarke_matrix*abc_vector;
end