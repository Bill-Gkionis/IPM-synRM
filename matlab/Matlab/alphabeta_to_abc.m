function abc = alphabeta_to_abc(alpha,beta)

Clarke_matrix = [2/3 0; -1/3 sqrt(3)/3; -1/3 -sqrt(3)/3];
alphabeta = [alpha, beta]';

abc = (3/2)*Clarke_matrix*alphabeta;
end