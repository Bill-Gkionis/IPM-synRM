function dq = alphabeta_to_dq(alpha,beta,angle)

theta = angle;
alphabeta = [alpha, beta]';
Park_matrix = [cosd(theta) sind(theta); -sind(theta) cosd(theta)];

dq = Park_matrix*alphabeta;
end