function alphabeta = dq_to_alphabeta(d,q,angle)

theta = angle;
dq = [d,q]';
Park_matrix = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];

alphabeta = Park_matrix*dq;
end