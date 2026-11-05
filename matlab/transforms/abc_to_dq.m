function dq = abc_to_dq(a,b,c,angle)

theta = angle;
abc = [a,b,c]';
Park_matrix = [cosd(theta) cosd(theta-120) cosd(theta-240); -sind(theta) -sind(theta-120) -sind(theta-240)];

dq = (2/3)*Park_matrix*abc;
end