function ed = EuclideanDistance(x1, x2)
% x1, x2 are vectors
% e.g. x1 = [1,1], x2 = [2,2]
% ed is the distance between the two points

ed = sqrt((x1(1) - x2(1))^2 + (x1(2) - x2(2))^2);

end