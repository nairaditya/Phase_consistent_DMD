function overlapM = FOVoverlap(strin)

% strin is an n x 1 structure of strin(i).x and strin(i).y of i = 1:n FOVs.
% Fields x and y must be meshgrid format
% overlapM is a n x n structure. overlapM(i,j).m is a matrix the size of
% FOV i with ones where it overlaps FOV j, zeros elsewhere, Empty where 
% i==j

for i = 1:length(strin) % i is the FOV in question
    for j = 1:length(strin) % j determines the FOV we look to overlap
        
        if i == j
            overlapM(i,j).m = [];
        else
            overlapM(i,j).m = zeros(size(strin(i).x));
            rx = [min(strin(j).x(:)),max(strin(j).x(:))];
            ry = [min(strin(j).y(:)),max(strin(j).y(:))];
            overlapM(i,j).m(strin(i).x >= rx(1) & strin(i).x <= rx(2) & strin(i).y >= ry(1) & strin(i).y <= ry(2)) = 1;
        end
    end
end