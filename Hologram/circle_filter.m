%R is Radius in pixels
%A is the input image
%x0 and y0 are the centre coordinates for the circle
function C = circle_filter(A, R, x0, y0)
S = size(A);
mask = zeros(S);
C = A;
for n = 1:1:S(1)
    for m = 1:1:S(2)
        if((n-y0)^2 + (m-x0)^2 > R^2)
            C(n,m)=0;
        end
    end
end
