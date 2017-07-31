function s = bisection(f,a,b)

x_left = a;
x_right = b;
j_count = 0;
while (1)
    j_count = j_count + 1;
    s = (x_left+x_right)/2;
    if f(s)<0
        if f(x_left)<0
            x_left = s;
            x_right =b;
        else
            x_left = a;
            x_right = s;
        end
    else
        if f(x_left)>0
            x_left = s;
            x_right =b;
        else
            x_left = a;
            x_right = s;
        end
    end
    if ((f(s)<0.000001) || (j_count>100))
        break
    end
end