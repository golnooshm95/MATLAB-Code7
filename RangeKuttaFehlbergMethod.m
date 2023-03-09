function G = RKF(a,b,alfa,Tole,hmax,hmin)
% Author: Golnoosh Morshedi
% Range-Kutta-Fehlberg Method
%   To approximate initial value problem y'=f(t,y), a<t<b and y(a)=alfa
%   You have to pass endpoint, initial value, tolerance, maximum and
%   minimum step
%   Example: y'=2y/t+(exp(t))      1<t<3            y(1)=0
%   RKF(1,3,0,10^-4,0.5,0.02)
clc
format long
t = a;
w = alfa;
h = hmax;
Flag = 1;
f = input('Enter IVP Function: ');
while(Flag==1)
    k1 = h*f(t,w);
    k2 = h*f((t+(1/4)*h),(w+(1/4)*k1));
    k3 = h*f((t+(3/8)*h),(w+((3/32)*k1)+((9/32)*k2)));
    k4 = h*f((t+(12/13)*h),(w+((1932/2197)*k1)-((7200/2197)*k2)+((7296/2197)*k3)));
    k5 = h*f((t+h),(w+((439/216)*k1)-((8)*k2)+((3680/513)*k3)-((845/4104)*k4)));
    k6 = h*f((t+(1/2)*h),(w-((8/27)*k1)+((2)*k2)-((3544/2565)*k3)+((1859/4104)*k4)-((11/40)*k5)));
    R = (1/h)*abs(((1/360)*k1)-((128/4275)*k3)-((2197/75240)*k4)+((1/50)*k5)+((2/55)*k6));
    if R <= Tole
        t = t+h;                         % Approximation accepted.
        w = w+((25/216)*k1)+((1408/2565)*k3)+((2197/4104)*k4)-((1/5)*k5);
        G = [t,w,h]
    end
    q = 0.84*(Tole/R)^(1/4);
    if q <= 0.1
        h = 0.1*h;
    else
        if q >= 4
            h = 4*h;
        else
            h = q*h;
        end
    end
    if h > hmax
        h = hmax;
    end
    if t >= b
        Flag = 0;
    else
        if t+h > b
            h = b-t;
        else
            if h < hmin
                Flag = 0;
                fprintf('Minimum h exceeded.');
                fprintf('Procedure completed unsuccessfully.');
            end
        end
    end
end


