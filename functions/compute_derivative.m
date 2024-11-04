function [y_prime]=compute_derivative(x,y,Model,disable_edge_points)
%%% computes the derivative: dy/dx
%%% delta x has to be constant !!!

%%%% for debugging use: x=[1:10]; y=x.^3 + x.^2 + 3; y1=3*x.^2 +2*x; y_prime=Vlad_compute_derivative(x,y,'5 points stencil',1);

if size(x,1)==1; dim=size(x,2); else; dim=size(x,1);end
    dx=diff(x);
    dx=dx(1);
    y_prime=zeros(size(x));
    switch Model
        case '2 points stencil'
            for n=1:dim
                if n==1                 % first point: forward differences
                    if disable_edge_points ;else; y_prime(n)=(y(n+1)-y(n))/dx; end
                elseif n==dim       % last point: backward differences
                    if disable_edge_points ;else; y_prime(n)=(y(n)-y(n-1))/dx; end
                else                    % middle points: central scheme
                    y_prime(n)=(y(n+1)-y(n-1))/(2*dx);
                end
            end
        case '5 points stencil'
            for n=1:dim
                % first two and last two points are disabled
                if n==1
                    if disable_edge_points ;else; y_prime(n)=(y(n+1)-y(n))/dx; end
                elseif n==2
                    if disable_edge_points ;else; y_prime(n)=(y(n+1)-y(n-1))/(2*dx); end
                elseif n==dim
                    if disable_edge_points ;else; y_prime(n)=(y(n)-y(n-1))/dx; end
                elseif n==dim-1
                    if disable_edge_points ;else; y_prime(n)=(y(n)-y(n-1))/dx; end
                else                     % middle points: central scheme
                    y_prime(n)=(-y(n+2) + 8*y(n+1) - 8*y(n-1) + y(n-2))/(12*dx);
                end
            end
    end
end