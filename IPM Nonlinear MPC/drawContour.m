function [] = drawContour(bp,x,dx)
%%%% This function is to draw the contour
%%%% Input: bp, x, dx
%%%%    bp - barrier problem
%%%%    x - primal variable
%%%%    dx - deviation of primal variable
%%%% Output: Na

    % generate contour plot
    % design variables at mesh points
    [x1,x2] = meshgrid(1:0.05:5,1:0.05:5);

    % equation 1
    eq1 = x1 .* x2 - 5;
    % equation 2
    eq2 = -x1.^2 - x2.^2 + 20;
    % objective
    obj = x2.*(5+x1);
    for i = 1:size(eq1,1)
        for j = 1:size(eq1,2)
            if (eq1(i,j)<=0.01||eq2(i,j)<=0.01||x1(i,j)>=4.95||x1(i,j)<=1.05||x2(i,j)>=4.95||x2(i,j)<=1.05)
                bobj(i,j) = 0;
            else
                bobj(i,j) = obj(i,j) - bp.mu * log(eq1(i,j)) ...
                    - bp.mu * log(eq2(i,j)) ...
                    - bp.mu * log(x1(i,j)-1) - bp.mu * log(x2(i,j)-1) ...
                    - bp.mu * log(5-x1(i,j)) - bp.mu * log(5-x2(i,j));
            end
        end
    end

    figure(1)
    hold off;
    [C,h] = contour(x1,x2,obj);
    clabel(C,h,'Labelspacing',250);
    hold on;
    lines = linspace(min(min(bobj))+1,max(max(bobj))-1,20);
    [C,h] = contour(x1,x2,bobj,lines);
    clabel(C,h,'Labelspacing',250);
    title('Problem Contour Plot');
    xlabel('x_1');
    ylabel('x_2');
    hold on;
    % solid lines to show constraint boundaries
    [C,h] = contour(x1,x2,eq1,[0,0.1],'r-','LineWidth',1);
    clabel(C,h,'Labelspacing',250);
    [C,h] = contour(x1,x2,eq2,[0,0.1],'b-','LineWidth',1);
    clabel(C,h,'Labelspacing',250);
    % plot step
    plot([x(1),x(1)+dx(1)],[x(2),x(2)+dx(2)],'k-','LineWidth',2);
    plot(x(1),x(2),'ro');
    plot(x(1)+dx(1),x(2)+dx(2),'bx');
    % show a legend
    legend('Objective','Barrier Problem','Constraint','New Step','Start','End');
 
end

