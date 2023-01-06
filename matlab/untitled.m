% Initialize the plot
x = -10:10;
y = x.^2;
plot(x,y,'o','MarkerFaceColor','blue')
title('Quadratic function')

while true
   % Select a point on the plot using the mouse
   [x,y] = ginput(1);

   % Find the index of the nearest point to the selected point
   [~,idx] = min((x-x').^2+(y-y').^2);

   % Change the color of the nearest point to red
   hold on
   plot(x(idx),y(idx),'o','MarkerFaceColor','red')
   hold off
end