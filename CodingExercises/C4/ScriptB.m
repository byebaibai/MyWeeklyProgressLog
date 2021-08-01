img = imread('../data/amsterdam.bmp'); 
figure(1);
imshow(img);
hold on;
line([385,348],[340,77],'linewidth', 3,'color','r');
plot(379, 480,'*', 'Color', 'm','MarkerSize', 10, 'linewidth', 2);
rgbMax = [-1 -1 -1];
for ri = 1:size(img, 1)
    for ci = 1:size(img, 2)
        pixel = img(ri, ci, 1:3);
        rgbMax(1) = max(rgbMax(1), pixel(1));
        rgbMax(2) = max(rgbMax(2), pixel(2));
        rgbMax(3) = max(rgbMax(3), pixel(3));
    end
end
disp(rgbMax)
title('RGB Original Image')
circle(379, 480, 100, rgbMax);

function circle(x, y, r, color)
    %x and y are the coordinates of the center of the circle
    %r is the radius of the circle
    %0.01 is the angle step, bigger values will draw the circle faster but
    %you might notice imperfections (not very smooth)
    ang=0:0.01:2*pi; 
    xp=r*cos(ang);
    yp=r*sin(ang);
    plot(x+xp,y+yp, 'Color', color/255, 'linewidth', 2);
end