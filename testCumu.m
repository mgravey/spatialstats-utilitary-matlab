maxTemp=[0 0;10 0; 0 20; 20 20]

maxTemp=[0 0; 240 0; 0 240]

im=[
    [0 0 0 0 0 0 0];
    [0 1 0 0 0 1 0];
    [0 0 0 0 0 0 0];
    [0 0 0 1 0 0 0];
    [0 0 0 0 0 0 0];
    [0 1 0 0 0 1 0];
    [0 0 0 0 0 0 0]]

im=imresize(im,40,'nearest');

val=cumulant2(im,maxTemp);

imagesc(flipud(val))
colormap(jet)
axis square