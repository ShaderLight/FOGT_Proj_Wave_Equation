%%Przydatne skróty

%komentowanie/odkomentowanie bloku tekstu
%ctrl+R
%ctrl+shift+R

prompt = {'Podaj rozmiar kwadratowej siatki (liczba punktów):'};
dimensionAnswer = inputdlg(prompt, 'Input', [1 50], {'100'});
dim = str2num(dimensionAnswer{1});

prompt = {'Podaj liczbę źródeł:'};
sourceAnswer = inputdlg(prompt, 'Input', [1 50], {'1'});
sourceCount = str2num(sourceAnswer{1});

figure
n = sourceCount;
coordinates = zeros(n,2);
hold on
grid on
axis equal
for i=1:n
axis([0 dim 0 dim])
set(gcf, 'Position',  [700, 200, 1000, 1000])
[x, y] = ginput(1);
coordinates(i,:) = [x, y];
plot(coordinates(:,1), coordinates(:,2), 'o');
axis([0 dim 0 dim])
set(gcf, 'Position',  [700, 200, 1000, 1000])
end
hold off
close
sourceCoords = round(coordinates);


% obstacle input
obstacleMat = ones(dim);
obstacleSize = round(dim/10);
obstacleRadius = round(obstacleSize/2);

prompt = {'Podaj liczbę przeszkód:'};
obstacleAnswer = inputdlg(prompt, 'Input', [1 50], {'0'});
obstacleCount = str2num(obstacleAnswer{1});

figure
n = obstacleCount;

coordinates = ones(n*obstacleSize^2,2);
hold on
grid on
axis equal
for i=1:n
axis([0 dim 0 dim])
set(gcf, 'Position',  [700, 200, 1000, 1000])
[x, y] = ginput(1);

xCandidates = [round(x)-obstacleRadius:round(x)+obstacleRadius];
yCandidates = [round(y)-obstacleRadius:round(y)+obstacleRadius];

coordCandidates = combvec(xCandidates, yCandidates)';

for k=1:length(coordCandidates)
    if coordCandidates(k, 1) < 0
        continue;
    elseif coordCandidates(k, 1) > dim
        continue;
    elseif coordCandidates(k, 2) < 0
        continue;
    elseif coordCandidates(k, 2) > dim
        continue;
    end

    coordinates(i*obstacleSize^2 + k,:) = [coordCandidates(k, 1), coordCandidates(k, 2)];
    plot(coordinates(i*obstacleSize^2 + k,1), coordinates(i*obstacleSize^2 + k,2), 'o');

end


%coordinates(i,:) = [x, y];
%plot(coordinates(:,1), coordinates(:,2), 'o');

axis([0 dim 0 dim])
set(gcf, 'Position',  [700, 200, 1000, 1000])
end
hold off
close

coordinates
obstacleCoords=unique(coordinates,'rows');


for i=1:length(obstacleCoords)
    obstacleMat(obstacleCoords(i, 1), obstacleCoords(i, 2)) = 0;
end



% Zmienne
dx = 0.1;
Lx = dim*dx;
Ly = Lx;
dy = dx;
nx = fix(Lx/dx);
ny = fix(Ly/dy);
x = linspace(0, Lx, nx);
y = linspace(0, Ly, ny);

T=60*8;

wn=zeros(nx, ny);
wnm1=wn; 
wnp1=wn;



CFL = 0.5;
c = 1;
dt = CFL*dx/c;

%Pętla idąca po czasie
t = 0;
set(gcf, 'Position',  [700, 200, 1000, 1000])
while(t < T)
    %odbija się od ścianek 
    wn(:,[1,end])=0;
    wn([1,end],:)=0;

    %nie odbija się od ścianek
     %wnp1(1,:) = wn(2,:) + ((CFL-1)/(CFL+1)) * (wnp1(2,:) - wn(1,:));
     %wnp1(end,:) = wn(end-1,:) + ((CFL-1)/(CFL+1)) * (wnp1(end-1,:) - wn(end,:));
     %wnp1(:,1) = wn(:,2) + ((CFL-1)/(CFL+1)) * (wnp1(:,2) - wn(:,1));
     %wnp1(:,end) = wn(:,end-1) + ((CFL-1)/(CFL+1)) * (wnp1(:,end-1) - wn(:,end));
    

    t=t+dt;
    wnm1=wn;
    wn=wnp1;
    
    for i=1:sourceCount
        wn(sourceCoords(i, 1), sourceCoords(i, 2)) = dt^2*20*sin(30*pi*t/20);
    end    
    
    wn = wn.*obstacleMat;


    for i=2:nx-1
        for j=2:ny-1
            wnp1(i,j) = 2*wn(i,j) - wnm1(i,j) + CFL^2*(wn(i+1,j) + wn(i,j+1) - 4*wn(i,j) + wn(i-1,j) + wn(i,j-1));
        end
    end


    set(gcf, 'name', 'Równanie fali')
    %strona wizualna
    clf;
    subplot(2,1,1);
    hold on;
    axis equal;
    imagesc(x, y, wn');
    colorbar;
    clim([-0.02 0.02])
    title(sprintf('t = %.2f', t));
    hold off;
    subplot(2,1,2);
    mesh(x, y, wn');
    colorbar;
    clim([-0.02 0.02])
    axis([0 Lx 0 Ly -0.05 0.05]);
    shg;
    pause(0.01);
end

