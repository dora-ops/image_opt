function [U,P] = OSP(HIM,p)
% OSP algorithm for endmember extraction.
% ------------------------------------------------------------------------------
% Input:   HIM : hyperspectral image cube [nrows x ncols x nchannels]
%          p   : desired number of endmembers to be extracted
%
% Output:  U   : set of extracted endmembers [nchannels x p] 
%          P   : spatial coordinates of the extracted endmembers (positions
%                rows x cols)
% 
% Copyright (2007) GRNPS group @ University of Extremadura, Spain. 

disp(' === Start OSP run ===')

% Obtener tiempo CPU actual | get current CPU time
t1=cputime;

% Obtener tamaño de la imagen (muestras,lineas,bandas) | get image size
[ns,nl,nb]=size(HIM);

if nb==1
    Uin_row(1,:,:) = HIM';
    HIM = Uin_row;
    [ns,nl,nb]=size(HIM);
end

% Visualizar imagen | visualize image
imagesc(mean(HIM,3)); colormap(gray); 
set(gca,'DefaultTextColor','black','xtick',[],'ytick',[],'dataaspectratio',[1 1 1])
po1 = get(gca,'position');

% Calculo del pixel (vector) con mayor intensidad en la imagen | 
% Calculate the pixel (vector) with major intensity in the image
max = 0;
for i = 1:ns
    for j = 1:nl
        r = squeeze(HIM(i,j,:));
        bright = r'*r;
        if bright > max
            max = bright;
            posx = i;
            posy = j;
        end
    end
end

% El pixel con mas intensidad es el pixel inicial del proceso |
% The pixel with more intensity is the initial pixel of the process
t0 = squeeze(HIM(posx,posy,:));

% Calculo de la matriz identidad | Generate the identity matrix.
I = eye(nb,nb);

% Inicializacion de la matriz de pixels puros |
% Initialization of the pure pixels matrix
U = [];
U = [U t0];

% Inicializacion de la matriz de posiciones |
% Initialization of the positions matrix
P = zeros(p,2);
P(1,1)=posx; P(1,2)=posy;

disp(sprintf('. found pixel @ coordinates x=%5d & y=%5d',posx,posy))

% Visualizacion de la posicion del primer pixel seleccionado |
% Visualization of the position of the first chosen pixel
drawnow;
%text(posy,posx,'o','Margin',1,'HorizontalAlignment','center','FontSize',22,'FontWeight','light','FontName','Garamond','Color','yellow');
text(posy,posx,'o','Color','yellow');

% Algoritmo OSP
for i = 1:p-1
    UC = U(:,1:i);
    % Calculo de la proyeccion ortogonal con respecto a los pixels
    % actualmente seleccionados. Esta parte puede sustituirse por cualquier
    % otra distancia |
    % Calculate the orthogonal projection with respect to the pixels at present chosen. 
    % This part can be replaced with any other distance
    PU = I-UC*pinv(UC'*UC)*UC';
    maximum = 0;
    % Calculo del pixel mas distinto a los ya seleccionados en funcion de la
    % proyeccion ortogonal (o cualquier otra distancia que se seleccione) |
    % Calculate the pixel most different from the already selected ones according to
    % the orthogonal projection (or any other distance selected)
    for n = 1:ns
        for m = 1:nl
            r = squeeze(HIM(n,m,:));
            result = PU*r;
            val = result'*result;
            if (val > maximum) 
                maximum = val;
                posx = n; posy = m;
            end
        end
    end
    % El siguiente pixel seleccionado es el mas diferente a los ya seleccionados |
    % The next chosen pixel is the most different from the already chosen ones
    ti = squeeze(HIM(posx,posy,:));
    % Mostrar posiciones de dicho pixel por pantalla |
    % Show positions of the above mentioned pixel at screen
    disp(sprintf('. found pixel @ coordinates x=%5d & y=%5d',posx,posy))
    % Almacenar posiciones en matriz de posiciones |
    % Store positions in the matrix of positions
    P(i+1,1)=posx; P(i+1,2)=posy;
    % Almacenar pixel en matriz de pixels puros |
    % Store positions in the matrix of pure pixels
    U = [U ti];
    % Visualizar pixel seleccionado en pantalla |
    % Visualize pixel selected on screen
    drawnow;
   % text(posy,posx,'o','Margin',1,'HorizontalAlignment','center','FontSize',22,'FontWeight','light','FontName','Garamond','Color','yellow');
    text(posy,posx,'o','Color','yellow');

end

% Obtener tiempo CPU actual | get current CPU time
t2=cputime;

% Mostrar tiempo total en ejecucion del algoritmo |
% Show total execution time of the algorithm
disp(sprintf('. Total CPU processing time .................... %6.3f [s]  ',(t2-t1)));
disp(' === Eng OSP ===');
