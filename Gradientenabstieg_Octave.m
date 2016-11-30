
function [xopt,fopt,niter,gnorm,dx] = grad_descent(varargin)
% Dieses m-file zeigt, wie man die Gradientenabstiegsmethode verwendet werden kann, um ein simples Optimierungsproblem zu lösen. 
% Große Schrittweiten (alpha > 0.32) können zur Instabilität des Algorithmus führen.


if nargin==0
    % Definition des Startpunktes.
    x0 = [1.6007 1.2504]';
elseif nargin==1
    % Falls nur ein Eingabe-Argument eingegeben wird, handelt es sich um eine Benutzerdefinierte Eingabe.
    x0 = varargin{1};
else
    error('Incorrect number of input arguments.')
end

% Abbruchtoleranz: Unterschied zum vorhergehenden Wert.
tol = 1e-6;

% Maximale Anzahl der Interationen.
maxiter = 10000;

% Minimale Veränderung.
dxmin = 1e-6;

% Schrittgröße ( 0.33 verursacht Instabilität, Werte zwischen 0 und 0.2 funktionieren gut).
alpha = 0.1;

% Vorbereitung der Norm des Gradienten ||g||, Optimierungsvektor, Iterationszähler und Änderungen.
gnorm = inf; x = x0; niter = 0; dx = inf;

% Definition der Funktion.
f = @(x1,x2) (4.45*exp(-((-2.24-x1)/x2).^2)-0.26).^2+(4.45.*exp(-((0.4-x1)/x2).^2)-1.74).^2+(4.45.*exp(-((1.92-x1)/x2).^2)-4.5).^2+(4.45.*exp(-((2.26-x1)/x2).^2)-3.3).^2+(4.45.*exp(-((4.1-x1)/x2).^2)-0.16).^2;

% Plot der Funktion-Kontour zur Visualisierung.
figure(1); clf; ezcontour(f,[-5 5 -5 5]); axis equal; hold on

% Umdefinierung des Funktionssyntax zur Verwendung in der Optimierung.
f2 = @(x) f(x(1),x(2));

% Gradientenabstiegs-Algorithmus.
while and(gnorm>=tol, and(niter <= maxiter, dx >= dxmin))
    % Berechnung des Gradienten.
    g = grad(x);
    gnorm = norm(g);
    % Schritt.
    xnew = x - alpha*g;
    % Schritt überprüfen.
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
    % Plotte den aktuellen Punkt.
    plot([x(1) xnew(1)],[x(2) xnew(2)],'ko-')
    refresh
    % Aktualisiere Abbruchdaten.
    niter = niter + 1;
    dx = norm(xnew-x);
    x = xnew;
    
end
xopt = x;
fopt = f2(xopt);
niter = niter - 1;
end

% Definition des Gradienten der Funktion.
function g = grad(x)
g = [(17.8.*(4.1-x(1)).*exp(-(4.1-x(1)).^2/x(2).^2).*(4.45.*exp(-(4.1-x(1)).^2/x(2).^2)-0.16))/x(2).^2+(17.8.*(2.26-x(1)).*exp(-(2.26-x(1)).^2/x(2).^2).*(4.45.*exp(-(2.26-x(1)).^2/x(2).^2)-3.3))/x(2).^2+(17.8.*(1.92-x(1)).*exp(-(1.92-x(1)).^2/x(2).^2).*(4.45.*exp(-(1.92-x(1)).^2/x(2).^2)-4.5))/x(2).^2+(17.8.*(0.4-x(1)).*exp(-(0.4-x(1)).^2/x(2).^2).*(4.45.*exp(-(0.4-x(1)).^2/x(2).^2)-1.74))/x(2).^2+(17.8.*(-x(1)-2.24).*exp(-(-x(1)-2.24).^2/x(2).^2).*(4.45.*exp(-(-x(1)-2.24).^2/x(2).^2)-0.26))/x(2).^2
    (17.8.*(4.1-x(1)).^2.*exp(-(4.1-x(1)).^2/x(2).^2).*(4.45.*exp(-(4.1-x(1)).^2/x(2).^2)-0.16))/x(2).^3+(17.8.*(2.26-x(1)).^2.*exp(-(2.26-x(1)).^2/x(2).^2).*(4.45.*exp(-(2.26-x(1)).^2/x(2).^2)-3.3))/x(2).^3+(17.8.*(1.92-x(1)).^2.*exp(-(1.92-x(1)).^2/x(2).^2).*(4.45.*exp(-(1.92-x(1)).^2/x(2).^2)-4.5))/x(2).^3+(17.8.*(0.4-x(1)).^2.*exp(-(0.4-x(1)).^2/x(2).^2).*(4.45.*exp(-(0.4-x(1)).^2/x(2).^2)-1.74))/x(2).^3+(17.8.*(-x(1)-2.24).^2.*exp(-(-x(1)-2.24).^2/x(2).^2).*(4.45.*exp(-(-x(1)-2.24).^2/x(2).^2)-0.26))/x(2).^3];
end

% Ausführen der Funktion.
[xopt,fopt,niter,gnorm,dx] = grad_descent