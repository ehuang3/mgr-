%% Generate height data.

f = 1; % frequency
p = zeros(2,1); % pixel location
s = zeros(2,1); % source location
n = 100; % number of spatial points
m = 20; % number of heights per point
x = zeros(m,n); % data
dx = 1; % space between points

n_p = 1000;
P = zeros(2,n_p);

m = 1000;
p = [1000; 0];
s = [500; 500];
c = [0; 0];
dh = f / m;
phase = zeros(1,m + 1);
j = 1;
for i = -m/2:m/2
    c(2) = dh * i;
    dist = norm(c - s) + norm(c - p);
    phase(j) = dist;
    j = j + 1;
end
plot(sin(2*pi*phase))

%% Generate random height field and response.
% Pixel locations.
n_px = 100;
n_py = 100;
px = linspace(-50, 50, n_px);
py = linspace(950, 1050, n_py);
[PX, PY] = meshgrid(px, py);
P = zeros(2,n_px*n_py);
P(1,:) = PX(:);
P(2,:) = PY(:);

% Laser location.
s = [500; 500];

n = 100;
f = 1;
dx = f;
x = linspace(-n*f/2, n*f/2, n);
h_max =  f/2;
h_min = -f/2;

% Random height field.
h = zeros(2, n);
for i = 1:n
    h(1,i) = x(i);
    h(2,i) = h_min + (h_max-h_min)*rand;
end

% Response.
R = zeros(length(P), 1);
for i = 1:length(P)
    for j = 1:length(h)
        dist = norm(h(:,j) - s) + norm(h(:,j) - P(:,i));
        R(i) = R(i) + sin(2*pi*dist);
    end
end

%% Generate A matrix.

n = 100;
f = 1;
dx = f;
x = linspace(-n*f/2, n*f/2, n);
h_max =  f/2;
h_min = -f/2;

m = 20;
A = zeros(length(P), n*m);

%% 

p1 = [-3; 0.9];
p2 = [ 5; 0.2];
c = [0; 1000];
d = [];
j = 1;
for i = -1000:1000
    c(1) = i;
    d(j) = norm(p1-c) - norm(p2-c);
    j = j + 1;
end
hold on
plot(cos(d))

%% 
t = 0:1/100:10-1/100;                     % Time vector
x = cos(2*pi*15*t + pi/3);      % Signal

y = fft(x);                               % Compute DFT of x
m = abs(y);                               % Magnitude
p = unwrap(angle(y));                     % Phase

f = (0:length(y)-1)*100/length(y);        % Frequency vector

subplot(2,1,1)
plot(f,m)
title('Magnitude')
ax = gca;
ax.XTick = [15 40 60 85];

subplot(2,1,2)
plot(f,p*180/pi)
title('Phase')
ax = gca;
ax.XTick = [15 40 60 85];

%% 

c = [0; 1000];
s = [500; 500];
p = [0; 0];
f = 1;
dist = norm(p-s) + norm(c-p);
exp(i*2*pi/f*dist);

%% Calculate electric field for a random height map.

f = 1;

% Camera locations.
n_c = 200;
c = zeros(2,n_c);
c(1,:) = linspace(-n_c/2, n_c/2, n_c);
c(2,:) = 1000;

% Laser location.
s = [500; 500];

% Random height map.
n_x = 100;
x = zeros(2, n_x);
x(1,:) = linspace(-10, 10, n_x);
h_max = f;
for i = 1:n_x
    x(2,i) = h_max * rand;
end

% Electric field.
E = calcE(c, s, x, f);

subplot(2,1,1)
plot(abs(E))

subplot(2,1,2)
plot(unwrap(angle(E)))
hold on
plot(angle(E));

%% Random height offset.
x(2,:) = x(2,:) + 0.1;
E = calcE(c,s,x,f);
subplot(2,1,1)
hold on
plot(abs(E))

subplot(2,1,2)
hold on
plot(unwrap(angle(E)))
plot(50 * angle(E));

%% Fourier transform.
y = fft(abs(E));
m = abs(y);                               % Magnitude
p = unwrap(angle(y));                     % Phase

f = (0:length(y)-1)*100/length(y);        % Frequency vector

subplot(3,1,3)
plot(abs(E))

subplot(3,1,1)
plot(f,m)
title('Magnitude')
ax = gca;
ax.XTick = [15 40 60 85];

subplot(3,1,2)
plot(f,p*180/pi)
title('Phase')
ax = gca;
ax.XTick = [15 40 60 85];