function [N x_tilda pNk] = eval_interpolator_c(tip, eps)
  a = 3;
  f = @(x) exp(a * cos(x))/(2 * pi * besseli(0, a));
  N = 1000;
  x_tilda = linspace(-pi, pi, N + 1);
  y_tilda = f(x_tilda);
  
  % calculul E(p(x), Nk)
  k = 2;
  while (1)
    h = 2 * pi / (N+1);
    Nk = 2^k;
    xk = linspace(-pi, pi, Nk);
    yk = f(xk);
    switch(tip)
      case 1 % interpolare Lagrange
        Nretur = 2^(k - 1);
        pNk = PolinomLagrange(x_tilda, xk, yk);
      case 2 % interpolare Newton
        Nretur = 2^(k - 1);
        pNk = PolinomNewton(x_tilda, xk, yk);
      case 3 % linear spline
        Nretur = 2^(k - 1);
        for i = 1 : N + 1
          pNk(i) = InterpolareLiniara(x_tilda(i), xk, yk);
        endfor
      case 4 % spline cubic natural
        Nretur = 2^(k - 1);
        pNk = C2Natural(x_tilda, xk, yk);
      case 5 % spline cubic tensionat
        Nretur = 2^(k - 1);
        pNk = C2Tensionat(x_tilda, xk, yk);
      case 6 % Fourier
        Nretur = Nk / 2 - 1;
        pNk = TFourier(Nk / 2, x_tilda, f);  
    endswitch  
    E(1) = sqrt(sum((y_tilda - pNk).^2) * h);
    if length(E) == 2
      if E(1) - E(2) > 0 % creste
        N = inf;
        return;
      else   % descreste
        if norm(E(1) - E(2)) < eps % se adinge eroare dorita
          N = Nretur;
          return;       
        endif
      endif
    endif
    E(2) = E(1);
    k++;
  endwhile
  
endfunction

function val = PolinomLagrange(b, x, y)
  n = length(x);
  nb = length(b);
  val = zeros(1, nb);
  
  % calcul valoare polinom Lagrange
  for k = 1 : n
    prod(1:nb) = y(k);   
    for i = 1 : n
      if i ~= k
        prod .*= (b - x(i)) / (x(k) - x(i));
      endif
    endfor
    val += prod;
  endfor
endfunction

function val = PolinomNewton(b, x, y)
  n = length(x);
  
  % calcul diferente divizate
  F(1:n) = y(1:n);
  for i = 1 : n - 1
    for j = n : -1 : i+1
      F(j) = (F(j-1) - F(j)) / (x(j - i) - x(j));    
    endfor
  endfor
  
  % calcul valoare polinom Newton
  nb = length(b);
  prod = ones(1, nb);
  val(1 : nb) = F(1);
  for i = 1 : n - 1
    prod .*= b - x(i);  
    val += prod * F(i + 1);
  endfor
endfunction

function val = InterpolareLiniara(b, x, y)
  n = length(x); 
  for i = 1 : n - 1
    % cauta intervalul in care se afla b
    if b >= x(i) && b <= x(i+1) 
      m = (y(i + 1) - y(i)) / (x(i+1) - x(i));
      n = (x(i+1) * y(i) - x(i) * y(i+1)) / (x(i+1) - x(i));
      val = m * b + n;
      return;
    endif
  endfor
  val = 0;
endfunction

function val = C2Natural(x0, x, y)
  % initializari
  n = length(x);
  a(1:n) = y(1:n);  
  h(1:n-1) = x(2:n) - x(1:n-1);

  % creare matrice tridiagonala
  ma(1:n-2) = h(1:n-2);
  ma(n-1) = 0;

  mb(1) = 1;
  mb(2:n-1) = 2 * (h(1:n-2) + h(2:n-1));
  mb(n) = 1; 

  mc(1) = 0;
  mc(2:n-1) = h(2:n-1);

  g(1) = 0; 
  g(2:n-1) = 3*(a(3:n) - a(2:n-1))/h(2:n-1) - 3*(a(2:n-1) - a(1:n-2))/h(1:n-2);
  g(n) = 0; 

  % rezolvare sistem
  c = Thomas(ma, mb, mc, g);

  % calculare coeficienti ramasi
  for i=1:n-1
    d(i) = (c(i+1) - c(i))/(3 * h(i));
    b(i) = (a(i+1) - a(i))/h(i) - h(i)/3 * (2*c(i) + c(i+1));
  endfor
  
  % calculare valori polinom
  nx = length(x0);
  for k = 1 : nx
    for i=1:n-1
      if x0(k) >= x(i) && x0(k) <= x(i+1)
        val(k) = a(i) + b(i) * (x0(k) - x(i)) + c(i) * (x0(k) - x(i))^2 + d(i)*(x0(k) - x(i))^3;
        break;
      endif
    endfor
  endfor
endfunction
 
function val = C2Tensionat(x0, x, y)
  % initializari
  n = length(x);
  a(1:n) = y(1:n);  
  h(1:n-1) = x(2:n) - x(1:n-1);
  
  % creare matrice tridiagonala
  ma(1:n-1) = h(1:n-1);

  mb(1) = 2 * h(1);
  mb(2:n-1) = 2 * (h(1:n-2) + h(2:n-1));
  mb(n) = 2 * h(n-1); 
  
  % calculare derivate in capete
  df1 = (y(2) - y(1)) / (x(2) - x(1));
  dfN =  (y(n) - y(n-1)) / (x(n) - x(n-1));
  
  g(1) =3 * (a(2) - a(1))/h(1) - 3 * df1;
  g(2:n-1) = 3 * (a(3:n) - a(2:n-1))/h(2:n-1) - 3 * (a(2:n-1) - a(1:n-2))/h(1:n-2);
  g(n) = 3 * dfN - 3 * (a(n)-a(n-1))/h(n-1);
 
  % rezolvare sistem
  c = Thomas(ma, mb, ma, g);
 
  % calculare coeficienti ramasi
  for i=1:n-1
    d(i) = (c(i+1) - c(i))/(3 * h(i));
    b(i) = (a(i+1) - a(i))/h(i) - h(i)/3 * (2*c(i) + c(i+1));
  endfor
  
  % calculare valori polinom
  nx = length(x0);
  for k = 1 : nx
    for i=1:n-1
      if x0(k) >= x(i) && x0(k) <= x(i+1)
        val(k) = a(i) + b(i) * (x0(k) - x(i)) + c(i) * (x0(k) - x(i))^2 + d(i)*(x0(k) - x(i))^3;
        break;
      endif
    endfor
  endfor
endfunction
 
function val = TFourier(m, x0, f)

  % initializare valori echidistante din interval
  for j = 1 : 2 * m - 1
    x(j) = -pi + ((j - 1)/m) * pi;
  endfor
  
  % calculare a si b
  for k = 1 : m
    suma = 0;
    for j = 1 : 2 * m - 1
      suma += f(x(j)) * cos((k - 1) * x(j));
    endfor 
    a(k) = suma / m;
    suma = 0;
    for j = 1 : 2 * m - 1
      suma += f(x(j)) * sin(k * x(j));
    endfor
    b(k) = suma / m;
  endfor
  suma = 0;
  for j = 1 : 2*m - 1
    suma += f(x(j)) * cos(m * x(j));
  endfor 
  a(m + 1) = suma / m;

  % calculare valori polinom trigonometric
  nx = length(x0);
  for i = 1 : nx
    suma = 0;
    for k = 1: m - 1
      suma += a(k + 1) * cos(k * x0(i)) + b(k) * sin(k * x0(i));
    endfor
    val(i) = (a(1) + a(m + 1) * cos(m * x0(i))) / 2 + suma;
  endfor
endfunction