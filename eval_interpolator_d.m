function [N x_r p_r] = eval_interpolator_d(tip, eps)  
  data = load("sunspot.dat");
  x_tilda = data(:, 1)';
  y_tilda = data(:, 2)';
  N = length(x_tilda);
 
  % calculul E(p(x), Nk)
  k = 2;
  Nretur = k ^ 2;
  while (1)
    h = 2 * pi / (N + 1);
    Nk = 2^k;
    if Nk > N
      Nk = N;
    endif
    xk = linspace(x_tilda(1), x_tilda(N), Nk);
    delta_x = (N - 1)/(Nk - 1);
    for i = 1 : Nk
      yk(i) = y_tilda(ceil((i - 1) * delta_x) + 1);
    endfor
    switch(tip)
      case 1 % interpolare Lagrange
        pNk = PolinomLagrange(x_tilda, xk, yk);
      case 2 % interpolare Newton
        pNk = PolinomNewton(x_tilda, xk, yk);
      case 3 % linear spline
        for i = 1 : N
          pNk(i) = InterpolareLiniara(x_tilda(i), xk, yk);
        endfor
      case 4 % spline cubic natural
        pNk = C2Natural(x_tilda, xk, yk);
      case 5 % spline cubic tensionat
        pNk = C2Tensionat(x_tilda, xk, yk);
      case 6 % Fourier 
        xk = linspace(x_tilda(1), x_tilda(N), Nk + 1);
        delta_x = (N - 1)/(Nk);
        yk = 0;
        for i = 1 : Nk + 1
          yk(i) = y_tilda(ceil((i - 1) * delta_x) + 1);
        endfor
        pNk = TFourier(x_tilda, xk, yk);
        Nk++;
    endswitch    
    
    % Caluclare eroare
    E(1) = max(abs(y_tilda - pNk));
    if E(1) == 0 % se atinge cel mai probabil pt valori extrem de mari
      N = inf;
      return
    endif
    if length(E) == 2 
      if (E(1) > max(yk)^3) % valorile cresc la infinit
        N = inf;
        return;
      endif
      if (E(1) - E(2)) == 0 
        if E(1) < 1  % valorile tind la 0
          N = Nretur;
          x_r = x_tilda;
          p_r = pNk;
        else N = inf; % valorile cresc la infinit
        endif
        return;
      endif
      if norm(E(1) - E(2)) < eps % oprire la eroarea dorita
       N = Nretur;
       x_r = x_tilda;
       p_r = pNk;
       return;
     endif
    endif
    
    % trecere la pasul urmator
    E(2) = E(1);
    x_r = x_tilda;
    p_r = pNk;
    k++;
    Nretur = Nk; 
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

function val = TFourier(x0, x, y)
  m = (length(x) - 1) / 2;
  
   % calculare a si b 
  for k = 1 : m
    suma = 0;
    for j = 1 : length(x)
      suma +=y(j) * cos((k - 1) * x(j));
    endfor 
    a(k) = suma / m;
    suma = 0;
    for j = 1 : length(x)
      suma +=y(j) * sin(k * x(j));
    endfor
    b(k) = suma / m;
  endfor
  suma = 0;
  for j = 1 : length(x)
    suma += y(j) * cos(m * x(j));
  endfor 
  a(m + 1) = suma / m;
  
  % calculare valori polinom trigonometric
  nx = length(x0);
  for i = 1 : nx
    suma = 0;
    for k = 1: m
      suma += a(k + 1) * cos(k * x0(i)) + b(k) * sin(k * x0(i));
    endfor
    val(i) = (a(1) + a(m + 1) * cos(m * x0(i))) / 2 + suma;
  endfor
endfunction