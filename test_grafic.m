function test_grafic(eps1, eps2)
  a = 3;
  f = @(x) exp(a * cos(x))/(2 * pi * besseli(0, a));
  N = 1000;
  x_tilda = linspace(-pi, pi, N + 1);
  y_tilda = f(x_tilda);
  
  %desenare primul subplot
  subplot (2, 1, 1)
  [n x1 y1] = eval_interpolator_c(1, eps1);
  [n x2 y2] = eval_interpolator_c(2, eps1);
  [n x3 y3] = eval_interpolator_c(3, eps1);
  [n x4 y4] = eval_interpolator_c(4, eps1);
  [n x5 y5] = eval_interpolator_c(5, eps1);
  [n x6 y6] = eval_interpolator_c(6, eps1);
  plot (x_tilda, y_tilda, 'k-', x1, y1, 'r-x', x2, y2, 'g-', x3, y3, 'b-*', x4, y4, 'y-o', x5, y5, 'm-', x6, y6, 'c-');
  legend("f", "lagrange", "newton", "linear spline", "natural", "cubic spline", "fourier");
  
  %desenare al doilea subplot
  subplot (2, 1, 2)
  data = load("sunspot.dat");
  x_tilda = data(:, 1)';
  y_tilda = data(:, 2)';
  [n x1 y1] = eval_interpolator_d(1, eps2);
  [n x2 y2] = eval_interpolator_d(2, eps2);
  [n x3 y3] = eval_interpolator_d(3, eps2);
  [n x4 y4] = eval_interpolator_d(4, eps2);
  [n x5 y5] = eval_interpolator_d(5, eps2);
  [n x6 y6] = eval_interpolator_d(6, eps2);
  plot (x_tilda, y_tilda, 'k-', x1, y1, 'r-x', x2, y2, 'g-', x3, y3, 'b-*', x4, y4, 'y-o', x5, y5, 'm-', x6, y6, 'c-');
  legend("fi", "lagrange", "newton", "linear spline", "natural", "cubic spline", "fourier");
  
endfunction