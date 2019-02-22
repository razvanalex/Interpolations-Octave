function val = test(eps1, eps2)
  % calculare N pentru fiecare interpolare
  for i = 1:6
    val(1, i) = eval_interpolator_c(i, eps1);
    val(2, i) = eval_interpolator_d(i, eps2);
  endfor
  
  % afisare matrice
  disp(val);
endfunction