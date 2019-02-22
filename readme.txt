---------------------------------------------------------------------
                 Tema 3 - Metode Numerice - Interpolari
---------------------------------------------------------------------

Autor:
  Smadu Razvan-Alexandru  315CB

Fisiere incluse:
  - eval_interpolator_c.m
  - eval_interpolator_d.m
  - test_grafic.m
  - test.m
  - Thomas.m (algoritmul din laborator)
  - readme.txt
  
README
  1. Evaluare continua
   1.1 Scurta descriere a implementarii
    In fisierul eval_interpolator_c.m au fost implementati cei 6 algoritmi 
    de interpolare (Lagrange, Newton, spline liniar, spline cubic natural,
    spline cubic tensionat si Fourier) care sunt aplicati asupra functiei
    de densitate de energie emisa in timpul maximului unui ciclu de 11 ani
    de soare. Functia este impartita practic in 3 sectiuni: determinarea 
    celor 1000 de puncte, calcularea erorii E si verificarea daca se 
    indeplinesc conditiile de convergenta. De asemenea sunt implementate si 
    functii pentru fiecare interpolare (unele luate din laborator altele 
    scrise pe baza informatiilor din curs), dintre care unele au fost adaptate
    pentru a primi ca intrare, pe langa nodurile de interpolare, un vector 
    de abcise in care sa se calculeze valorile, pentru a scadea timpul de 
    executie pentru unii algoritmi. Pentru interpolarea cu spline cubic 
    tensionat si natural, rezolvarea sistemul se face cu ajutorul algoritmului
    Thomas, gasit in fisierul cu acelasi nume.
    
   1.2 Interpretarea rezultatelor obtinute
    Pentru eps=1 toate metodele converg, insa aproximarea nu este cea mai buna.
    Cea mai buna aproximare pare sa o dea interpolarea trigonometrica (Fourier),
    cele trei interpolari prin spline-uri dau aproximativ acelasi rezultat 
    (arata ca functia de interpolare spline liniar), iar Lagrange si Newton 
    dau tot cam acelasi rezultat impreuna, dar aproximarea nu este foarte buna
    la capete unde se observa fenomenul Runge. Toti algoritmii converg la 
    fel de repede.
    
    Pentru eps=0.1 la interpolarile Newton si Lagrange se observa la capete 
    fenomenul Runge, dar in interiorul domeniului aproximarea este foarte
    buna. Insa aproximari mai bune, din punct de vedere al intregului domeniu,
    le confera celelalte interpolari. Inca o data interpolarea Fourier are
    cel mai bun rezultat (aproximarea se realizeaza la o precizide de ~0.00001).
    La interpolarile cu spline-uri sa obtinut cam acelasi rezultat. 
    Interpolarile Newton si Lagrange nu converg, iar viteza cea mai mare de 
    convergenta o are interpolarea Fourier.
    
    Pentru eps=0.001 polinoamele de interpolare sunt foarte apropiate de functia
    reala. La interpolarile Newton si Lagrange din nou se observa fenomenul 
    Runge la capetele intervalului. Cea mai buna aproximare este data 
    de interpolarea trigonometrica. Interpolarile Newton si Lagrange nu converg,
    iar viteza cea mai mare de convergenta o are interpolarea Fourier.
    
   1.3 Concluzii trase in urma obtinerii acestor rezultate
    Cel mai bun polinom de interpolare pentru aceasta functie este Fourier,
    care are atat o viteza mare de convergenta cat si ofera o aproximare foare
    buna. Polinomul Lagrange si polinomul Newton ofera o aproximare buna in 
    interiorul intervalului, insa la capete, din cauza fenomenului Runge, 
    aproximarea nu e cea mai buna. Celelalte polinoame nu ofera o aproximare 
    chiar asa de buna fata de cele mentionate pana acum, dar pastreaza destul 
    de bine forma functiei, cu cat eps este mai aproape de 0.
   
  2. Evaluare discreta
    Avand in vedere ca in acest caz functia continua a fost inlocuita cu o
    functie discreta, au aparut neconcordante cu formulele utilizate.
    La functia continua se putea cunoaste functia in orice punct; la
    aceasta functie discreta, sunt cunoscute un numar finit de puncte.
    In concluze a trebuit sa aproximez cumva si valorile dintre puncte (pentru
    a putea aplica aceeasi abordare ca cea de la evaluarea contiuna). Solutia
    pe care am aplicat-o este de a considera ca in fiecare an, numarul de pete
    solare este acelasi, iar modificarea acestui numar sa se faca doar
    la inceputul anului. Cu alte cuvinte:
    
                                 +-
                                 | f(x1),     x1 din [1700, 1701)
                         f(x) = <  f(x2),     x2 din [1701, 1702)
                                 | ...
                                 | f(x300), x300 din [2000, 2001)
                                 +-
    
    Datoria acestei aproximari, am fost nevoit sa limitez valoarea maxima a lui
    Nk intrucat, pentru valori mai mari de 300, aproximarea nu ar mai fi mai
    buna decat in cazul pt 300. O alta modificare a aparut si in cazul 
    interpolarii trigonometrice in care functia TFourier primea ca argument 
    o functie; in cazul discret primeste doar nodurile de interpolare.
    Functia primeste un vector a carei dimensiune este un numar impar de forma
    2 * n - 1 pentru valorile lui x. De asemenea, a trebuit sa fie modificata si
    modalitatea de evaluare a erorii, intrucat in metoda de la varianta discreta 
    variatia erorii nu mai este descrescatoare si sa tinda la 0. Atunci am 
    hotarat sa folosesc teorema de aproximare a lui Weierstass care spune ca 
    
                    lim(     max  |f(x) - P (x)| ) = 0
                   n->inf  a<=x<=b         n
                   
    Eroarea pe care o calculez este termenul din limita.
    De asemenea, in urma unor observatii din teste manuale, in cazul de 
    convergenta, eroarea este foare aproape de 0, iar E(Nk) = E(Nk+1).
    Pentru a diferentia cazurile de divergenta cu cele de convergenta
    dar pentru un numar mai mare de puncte, am limitat eroarea la cubul
    celei mai mari valori ale functiei (intrucat in caz de divergenta,
    erorile tind monoton la infinit). In caz de convergenta nu exista 
    o anume monotonie, dar valorile tind la 0. 
    
    Nota: Aplicand un caz continuu (functii concrete precum x^2) pe cazul
    discret (lunad un numar finit de puncte) s-a obtinut acelasi rezulat
    ca in cazul continuu.
    
  3. Interpretarea graficelor pentru functia discreta si convergenta lor
    Pentru testul test_grafic(0.01, 0) (am ales aceste valori pentru ca ele
    releva o aproximare destul de buna pentru unii algoritmi de interpolare),
    la interpolarea discreta, cele mai bune rezultate sunt obtinute prin 
    interpolarea folosind splile-uri. Interpolarea Newton si interpolarea
    Lagrange nu aproximeaza deloc bine functia discreta, iar interpolarea
    Fourier da o functie (aparent periodica) dar nu este de aceeasi "faza" si
    "amplitudine" cu functia discreta data. 
    Ruland test(0.01, 0) se poate observa ca, precum in cazul continuu, in 
    cazul discret, interpolarile Newton si Lagrange nu converg. La 
    interpolarile folosind spline, cele 3 metode converg la fel. In cazul
    Fourier, se observa ca functia discreta nu converge (se poate observa
    si din grafic acest lucru).
    
    Matricea obtinuta pentru acest test:
               _                                   _
              |  Inf   Inf    32    32    32    15  |
              |_ Inf   Inf   300   300   300   Inf _|
    
    Pentru testul test(1, 250) se observa ca algoritmii converg, insa nu
    au si precizia dorita. Niciuna dintre metode nu aproximeaza bine functia
    discreta.
    
    Matricea obtinuta pentru acest test: 
               _                       _
              |  4   4   4   4   4   3  |
              |_ 4   4   4   4   4   5 _|    
    
       
  4. Observatii si concluzii
    - Diferentele divizate au fost calculate eficient din punct de vedere al 
    memoriei, folosind un singur vector
    - La functiile C2Natural() si C2Tensionat() ma, mb si respectiv mc  
    reprezinta valorile diagonalelor matricelor prin care se calculeaza c.
    - Apelul functiei din test.m se face test(eps1, eps2) unde eps1 si eps2
    reprezinta erorile dorite (pentru eps1 se recomanda valori sub 1; pentru
    eps2 valori intre 0 si 300).
    - Apelul functiei din test_grafic.m se face la ca test.m 
    (test_grafic(eps1, eps2)) unde eps1 si eps2 au acelasi rol.
    - Interpolarile Newton si Lagrage nu sunt relevante pentru un sir discret
    de puncte precum cele din fisierul sunspot.dat
    - Am folosti mai multe functii si mai multe fisiere, intrucat pe forumul
    asociat temei s-a specificat ca este permis acest lucru
    
  Bibilografie pentru aspectele teoretice:
  [1] https://en.wikipedia.org/wiki/Stone%E2%80%93Weierstrass_theorem 
  [2] https://en.wikipedia.org/wiki/Runge%27s_phenomenon

