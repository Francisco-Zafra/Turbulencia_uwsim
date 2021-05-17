# Turbulencia_uwsim

Se ha añadido un suelo:
  - Profundidad del suelo: --floor-depth
  - Indice n: --floor-index (por defecto 1.5)

Particulares de la turbulencia:
  - Parametros de entrada actuales: 
    - Variacion en z de las transiciones: --varZ
    - Variacion del n del agua: --var_n_water
    - Maximo valor de theta proporcional a PI/2: --boundary_max_theta
    - Numero de capas: -L
  - Cambio de indices n y posicion de transiciones para cada bloque (10e5 fotones)
  - Valores aleatorios:
    - Theta: Uniforme
    - Phi: Uniforme
    - Variacion de n: Uniforme
    - Variacion de z: Uniforme
