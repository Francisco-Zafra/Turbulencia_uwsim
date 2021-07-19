# Turbulencia_uwsim

Pantallas de fase:
  - Resolucion en mm (no es parametro, se cambia en el codigo)
  - Parametro --phase_json_x para introducir el archivo de derivadas parciales x
  - Parametro --phase_json_y para introducir el archivo de derivadas parciales y
  - Si no se indican alguno de los dos parametros, no se realizan pantallas de fase

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
