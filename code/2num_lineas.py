#!/usr/bin/env python
import pandas as pd
import numpy as np

vector = open("results/visualizar_datos/temporal.txt", "r")

mi_lista = []

for linea in vector:
    mi_lista.append(len(linea) -1)
    
array = np.array(mi_lista)

maximo = array.max()
minimo = array.min()

for valor in [maximo, minimo]:
    print(valor)
