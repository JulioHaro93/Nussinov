#!/usr/bin/env python
#_*_ coding: utf8 _*_


import random
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw
import time

with open('secuencia.txt', 'r') as archivo:
    contenido = archivo.read()

print(contenido)
def generar_cadena_genetica(longitud):
    bases = ['A', 'T', 'C', 'G']
    cadena = ""
    cadena += "ATG"
    for _ in range(longitud - 6):  # Restamos 6 para tener un múltiplo de 3
        cadena += random.choice(bases)

    cadena += random.choice(["TAA", "TAG", "TGA"])

    return list(cadena)

longitud = 30
cadena_genetica = generar_cadena_genetica(longitud)

print(f"Cadena genética generada: {cadena_genetica}")


def traduceCodigoaARN(cadena_genetica):
    i=0
    
    for letra in cadena_genetica:
        if(cadena_genetica[i]=='T'):
            cadena_genetica[i]='U'
            i+=1
        else:
            cadena_genetica[i]=cadena_genetica[i]
            i+=1
    cadenaARN=''.join(cadena_genetica)
    return cadenaARN


tiempoInicio = time.time()
def nussinov(secuencia):
    #Largo de la secuencia, esto considera las posiciones del arreglo, más no el tamaño del ARN, el tamaño del ARN es n-1

    n = len(secuencia)
    print(n)
    #Esto genera una matriz nxn llena de ceros que se inicia de esta manera para después ir sumando las coincidencias de nucleótidos apareados
    score = np.zeros((n, n), dtype=int)
    #genera un ciclo tomando en cuenta el tamaño del ARN
    for i in range(n - 1, -1, -1):
        for j in range(i + 1, n):
            noApareados = score[i + 1][j]
            #AQuí inicia el código del algoritmo de nussinov, el cuál va iterando en el rango i, j
            #i,j son los valores de los nucleótidos, j+1 es el siguiente nucleótido, de tal manera que se pueden ir comparando y llenando la matriz
            for k in range(i, j):
                #Esta validación permite que se comparen las bases que se pueden aparear
                if (secuencia[k], secuencia[j]) in [('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')]:
                    #Tenemos dos partes, si se aparean, entonces se suma + en la posición i de la matriz y en [k-1] mas [k+1] y [j-1]
                    apareados = 1 + score[i][k - 1] + score[k + 1][j - 1]
                    #y se agrega el valor de los no apareados con el valor máximo entre no apareados y apareados para continuar el bucle
                    noApareados = max(noApareados, apareados)
            #una vez realizada la suma, se establece el valor en la matriz en la posición [i][j] el valor máximo entre el valor[i+1] o el valor de los no apareados
            score[i][j] = max(score[i + 1][j], noApareados)
    valorOptimizacion = score[0][n-1]
    return score, valorOptimizacion

tiempoFinalNussinov = time.time()

tiempoTranscurrido = tiempoFinalNussinov-tiempoInicio
def visualizar_matriz(matriz,secuencia, nombre_archivo):
    altura, ancho = matriz.shape
    img = Image.new('RGB', (ancho * 20, altura * 20), color='white')
    draw = ImageDraw.Draw(img)

    # Dibujar los valores de la matriz
    for i in range(altura):
        for j in range(ancho):
            valor = str(matriz[i, j])
            x = j * 20
            y = i * 20
            draw.text((x, y), valor, fill='black')

    # Etiquetas de las bases nitrogenadas
    bases_nitrogenadas = ['A', 'U', 'G', 'C']
    for i, base in enumerate(secuencia):
        x = i * 20 + 10
        y = altura * 20 + 5
        draw.text((x, y), base, fill='black', font=None)

    # Guardar la imagen
    img.save(nombre_archivo)
    

def visualizar_matriz_con_colores(matriz, nombre_archivo):
    fig, ax = plt.subplots()
    cax = ax.matshow(matriz, cmap='viridis', origin='upper')
    fig.colorbar(cax)
    
    for i in range(matriz.shape[0]):
        for j in range(matriz.shape[1]):
            ax.text(j, i, str(matriz[i, j]), va='center', ha='center')

    plt.title("Matriz de Puntuación")
    plt.xlabel("Posición Final")
    plt.ylabel("Posición Inicial")
    plt.savefig(nombre_archivo)
    plt.close()

def puntoParentesis(matriz):

    n = len(matriz)
    estructura = ['.'] * n  # Inicializar la estructura con puntos

    def encontrar_parejas(i, j):
        nonlocal estructura
        if i < j and matriz[i][j] == matriz[i + 1][j]:
            estructura[i] = '.'
            estructura[j] = '.'
            encontrar_parejas(i + 1, j)
        elif i < j and matriz[i][j] == matriz[i][j - 1]:
            estructura[i] = '.'
            estructura[j] = '.'
            encontrar_parejas(i, j - 1)
        elif i < j:
            estructura[i] = '('
            estructura[j] = ')'
            encontrar_parejas(i + 1, j - 1)

    # Iniciar el proceso desde la esquina inferior izquierda
    encontrar_parejas(0, n - 1)

    return ''.join(estructura)


secuenciaARN = traduceCodigoaARN(cadena_genetica)
print("ARN:"+secuenciaARN)
result = nussinov(str(secuenciaARN))
print("Tiempo transcurrido para la finalización del algoritmo de Nussinov: \n"+str(tiempoTranscurrido))
visualizar_matriz_con_colores(result[0], "matriz_puntuacion_colores.jpg")
visualizar_matriz(result[0],secuenciaARN, "matriz_puntuacion.jpg")
print(f"La puntuación de la estructura secundaria óptima es: {result[1]}")

print("notación punto paréntesis de la matriz de Nussinov\n"+secuenciaARN +"\n"+puntoParentesis(result[0]))
