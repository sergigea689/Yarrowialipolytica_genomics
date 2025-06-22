# Script for the analysis with fast_LMM

# Importar el algoritmo
import numpy as np
from fastlmm.association import single_snp
from fastlmm.util import example_file
import os
import pylab
import fastlmm.util.util as flutil
import pandas as pd
from fastlmm.util.stats import plotp

# Rutas a los archivos
bed_fn = 'fastlmm_input'  # Basename of input snps in plink format
pheno_fn = 'pheno1.txt'     # Archivo con los fenotipos
cov_fn = 'covariables_final.txt' 

# Verificar si los archivos de salida ya existen
if os.path.exists('manhattan_plot.png') and os.path.exists('qq_plot.png') and os.path.exists('results_df.csv'):
    print("Los archivos ya existen. El script ya ha sido ejecutado anteriormente.")
else:
    # Ejecutar FaST-LMM solo si los archivos no existen
    print("Ejecutando el análisis...")
    results_df = single_snp(bed_fn, pheno_fn, covar=cov_fn,  count_A1=False)
    
    # Verificando el DataFrame de resultados
    print("Análisis completado, comenzando la representación de gráficos...")

    # Graficar el gráfico de Manhattan
    pylab.rcParams['figure.figsize'] = (10.0, 8.0)  # Configuración del tamaño de la figura
    flutil.manhattan_plot(results_df[["Chr", "ChrPos", "PValue"]].values, pvalue_line=1e-5, xaxis_unit_bp=False)

    # Guardar el gráfico de Manhattan
    print("Guardando gráfico de Manhattan...")
    pylab.savefig('manhattan_plot.png')
    pylab.show()

    # Graficar el gráfico QQ
    print("Generando gráfico QQ...")
    plotp.qqplot(results_df["PValue"].values, xlim=[0, 5], ylim=[0, 5])

    # Guardar el gráfico QQ
    print("Guardando gráfico QQ...")
    pylab.savefig('qq_plot.png')

    # Mostrar el gráfico
    pylab.show()

    # Guardar los resultados en un archivo CSV
    print("Guardando resultados en CSV...")
    results_df.to_csv('results_df.csv', index=False)

    print("Análisis completado y archivos guardados.")
