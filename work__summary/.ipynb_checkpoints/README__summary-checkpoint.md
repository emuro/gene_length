# Resumen en 3 bloques: 
- Genes
- Proteínas 
- Transición de fase:

## 1. GENES

([1](./c1_1_gene_lognormDist.ipynb)) La distribución que mejor ajusta a las longitudes de los genes codificantes de proteínas de las **33.630** entradas/assemblies de ENSEMBL ($33,.630 \supset 21,922$ taxa diferentes) que analizamos es la lognormal, (Fig 1). En la literatura no hemos encontrado interés por estas distribuciones y apenas hay reportado nada.

([2](./c1_2_gene_taylor.ipynb)) Cuando representamos la longitud media Lg frente a la varianza de las 33.460 especies($33.630 \supset 33,460 $ con taxa además validada por NCBI) vemos que existe una relación potencial (una ley de Taylor, Figura 3).

_Nota: Mus_caroli es eliminado también porque su varianza tiene un valor tan elevado que está fuera de rango._

(3) Además sobre esta recta en log-log, los clados se ordenan según su aparición evolutiva: bacteria, archaea, fungi, plants, metazoans, vertebrata (los protistas aparecen desparramaos como el cajón de sastre que son). Si miramos con algo más de detalle en un clado como el de vertebrata, vemos que se cumple también para la secuencia peces, pájaros y primates.

(4) Así, que en general, a mayor longitud media de genes o varianza en la distribución de una especie más reciente es evolutivamente hablando. De donde podemos deducir que tanto la media como la varianza de las distribuciones han crecido a lo largo del tiempo, desde las medias y varianzas más pequeñas de las bacterias y arqueas hasta las de los vertebrados.

(6) De lo anterior inferimos que el proceso dinámico que ha generado las distribuciones de las longitudes de los genes ha sido capaz de mantener la lognormalidad mientras crecían en media y varianza estableciendo entre ambos parámetros una ley potencial que se ha mantenido a lo largo de toda la evolución.

(7) Atendiendo a lo que conocemos sobre el crecimiento de los genes proponemos un sencillo modelo multiplicativo que da cuenta tanto de esta estabilidad lognormal como de la ley de Taylor.

## 2. PROTEÍNAS

([1](./c2_1_prot_lognormDist.ipynb)) De los (19.854) proteomas de referencia de Uniprot utilizamos las **9.915** ($\subset 19.854$) especies que corresponden a los clados analizados en los apartados anteriores: bacteria, archaea, protists, fungi, plants, metazoans y vertebrata. Es decir, no tenemos en cuenta los 9915 proteomas de virus.

De nuevo encontramos para las **9.915** especies, sin contar los virus, que las distribuciones (longitud de sus proteínas) son lognormales (Figura 1). Se habían ajustado solo unas pocas en la literatura reportando lognormalidad también.

([2.a](./c2_2a_prot_taylor.ipynb)) De nuevo se cumple la misma ley de Taylor (Figura 4).

([2.b](./c2_2b_prot_taylor.ipynb)) Igualmente, cuando consideramos las **6.521** especies resultantes después de hacer una asociación de especies con las de ensembl ($ 6.521 \subset 9.915$).  
Nota: A pesar de filtrar tantas species con diferentes identificadores taxonómicos (3.394 menos ahora), la línea de regresión resultante es prácticamente la misma.

(3) La longitud media de las proteínas $L_p$ y la varianza ($\sigma_p^2$) siguen clusterizando (clados), pero ahora no parecen crecer siguiendo el orden evolutivo, aunque sí permiten observar el sistema de tres dominios con cierto solapamiento (bacterias, arqueas, eukariotas) e incluso los principales clados de eukariotas, con solapamiento, pero el orden evolutivo no se mantiene. 

(4) Todo esto sugiere que la longitud media de los genes funciona mejor como proxy de _"complejidad"_, como han propuesto algunos autores, que la longitud media de las proteínas, como han apuntado otros.


## 3. LONGITUD MEDIA UMBRAL

([1](./c3_1_longitud_media_umbral.ipynb)) Si representamos $L_g$ frente a $L_p$, para las **6521** especies de las que tenemos datos de genes y proteínas cruzadas (Figura 6), vemos que hasta un valor umbral $L_g = L_c$, la longitud media de los genes es 3 veces más grande que la de las proteínas (y si representamos las varianzas hasta un valor $\sigma^2_g=\sigma^2_c$, la varianza de los genes es 9 veces mayor). Las distribuciones están simplemente escaladas.

(2) Eso no ocurre a partir del valor umbral Lc (o Var_c), donde la longitud media de las proteínas se mantiene con fluctuaciones que van amortiguándose alrededor del valor Lc/3 (o Var_c/9). Estas fluctuaciones explican por qué en la gráfica de la ley de Taylor de las proteínas a partir de Lc/3 no colapsan todas las especies en un punto y no siguen el orden evolutivo.

(3) Estimamos Lc y a partir de su valor y el modelo multiplicativo (con dos puntos de referencia) concluimos que corresponde aproximadamente al momento en el que aparece la célula eucariota en la evolución.

(4) Efectivamente, Lc separa bastante bien los datos en procariotas y eucariotas.

(5) Para procariotas (bacterias y arqueas) cuyos genes están constituidos solo por exones, simplemente se traducen 3 a 1 las bases a aminoácidos. Pero para eucariotas, a pesar de que sigue creciendo en tamaño medio de genes, no crece el tamaño medio de las proteínas.

(6) Concluimos que a partir de este tamaño o momento el mecanismo de crecimiento de los genes sigue siendo el mismo, pero incorpora casi exclusivamente intrones que no se traducen en proteínas.


# 4. TRANSICIÓN DE FASE

(1) Si representamos los datos Lg vs rho (la densidad de intrones) la gráfica (Figura 7) se asemeja a una transición de fase de segundo orden, con Lg actuando de parámetro de control y Lc de punto crítico separando dos fases: la exónica donde el parámetro de orden es nulo (ausencia de exones) y la intrónica donde el parámetro de orden crece con Lg (y solo se incorporan intrones).

(3) Como la varianza y la entropía de rho muestran cualitativamente máximos alrededor del punto crítico Lc, pero resultan medidas ruidosas, proponemos una medida de diversidad de interpretación semejante, que claramente se hace máxima en el punto crítico como se esperaría de una transición de fase.

(2) Esta transición es de naturaleza computacional, describe el paso de organismos procariotas donde los genes se transcriben en su totalidad a proteínas formando una red genética autorregulada, a organismos eucariotas y muticelulares donde se requiere un salto cualitativo de capacidad de cómputo para regular redes de redes.

(3) Los intrones capacitaron a los genes para orquestar la complejidad que entrañan los organismos eucariotas y pluricelulares.
