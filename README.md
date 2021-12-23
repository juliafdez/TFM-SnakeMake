# TFM-SnakeMake
workflow de análisis bioinformático para datos de metilación generados con la tecnología de BS-seq

Este proyecto es el trabajo final del master en bioinformática y bioestadística. Impartido por la Universidad Oberta de Catalunya y la Universidad de Barcelona.

Consta de un workflow realizado con Snakemake para obtener información sobre metilación. El workflow consta de los siguientes pasos y tecnologías:

1. Control de calidad de las secuencias - FastQC
2. Trimming de las secuencias - Trim Galore!
3. Control de calidad del trimming - FastQC
4. Comparación de calidad de las secuencias antes y después del trimming - MultiQC
5. Preparación del genoma de referencia - Bismark
6. Alineamiento contra el genoma de referencia - Bismark
7. Deduplicación - Bismark
8. Ordenar (opcional. Inactivo en el workflow) - Samtools
9. Extraccion de metilación - Bismark
10. Reporte del alineamiento y la metilación - Bismark
11. Extracción de los DML/DMR - Paquete DSS (R)

Para poder realizaar la ejecución del workflow se necesita realizar la instalacion de conda.
Una vez instalado, instalar mediante conda los paquetes fastqc, multiqc y bismark (el cual contiene bismark_align, bismark_methylation_extractor, bismark2report)
