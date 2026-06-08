import pandas as pd
import argparse

def process_csv(input_file, output_file, window_size, step_size):
    # Leer archivo CSV delimitado por tabulador
    df = pd.read_csv(input_file, sep='\t')

    # Agrupar por cromosoma
    grouped = df.groupby('Chromosome')

    results = []

    for chrom, group_df in grouped:
        # Ordenar por posición
        group_df = group_df.sort_values(by='Position')

        # Empezar siempre en 1 aunque no esté en los datos
        min_pos = 1
        max_pos = group_df['Position'].max()

        # Iterar ventanas desde min_pos hasta max_pos
        for start in range(min_pos, max_pos + 1, step_size):
            end = start + window_size - 1

            # Filtrar posiciones dentro de la ventana
            window = group_df[(group_df['Position'] >= start) & (group_df['Position'] <= end)]

            if window.empty:
                print(f"Bin vacío omitido: {chrom} {start}-{end}")
                continue

            # Contar ceros por individuo
            counts = window.iloc[:, 3:].eq(0).sum()

            # Cantidad de sitios en la ventana
            n_sites = window.shape[0]

            # Guardar resultado
            result_row = [chrom, start, end, n_sites] + list(counts)
            results.append(result_row)

    # Crear DataFrame con los resultados
    result_df = pd.DataFrame(
        results,
        columns=['Chromosome', 'Start', 'End', 'Sites_in_bin'] + list(df.columns[3:])
    )

    # Guardar en archivo CSV
    result_df.to_csv(output_file, index=False, header=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Agrupa datos en bins por posición real del VCF.')
    parser.add_argument('input_file', help='Archivo CSV de entrada')
    parser.add_argument('output_file', help='Archivo CSV de salida')
    parser.add_argument('--window_size', type=int, default=100, help='Tamaño de la ventana (bp)')
    parser.add_argument('--step_size', type=int, default=100, help='Paso entre ventanas (bp)')
    args = parser.parse_args()

    process_csv(args.input_file, args.output_file, args.window_size, args.step_size)
