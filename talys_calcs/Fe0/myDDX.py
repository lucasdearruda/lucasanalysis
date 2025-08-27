#!/usr/bin/env python3
import sys, os, re
import pandas as pd
import matplotlib.pyplot as plt

# Mapeamento do final do arquivo para isótopo
deg_to_isotope = {
    '054': '54Fe',
    '056': '56Fe',
    '057': '57Fe',
    '058': '58Fe'
}
def parse_isotopes(output_file):
    """Lê o output do TALYS e retorna um dicionário {isotope: fraction}"""
    isotopes = {}
    if not os.path.exists(output_file):
        print(f"Error: output file '{output_file}' not found")
        return isotopes

    with open(output_file, 'r') as f:
        lines = f.readlines()

    start = False
    for line in lines:
        line_strip = line.strip()
        # inicia leitura após encontrar "multi-isotope case"
        if 'multi-isotope case' in line_strip:
            start = True
            continue
        if start:
            # ignora linhas vazias ou com cabeçalhos
            if line_strip == '' or 'Isotope Abundance' in line_strip:
                continue
            # captura linhas do tipo "54Fe   0.058450"
            match = re.match(r'(\d+Fe)\s+([0-9\.Ee+-]+)', line_strip)
            if match:
                isotopes[match[1]] = float(match[2])
            else:
                # para no fim da tabela
                break
    return isotopes

def load_pddx_file(filename):
    """Lê um arquivo pddx em DataFrame"""
    df = pd.read_csv(filename, sep=r'\s+', comment='#',
                 names=['E', 'xs', 'Direct', 'Preeq', 'Multiple_preeq', 'Compound'])
   
    return df
def main():
    if len(sys.argv) != 4:
        print("Usage: checkDDX.py <angle> <particle> <ENN>")
        sys.exit(1)

    angle, particle, ENN = sys.argv[1], sys.argv[2], sys.argv[3]
    output_file = 'output'  # arquivo TALYS principal

    # Passo 1: ler abundâncias isotópicas
    isotopes = parse_isotopes(output_file)
    if not isotopes:
        print("No isotopes found. Check the output file format!")
        sys.exit(1)
    print("Isotopes found:", isotopes)

    # Passo 2: carregar arquivos de cada isótopo e ponderar
    df_total = None
    for deg_code, isotope in deg_to_isotope.items():
        if isotope not in isotopes:
            continue
        frac = isotopes[isotope]
        ENN_float = float(ENN)
        ENN_str = f"{ENN_float:08.3f}"  # 0025.000
        angle_str = f"{int(angle):03d}.0"  # 020.0
        filename_pattern = f"{particle}ddxE{ENN_str}A{angle_str}.deg.{deg_code}"
        if not os.path.exists(filename_pattern):
            print(f"Warning: file {filename_pattern} not found, skipping.")
            continue
        df_iso = load_pddx_file(filename_pattern)
        df_iso.iloc[:, 1:] *= frac  # pondera todas as colunas exceto energia
        if df_total is None:
            df_total = df_iso.copy()
        else:
            df_total.iloc[:, 1:] += df_iso.iloc[:, 1:]

    if df_total is None:
        print("No data loaded. Exiting.")
        sys.exit(1)

    # Passo 3: plot
    plt.figure(figsize=(8,5))
    plt.plot(df_total['E'], df_total['xs'], label='Etot')
    plt.plot(df_total['E'], df_total['Direct'], '--', label='Direct')
    plt.plot(df_total['E'], df_total['Preeq'], '--', label='Preequilibrium')
    plt.plot(df_total['E'], df_total['Multiple_preeq'], '--', label='Multiple Preeq')
    plt.plot(df_total['E'], df_total['Compound'], '--', label='Compound')
    plt.xlabel('E [MeV]')
    plt.ylabel('d²σ/dΩdE [mb/MeV/sr]')
    plt.title(f'{particle}-emission at {ENN} MeV, {angle} deg (Fe natural)')
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()

if __name__ == "__main__":
    main()

