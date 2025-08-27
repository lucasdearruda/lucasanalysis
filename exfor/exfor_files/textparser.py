import csv
import re
import sys

def parse_to_csv(input_file, output_file, extra_column=None):
    with open(input_file, 'r') as file:
        lines = file.readlines()
    
    # Lista para armazenar os dados processados
    parsed_data = []
    
    # Itera pelas linhas do arquivo
    for line in lines:
        # Remove os espaços no começo e no final
        line = line.strip()

        # Verifica se a linha começa com número (é uma linha de dados válida)
        if re.match(r'^\d', line):
            # Divide a linha com base nos espaços e remove qualquer espaço extra
            columns = re.split(r'\s+', line)
            
            # Se o argumento extra_column foi fornecido, adiciona como primeira coluna
            if extra_column is not None:
                columns.insert(0, str(extra_column))
            
            # Adiciona os dados na lista
            parsed_data.append(columns)
    
    # Escreve os dados no arquivo CSV
    with open(output_file, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        
         # Escreve o cabeçalho no CSV
        header = ["#EN", "E", "ANG", "DATA", "ERR-S"]
        writer.writerow(header)

        # Escreve as linhas no CSV
        writer.writerows(parsed_data)

    print(f"Arquivo CSV gerado: {output_file}")

# Função principal
if __name__ == "__main__":
    # Verifica se os parâmetros de linha de comando foram fornecidos corretamente
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Uso: python textparser.py <input_file.txt> <output_file.csv> [<extra_column_value>]")
        print("\nO parâmetro extra_column na sua função parse_to_csv serve para inserir um valor extra como uma nova coluna no início de cada linha de dados extraída do arquivo de entrada.")
        print("Exemplo: python textparser.py input.txt output.csv 42.0, vai inserir o valor 42.0 como a primeira coluna em cada linha do CSV gerado.")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Verifica se um valor extra foi fornecido como argumento
    extra_column = None
    if len(sys.argv) == 4:
        try:
            extra_column = float(sys.argv[3])  # Converte para float se possível
        except ValueError:
            print("Erro: o valor extra fornecido não é um número válido.")
            sys.exit(1)

    # Chama a função para processar o arquivo e gerar o CSV
    parse_to_csv(input_file, output_file, extra_column)
