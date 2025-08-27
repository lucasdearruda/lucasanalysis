import csv
import re
import sys

def parse_to_csv(input_file, output_file):
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
            
            # Adiciona os dados na lista
            parsed_data.append(columns)
    
    # Escreve os dados no arquivo CSV
    with open(output_file, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        
        # Escreve as linhas no CSV
        writer.writerows(parsed_data)

    print(f"Arquivo CSV gerado: {output_file}")

# Função principal
if __name__ == "__main__":
    # Verifica se os parâmetros de linha de comando foram fornecidos corretamente
    if len(sys.argv) != 3:
        print("Uso: python textparser.py <input_file.txt> <output_file.csv>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Chama a função para processar o arquivo e gerar o CSV
    parse_to_csv(input_file, output_file)
