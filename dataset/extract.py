from bs4 import BeautifulSoup
import os
import pandas as pd
import json

def extract_info(html_file):
    with open(html_file, 'r', encoding='utf-8') as file:
        html_content = file.read()
    
    soup = BeautifulSoup(html_content, 'html.parser')
    table_row = soup.find('tr', class_='exp')
    
    if table_row == None:
        return None
    else:
        script_tags = soup.find_all('script', type='application/ld+json')
        if len(script_tags) >= 2:
            json_data = script_tags[1].string
            json_data = json.loads(json_data)
            json_data.pop('@context')
            json_data.pop('@type')
            
            value_cell = table_row.find_all('td', class_='right-nowrap')
            column_name = ("Column type", "Active phase", "I")
            value_cell = {column_name[i]: value_cell[i].text.strip() for i in range(len(value_cell))}
            
            json_data.update(value_cell)
            
            return json_data
        
        else:
            return None
        
    
    
def process_html_files(directory):
    extracted_data = []
    
    for file_name in os.listdir(directory):
        if file_name.endswith('.html'):
            file_path = os.path.join(directory, file_name)
            extracted_info = extract_info(file_path)
            if extracted_info:
                extracted_data.append(extracted_info)
    
    return extracted_data

def main():
    directory_path = 'dataset/testextract'
    extracted_data = pd.DataFrame(process_html_files(directory_path))
    extracted_data['I'] = extracted_data['I'].apply(lambda x: int(x[:-1]))
    extracted_data.to_csv('dataset/extracted_data.csv', index=False)
    
def test():
    extract_info('dataset/testextract/TE1.html')

main()