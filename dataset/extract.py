from bs4 import BeautifulSoup
import os
import pandas as pd
import json

def extract_info(html_path: str) -> pd.DataFrame:
    with open(html_path, 'r', encoding='utf-8') as file:
        html_content = file.read()
    
    soup = BeautifulSoup(html_content, 'html.parser')
    tables = soup.find_all('table', class_='data')
    
    column_name = ("Column type", "Active phase", "I")
        
    if tables: #tables is not empty
        full_joined_table = pd.DataFrame()
        for table in tables:
            table_name = table.get('aria-label')
            table_rows = table.find_all('tr', class_='exp')
            agg_table_rows = pd.DataFrame()
            for table_row in table_rows:
                value_cell = table_row.find_all('td', class_='right-nowrap')
                value_cell = {column_name[i]: value_cell[i].text.strip() for i in range(len(value_cell))}
                value_cell['I'] = int(value_cell['I'][:-1])
                value_cell = pd.DataFrame(value_cell, index=[0])
                agg_table_rows = pd.concat([agg_table_rows, value_cell], ignore_index=True)
            agg_table_rows = agg_table_rows.groupby(['Column type', 'Active phase']).agg({'I': 'mean'}).reset_index()
            agg_table_rows['Details'] = table_name
            full_joined_table = pd.concat([full_joined_table, agg_table_rows], ignore_index=True)

        script_tags = soup.find_all('script', type='application/ld+json')
        if len(script_tags) >= 2:
            json_data = script_tags[1].string
            json_data = json.loads(json_data)
            full_joined_table['name'] = json_data['name']
            full_joined_table['molecularFormula'] = json_data['molecularFormula']
            full_joined_table['molecularWeight'] = json_data['molecularWeight']
            
            if 'monoisotopicMolecularWeight' not in json_data:
                full_joined_table['monoisotopicMolecularWeight'] = None
            else:
                full_joined_table['monoisotopicMolecularWeight'] = json_data['monoisotopicMolecularWeight']
                
            if 'inChI' not in json_data:
                full_joined_table['inChI'] = None
            else:
                full_joined_table['inChI'] = json_data['inChI']
                
            if 'inChIKey' not in json_data:
                full_joined_table['inChIKey'] = None
            else:
                full_joined_table['inChIKey'] = json_data['inChIKey']
                
            return full_joined_table
        else:
            return pd.DataFrame()
    else:
        return pd.DataFrame()   
    
    
def process_html_files(directory: str) -> pd.DataFrame:
    extracted_data = pd.DataFrame()
    
    for file_name in os.listdir(directory):
        if file_name.endswith('.html'):
            file_path = os.path.join(directory, file_name)
            extracted_info = extract_info(file_path)
            if extracted_info.empty == False:
                extracted_data = pd.concat([extracted_data, extracted_info], ignore_index=True)
    
    return extracted_data

def main():
    directory_path = 'dataset/testextract'
    extracted_data = process_html_files(directory_path)
    name = extracted_data.pop('name')
    extracted_data.insert(0, 'Name', name)
    extracted_data.to_csv('dataset/extracted_data.csv', index=False)
    
def test():
    print(extract_info('dataset/testextract/R2.html'))

main()


