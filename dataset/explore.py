import pandas as pd

def separate_detail_column():
    
    combined_data = pd.read_csv('/Users/ponynie/Developer/Python_Code/GCRI Prediction/dataset/data-combined.csv')

    # Count the occurrences of unique values in the "Details" column
    details_counts = combined_data["Details"].value_counts()

    # Print the counts
    print("Unique values in the 'Details' column and their counts:")
    for detail, count in details_counts.items():
        print(f"{detail}: {count}")

    df = combined_data

    # Define regular expressions to extract RI Type, Column Type, and Temperature Mode
    ri_type_pattern = r"([\w'\s]+) RI"
    column_type_pattern = r", (non-polar|polar) column"
    temperature_mode_pattern = r", (isothermal|temperature ramp|custom temperature program)"

    # Extract RI Type, Column Type, and Temperature Mode from the "Details" column
    df["RI Type"] = df["Details"].str.extract(ri_type_pattern, expand=False)
    df["Phase Polarity"] = df["Details"].str.extract(column_type_pattern, expand=False)
    df["Temperature Mode"] = df["Details"].str.extract(temperature_mode_pattern, expand=False)

    # Drop the "Details" column
    df = df.drop("Details", axis=1)

    # Save the modified DataFrame to a new CSV file
    df.to_csv('/Users/ponynie/Developer/Python_Code/GCRI Prediction/dataset/data.csv', index=False)

def qurries():
    # Read the CSV data
    df = pd.read_csv('/Users/ponynie/Developer/Python_Code/GCRI Prediction/dataset/data.csv')

    # Filter the DataFrame based on the specified conditions
    filtered_df = df[(df['RI Type'] == 'Van Den Dool and Kratz') & (df['Phase Polarity'] == 'polar')]

    # Aggregate the duplicates based on the specified columns
    aggregated_df = filtered_df.groupby(['Name', 'molecularFormula', 'inChI', 'inChIKey'])['I'].sum().reset_index()

    # Print the aggregated DataFrame
    aggregated_df.head(20)
    
    len(aggregated_df)