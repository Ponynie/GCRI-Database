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

def queries():
    # Read the CSV data
    df = pd.read_csv('/Users/ponynie/Developer/Python_Code/GCRI Prediction/dataset/data.csv')

    # Filter the DataFrame based on the specified conditions
    filtered_df = df[(df['RI Type'] == "Kovats'") & (df['Phase Polarity'] == "non-polar")]

    # Aggregate the duplicates based on the specified columns
    aggregated_df = filtered_df.groupby(['Name', 'molecularFormula', 'inChI', 'inChIKey'])['I'].mean().reset_index()

    # Print the aggregated DataFrame
    aggregated_df.head(20)
    
    len(aggregated_df)

def queries(**kwargs):
    """
    Filters and aggregates the CSV data based on the specified conditions.

    Args:
        **kwargs: Key-value pairs for filtering the DataFrame.

    Returns:
        pandas.DataFrame: The aggregated DataFrame.
    """
    # Read the CSV data
    data_path = '/Users/ponynie/Developer/Python_Code/GCRI Prediction/dataset/data.csv'
    df = pd.read_csv(data_path)

    # Create a dictionary mapping keyword arguments to column names
    column_mapping = {
        'RI_Type': 'RI Type',
        'Phase_Polarity': 'Phase Polarity',
        'Temperature_Mode': 'Temperature Mode',
        'Active_Phase': 'Active phase',
        'Column_Type': 'Column type',
    }

    # Filter the DataFrame based on the specified chromatographic conditions
    filter_conditions = pd.Series([True] * len(df), index=df.index)
    for arg_name, value in kwargs.items():
        column_name = column_mapping.get(arg_name, arg_name)
        if column_name in df.columns:
            if '|' in str(value):
                filter_conditions &= df[column_name].isin(value.split('|'))
            else:
                filter_conditions &= (df[column_name] == value)

    filtered_df = df[filter_conditions]

    # Define the constant group_by_columns (unique identifiers for each compound)
    group_by_columns = ['Name', 'molecularFormula', 'inChI', 'inChIKey']

    # Aggregate the duplicates based on retention index (I) column
    aggregate_column = 'I'
    aggregated_df = filtered_df.groupby(group_by_columns)[aggregate_column].mean().reset_index()

    return aggregated_df