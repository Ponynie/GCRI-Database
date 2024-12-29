import pandas as pd
from queries import queries
from queries import ri_types, phase_types, temperature_modes

def export_data_counts():
    results = []  # List to store data for CSV export
    s = 0
    for ri_type in ri_types:
        for phase_type in phase_types:
            for temperature_mode in temperature_modes:
                count = len(queries(cleaned=False, RI_Type=ri_type, Phase_Polarity=phase_type, Temperature_Mode=temperature_mode))
                print("-----------------------------------------------------------")
                print(f"RI Type: {ri_type}, Phase Polarity: {phase_type}, Temperature Mode: {temperature_mode}")
                print(count)
                print("-----------------------------------------------------------")
                
                # Append data to results list
                results.append({
                    "RI Type": ri_type,
                    "Phase Polarity": phase_type,
                    "Temperature Mode": temperature_mode,
                    "Count": count
                })
                
                s += count
    print(s)
    
    # Convert results to a DataFrame and export to CSV
    results_df = pd.DataFrame(results)
    results_df.to_csv('counts_summary.csv', index=False)
    print("Counts exported to 'counts_summary.csv'")
    
def export_data_and_average():
    results = []  # List to store data for the summary CSV export
    
    for ri_type in ri_types:
        for phase_type in phase_types:
            for temperature_mode in temperature_modes:
                # Get the filtered data for each combination
                filtered_df = queries(cleaned=False, RI_Type=ri_type, Phase_Polarity=phase_type, Temperature_Mode=temperature_mode)
                
                # Calculate the average 'I' value if the filtered data is not empty
                if not filtered_df.empty:
                    average_value = filtered_df['I'].mean()
                    
                    # Save the filtered data to a separate CSV file
                    filename = f"{ri_type}_{phase_type}_{temperature_mode}.csv".replace(" ", "_").replace("'", "")
                    filtered_df.to_csv(filename, index=False)
                    print(f"Saved CSV for {ri_type}, {phase_type}, {temperature_mode} as '{filename}'")
                else:
                    average_value = None  # or you could set it to 0 or 'N/A' if preferred
                
                # Append data to results list for the summary file
                results.append({
                    "RI Type": ri_type,
                    "Phase Polarity": phase_type,
                    "Temperature Mode": temperature_mode,
                    "Average I": average_value
                })
    
    # Convert results to a DataFrame and export the summary to CSV
    results_df = pd.DataFrame(results)
    results_df.to_csv('average_summary.csv', index=False)
    print("Average values summary exported to 'average_I_summary.csv'")
    


