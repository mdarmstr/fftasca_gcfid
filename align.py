import os
import pandas as pd
import matplotlib.pyplot as plt

def create_vector(filename):
    vector = [0, 0, 0, 0]  # Initialize vector with zeros
    
    # Condition for "Blank"
    if "Blank" in filename:
        return vector
    
    # Condition for "24h" and "72h"
    if "24h" in filename:
        vector[0] = 1
    elif "72h" in filename:
        vector[0] = 2
    
    # Condition for 5th entry from the end
    if filename[-6] == 'N':
        vector[1] = 1
    elif filename[-6] == 'P':
        vector[1] = 2
    elif filename[-6] == 'W':
        vector[1] = 3
    
    # Condition for 4th entry from the end
    if filename[-5] == 'm':
        vector[2] = 1

    vector[3] = int(filename[-4])
    
    return vector

# Define the directory containing the files
directory = 'C:/Users/Michael Armstrong/Documents/gcfid_data/MTBLS2277/FILES/DERIVED_FILES/XY_files_gland'

# Get a list of all files in the directory
files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

# Initialize an empty DataFrame to store the final merged data
merged_df = pd.DataFrame()

# Iterate over the files and read them into DataFrames
for file in files:
    # Construct the full file path
    file_path = os.path.join(directory, file)
    
    # Read the file into a DataFrame
    df = pd.read_csv(file_path, delimiter='\t', header=None, names=['x', file])
    
    # Merge the current DataFrame with the merged_df on the 'x' column
    if merged_df.empty:
        merged_df = df
    else:
        merged_df = pd.merge(merged_df, df, on='x', how='outer')

# Sort the merged DataFrame by the 'x' column
merged_df = merged_df.sort_values(by='x')

# Reset the index of the merged DataFrame
merged_df = merged_df.reset_index(drop=True)

# Create a dictionary to store the vectors for each column
vectors = {}

# Iterate over the column headers (excluding the 'x' column)
for col in merged_df.columns[1:]:
    vectors[col] = create_vector(col)

# Convert the dictionary to a DataFrame for better visualization
vectors_df = pd.DataFrame.from_dict(vectors, orient='index', columns=['Position 1', 'Position 2', 'Position 3', 'Position 4'])

# Display the vectors DataFrame
print(vectors_df)

# Display the merged DataFrame
print(merged_df)

merged_df.iloc[:,1:].to_csv('X_data.csv',index=False, header=False)
vectors_df.to_csv('Y_data.csv',index=False, header=False)