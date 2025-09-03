import pandas as pd
import itertools
import numpy as np

# Define the FT sensor maximum values (converted to base units)
sensor_limits = {
    'Fx': 5000,    # N (was 5 kN)
    'Fy': 5000,    # N (was 5 kN)
    'Fz': 25000,   # N (was 25 kN)
    'Mx': 350,     # Nm (was 0.35 kNm)
    'My': 350,     # Nm (was 0.35 kNm)
    'Mz': 250      # Nm (was 0.25 kNm)
}

# Create combinations with -1, 0, +1 levels
# -1 = negative max, 0 = zero, +1 = positive max
levels = [-1, 0, 1]

# Separate force and torque axes
force_axes = ['Fx', 'Fy', 'Fz']
torque_axes = ['Mx', 'My', 'Mz']

# Generate ONLY force combinations (3^3 = 27 combinations)
# Torques are zero
force_combinations = list(itertools.product(levels, repeat=3))

# Generate ONLY torque combinations (3^3 = 27 combinations)  
# Forces are zero
torque_combinations = list(itertools.product(levels, repeat=3))

# Create the test matrix
test_data = []
test_id = 1

# Add all force combinations
for force_combo in force_combinations:
    test_case = {'Test_ID': test_id}
    
    # Apply force values
    for i, axis in enumerate(force_axes):
        test_case[axis] = force_combo[i] * sensor_limits[axis]
    
    # Set torques to zero
    for axis in torque_axes:
        test_case[axis] = 0
    
    test_data.append(test_case)
    test_id += 1

# Add all torque combinations
for torque_combo in torque_combinations:
    test_case = {'Test_ID': test_id}
    
    # Set forces to zero
    for axis in force_axes:
        test_case[axis] = 0
    
    # Apply torque values  
    for i, axis in enumerate(torque_axes):
        test_case[axis] = torque_combo[i] * sensor_limits[axis]
    
    test_data.append(test_case)
    test_id += 1

# Create DataFrame
df = pd.DataFrame(test_data)

# Add magnitude calculations for analysis
df['Force_Magnitude'] = np.sqrt(df['Fx']**2 + df['Fy']**2 + df['Fz']**2)
df['Torque_Magnitude'] = np.sqrt(df['Mx']**2 + df['My']**2 + df['Mz']**2)

# Reorder columns
column_order = ['Test_ID', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz', 'Force_Magnitude', 'Torque_Magnitude']
df = df[column_order]

# Save to CSV
filename = 'ft_sensor_combinations.csv'
df.to_csv(filename, index=False, float_format='%.0f')

# Display summary
print(f"=== FT SENSOR COMBINATIONS ===")
print(f"Force combinations: 27 (3^3)")
print(f"Torque combinations: 27 (3^3)")
print(f"Total combinations: {len(df)}")
print()

print(f"OUTPUT FILE: {filename}")
print()

print("FIRST 15 COMBINATIONS:")
print(df.head(15))
print()

print("LAST 10 COMBINATIONS:")
print(df.tail(10))
print()

print("SUMMARY:")
print(f"Force values: ±{sensor_limits['Fx']} N, ±{sensor_limits['Fy']} N, ±{sensor_limits['Fz']} N")
print(f"Torque values: ±{sensor_limits['Mx']} Nm, ±{sensor_limits['My']} Nm, ±{sensor_limits['Mz']} Nm")
print()

print("FORCE MAGNITUDE RANGE:")
force_rows = df[(df['Mx'] == 0) & (df['My'] == 0) & (df['Mz'] == 0)]
force_mags = force_rows['Force_Magnitude']
print(f"Min: {force_mags.min():.0f} N, Max: {force_mags.max():.0f} N")
print()

print("TORQUE MAGNITUDE RANGE:")
torque_rows = df[(df['Fx'] == 0) & (df['Fy'] == 0) & (df['Fz'] == 0)]
torque_mags = torque_rows['Torque_Magnitude'] 
print(f"Min: {torque_mags.min():.0f} Nm, Max: {torque_mags.max():.0f} Nm")