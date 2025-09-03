import numpy as np
import pandas as pd
import math
import os
import glob

class TransformationCalculator:
    def __init__(self, radius=1.0):
        self.radius = radius
        self.reference_triangle = self.create_reference_triangle()
        
    def create_reference_triangle(self):
        """Create a perfect equilateral triangle with given radius in XY plane"""
        # Three points of equilateral triangle, centered at origin
        angle_offset = 2 * math.pi / 3  # 120 degrees
        points = []
        
        for i in range(3):
            angle = i * angle_offset
            x = self.radius * math.cos(angle)
            y = self.radius * math.sin(angle)
            z = 0.0
            points.append([x, y, z])
            
        return np.array(points)
    
    def calculate_center_and_normal(self, points):
        """Calculate center and normal vector from 3 points"""
        points = np.array(points)
        
        # Center is average of points
        center = np.mean(points, axis=0)
        
        # Calculate normal using cross product
        v1 = points[1] - points[0]
        v2 = points[2] - points[0]
        normal = np.cross(v1, v2)
        
        # Normalize
        normal_length = np.linalg.norm(normal)
        if normal_length > 0:
            normal = normal / normal_length
        else:
            normal = np.array([0, 0, 1])  # Default to Z-up
            
        return center, normal
    
    def create_transformation_matrix(self, center, normal):
        """Create 4x4 transformation matrix from center and normal"""
        # Create rotation matrix to align Z-axis with normal
        z_axis = np.array([0, 0, 1])
        
        if np.allclose(normal, z_axis):
            # Already aligned
            rotation = np.eye(3)
        elif np.allclose(normal, -z_axis):
            # Opposite direction, rotate 180 degrees around X
            rotation = np.array([[ 1,  0,  0],
                               [ 0, -1,  0],
                               [ 0,  0, -1]])
        else:
            # General case: rotate z_axis to normal
            axis = np.cross(z_axis, normal)
            axis = axis / np.linalg.norm(axis)
            
            cos_angle = np.dot(z_axis, normal)
            sin_angle = np.linalg.norm(np.cross(z_axis, normal))
            
            # Rodrigues' rotation formula
            K = np.array([[     0, -axis[2],  axis[1]],
                         [ axis[2],      0, -axis[0]],
                         [-axis[1],  axis[0],     0]])
            
            rotation = np.eye(3) + sin_angle * K + (1 - cos_angle) * np.dot(K, K)
        
        # Create 4x4 transformation matrix
        transform = np.eye(4)
        transform[0:3, 0:3] = rotation
        transform[0:3, 3] = center
        
        return transform
    
    def calculate_relative_transformation(self, points_a1, points_a2):
        """Calculate transformation from A1 coordinate system to A2"""
        # Calculate centers and normals
        center_a1, normal_a1 = self.calculate_center_and_normal(points_a1)
        center_a2, normal_a2 = self.calculate_center_and_normal(points_a2)
        
        # Create transformation matrices
        transform_a1 = self.create_transformation_matrix(center_a1, normal_a1)
        transform_a2 = self.create_transformation_matrix(center_a2, normal_a2)
        
        # Relative transformation: A1^-1 * A2
        transform_a1_inv = np.linalg.inv(transform_a1)
        relative_transform = np.dot(transform_a2, transform_a1_inv)
        
        return relative_transform
    
    def get_column_names(self, df):
        """Extract the 3 node column names from DataFrame (generic for any sensor)"""
        columns = df.columns.tolist()
        
        # Group columns by node name
        nodes = {}
        for col in columns:
            if col.endswith('_X') or col.endswith('_Y') or col.endswith('_Z'):
                node_name = col[:-2]  # Remove _X, _Y, or _Z
                if node_name not in nodes:
                    nodes[node_name] = []
                nodes[node_name].append(col)
        
        # Should have exactly 3 nodes
        if len(nodes) != 3:
            raise ValueError(f"Expected 3 nodes, found {len(nodes)}: {list(nodes.keys())}")
        
        # Sort to ensure consistent order
        node_names = sorted(nodes.keys())
        return node_names
    
    def extract_displacements(self, df, row_index, node_names):
        """Extract XYZ displacements for given row and nodes"""
        displacements = []
        for node in node_names:
            x = df.iloc[row_index][f'{node}_X']
            y = df.iloc[row_index][f'{node}_Y']
            z = df.iloc[row_index][f'{node}_Z']
            displacements.append([x, y, z])
        return np.array(displacements)
    
    def process_sensor_pair(self, file1, file2, sensor_name):
        """Process a single sensor pair (e.g., A1 + A2)"""
        print(f"\nProcessing sensor {sensor_name}...")
        
        # Read CSV files
        df1 = pd.read_csv(file1)
        df2 = pd.read_csv(file2)
        
        print(f"  {sensor_name}1: {len(df1)} rows")
        print(f"  {sensor_name}2: {len(df2)} rows")
        
        # Get node names for each sensor
        nodes1 = self.get_column_names(df1)
        nodes2 = self.get_column_names(df2)
        
        print(f"  {sensor_name}1 nodes: {nodes1}")
        print(f"  {sensor_name}2 nodes: {nodes2}")
        
        # Ensure same number of rows
        min_rows = min(len(df1), len(df2))
        
        transformations = []
        
        for i in range(min_rows):
            # Extract displacements
            displacements1 = self.extract_displacements(df1, i, nodes1)
            displacements2 = self.extract_displacements(df2, i, nodes2)
            
            # Calculate actual positions
            actual1 = self.reference_triangle + displacements1
            actual2 = self.reference_triangle + displacements2
            
            # Calculate relative transformation
            rel_transform = self.calculate_relative_transformation(actual1, actual2)
            transformations.append(rel_transform)
        
        print(f"  Generated {len(transformations)} transformation matrices")
        return transformations
    
    def write_transformations(self, transformations, filename, sensor_name):
        """Write transformation matrices to file in C++ readable format"""
        with open(filename, 'w') as f:
            f.write(f"# Transformation matrices for sensor {sensor_name} ({sensor_name}2 relative to {sensor_name}1)\n")
            f.write(f"# Number of matrices: {len(transformations)}\n")
            f.write(f"MATRICES {len(transformations)}\n")
            
            for i, matrix in enumerate(transformations):
                f.write(f"MATRIX {i}\n")
                for row in matrix:
                    f.write(" ".join([f"{val:.12f}" for val in row]) + "\n")
                f.write("\n")
        
        print(f"  Written to {filename}")
    
    def process_all_sensors(self, data_folder="data"):
        """Process all sensor pairs in the data folder"""
        if not os.path.exists(data_folder):
            print(f"Error: Data folder '{data_folder}' not found!")
            return
        
        # Find all sensor pairs
        sensors = ['A', 'B', 'C']
        
        for sensor in sensors:
            file1 = os.path.join(data_folder, f"{sensor}1Sample.csv")
            file2 = os.path.join(data_folder, f"{sensor}2Sample.csv")
            
            if os.path.exists(file1) and os.path.exists(file2):
                try:
                    # Process this sensor pair
                    transformations = self.process_sensor_pair(file1, file2, sensor)
                    
                    # Write output file
                    output_file = f"transformations_{sensor}.txt"
                    self.write_transformations(transformations, output_file, sensor)
                    
                except Exception as e:
                    print(f"  Error processing sensor {sensor}: {e}")
            else:
                missing = []
                if not os.path.exists(file1):
                    missing.append(f"{sensor}1Sample.csv")
                if not os.path.exists(file2):
                    missing.append(f"{sensor}2Sample.csv")
                print(f"Skipping sensor {sensor}: Missing files {missing}")
        
        print(f"\nProcessing complete!")
        print("Generated files:")
        for sensor in sensors:
            output_file = f"transformations_{sensor}.txt"
            if os.path.exists(output_file):
                print(f"  - {output_file}")

if __name__ == "__main__":
    calculator = TransformationCalculator(radius=1.0)
    
    # Process all sensors in the data folder
    calculator.process_all_sensors("data")